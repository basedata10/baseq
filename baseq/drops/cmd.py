import click, sys
from baseq.fastq.sample_file import check_sample_files
import multiprocessing as mp
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """Run DropSeq Data Analysis
    """
    pass

from baseq.drops.apa.cmd import *
from baseq.drops.cellranger.cmd import *

@cli.command(short_help="Main cmd for inDrop/Drop-Seq/10X")
@click.option('--config', default="", help="config path")
@click.option('--genome', '-g', help="human/mouse/mixed")
@click.option('--protocol', '-p', type=click.Choice(['10X', 'indrop', 'dropseq']))
@click.option('--cells', default=5000, help='Expected number of cells')
@click.option('--minreads', default=10000, help='Minimum reads for a barcode')
@click.option('--name', '-n', default='sample', help="sample name")
@click.option('--fq1', '-1', default='', help='fastq1 path')
@click.option('--fq2', '-2', default='', help='fastq2 path (optional)')
@click.option('--top_million_reads', default=1000, help='Use the top X M reads (1000M for default)')
@click.option('--parallel', default='4', help='Process N splitted at the same time...')
@click.option('--dir', '-d', default='./', help='Folder (./)')
@click.option('--step', default='all', help='all/star/')

def run_pipe(config, genome, protocol, cells, minreads, name, fq1, fq2, dir, top_million_reads, step, parallel):
    """
    Run the Total Pipeline...

    #. Read the barcode counts files;
    #. Correct the barcode with 1bp mismatch;
    #. Stats the mismatch barcode reads and sequences;
    #. Determine wheather mutate on the last base (show A/T/C/G with similar ratio at the last base);
    #. Filter by whitelist;
    #. Filter by read counts (>=min_reads);
    #. Print the number of barcode and reads retained after each steps.

    Pipeline:
    ::
        baseq-Drop runpipe -1 1.fq.gz -2 2.fq.gz -n test

    Return:
        XXXX
    """

    print('Start Processing ...')
    from baseq.drops.barcode.count import count_barcodes
    from baseq.drops.barcode.stats import valid_barcode
    from baseq.drops.barcode.split import split_16
    dir = os.path.abspath(os.path.join(dir, name))
    bc_counts = os.path.join(dir, "barcode_count_{}.csv".format(name))
    bc_stats = os.path.join(dir, "barcode_stats_{}.csv".format(name))
    bc_splits_dir = os.path.join(dir, "barcode_splits")
    align_dir = os.path.join(dir, "star_align")
    tagging_dir = os.path.join(dir, "read_tagging")
    tpm_table = os.path.join(dir, "UMIs.{}.txt".format(name))

    from itertools import product
    barcode_prefix = [x[0] + x[1] for x in list(product('ATCG', repeat=2))]

    dirs = [dir, align_dir, tagging_dir, bc_splits_dir]
    for dir in dirs:
        if not os.path.exists(dir):
            os.mkdir(dir)

    if config:
        if not os.path.exists(config):
            sys.exit("[error] Config file not exists")
        os.environ["BASEQCFG"] = config

    #count barcode
    if step in ["all", "count"]:
        print("[info] Counting the barcodes ...")
        count_barcodes(fq1, bc_counts, protocol, 30, int(top_million_reads))

    #aggregate
    if step in ["all", "stats"]:
        print("[info] Aggregating the barcodes errors ...")
        valid_barcode(protocol, bc_counts, max_cell=cells, min_reads=minreads, output=bc_stats)

    #barcode split
    if step in ["all", "split"]:
        print("[info] Split the barcode ...")
        split_16(name, protocol, bc_stats, fq1, fq2, bc_splits_dir, int(top_million_reads))

    #run alignment
    if step in ["all", "star"]:
        from baseq.drops.run_star import run_star_multiple
        run_star_multiple(bc_splits_dir, align_dir, name, genome, parallel=2)

    #run reads tagging
    if step in ["all", "tagging"]:
        from baseq.drops.tag_gene import tagging_reads
        print('[info] Tagging the reads to genes...')
        pool = mp.Pool(processes=int(parallel))
        for bc in barcode_prefix:
            bamfile = os.path.join(align_dir, "{}/{}.sort.bam".format(bc, bc))
            outfile = os.path.join(tagging_dir, "tagging.{}.txt".format(bc))
            pool.apply_async(tagging_reads, (genome, bamfile, outfile,))
        pool.close()
        pool.join()

    #run Table aggragation
    if step in ["all", "table"]:
        from baseq.drops.aggregate import read_barcode_gene_file, write_to_table
        pool = mp.Pool(processes=int(parallel))
        result = []
        for bc in barcode_prefix:
            filepath = os.path.join(tagging_dir, "tagging.{}.txt".format(bc))
            result.append(pool.apply_async(read_barcode_gene_file, (filepath,)))
        pool.close()
        pool.join()
        from itertools import chain
        barcodes_all = [x.get()[0] for x in result]
        barcodes_lists = list(chain(*barcodes_all))
        exp = {}
        UMIs_all = [x.get()[1] for x in result]
        for UMI in UMIs_all:
            for gene in UMI:
                if gene in exp:
                    exp[gene].update(UMI[gene])
                else:
                    exp[gene] = UMI[gene]
        write_to_table(barcodes_lists, exp, tpm_table)

@cli.command(short_help="Barcode counting")
@click.option('--name', '-n', default='sample', help="sample name")
@click.option('--protocol', '-p', type=click.Choice(['10X', 'indrop', 'dropseq']))
@click.option('--fq1', '-1', default='', help='fastq1, reads with barcode/UMI')
@click.option('--dir', '-d', default='./', help='Folder of output (./)')
def barcode_counting(name, protocol, fq1, dir):
    print('Start Processing inDrop Results')
    samples = check_sample_files("", name, fq1)
    if samples == []:
        sys.exit("[error] No valid sample, Exit.")
    from baseq.drops.barcode.count import count_barcodes
    outpath = os.path.join(dir, "barcode_count.{}.csv".format(name))
    count_barcodes(fq1, outpath, protocol, 20)

@cli.command(short_help="Barcode correlation, filter low abundance")
@click.option('--name', '-n', default='sample', help="sample name")
@click.option('--protocol', '-p', type=click.Choice(['10X', 'indrop', 'dropseq']))
@click.option('--bcfile', '-b', default='', help='Barcode count file')
@click.option('--minreads', default=2000, help='Minimum reads for a barcode')
@click.option('--dir', '-d', default='./', help='Folder of output (./)')
def barcode_filter(name, protocol, bcfile, minreads, dir):
    print('Start Processing inDrop Results')
    from baseq.drops.barcode.stats import valid_barcode
    outpath = os.path.join(dir, "barcode_stats.{}.csv".format(name))
    valid_barcode(protocol = protocol, barcode_count=bcfile, min_reads=int(minreads), output=outpath)

@cli.command(short_help="Split Barcode")
@click.option('--name', '-n', default='sample', help="sample name")
@click.option('--protocol', '-p', type=click.Choice(['10X', 'indrop', 'dropseq']))
@click.option('--bcstats', '-b', help='Files generated by barcode stats step')
@click.option('--fq1', '-1', default='', help='fastq1, reads with barcode/UMI')
@click.option('--fq2', '-2', default='', help='fastq2')
@click.option('--minreads', default=2000, help='Minimum reads for a barcode')
@click.option('--maxcell', default=10000, help='Max cell number')
@click.option('--dir', '-d', default='./', help='Folder of output (./)')
def barcode_split(name, protocol, bcstats, fq1, fq2, minreads, maxcell, dir):
    from baseq.drops.barcode.split import split_16
    print('[info] Split the cell barcodes ...')
    split_16(name, protocol, bcstats, fq1, fq2, minreads, maxcell, dir)


@cli.command(short_help="Split Barcode")
@click.option('--name', '-n', default='sample', help="sample name")
@click.option('--protocol', '-p', type=click.Choice(['10X', 'indrop', 'dropseq']))
@click.option('--bcstats', '-b', help='Files generated by barcode stats step')
@click.option('--fq1', '-1', default='', help='fastq1, reads with barcode/UMI')
@click.option('--fq2', '-2', default='', help='fastq2')
@click.option('--minreads', default=2000, help='Minimum reads for a barcode')
@click.option('--maxcell', default=10000, help='Max cell number')
@click.option('--dir', '-d', default='./', help='Folder of output (./)')
def barcode_split_fast(name, protocol, bcstats, fq1, fq2, minreads, maxcell, dir):
    from baseq.drops.barcode.split_fast import barcode_splits_all
    print('[info] Split the cell barcodes ...')
    barcode_splits_all(name, protocol, bcstats, fq1, fq2, minreads, maxcell, dir)

@cli.command(short_help="Genome position to gene name")
@click.option('--genome', '-g', default='hg38', help="Genome:hg38...")
@click.option('--bam', '-b', help='Path to the bam file')
@click.option('--out', '-o', default='./', help='Output file path')
def reads_tagging(bam, genome, out):
    from baseq.drops.tag_gene import tagging_reads
    print('[info] Tagging the reads ...')
    tagging_reads(genome, bam, out)