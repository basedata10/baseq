import click, os, sys
from baseq.rna import cli
from baseq.fastq.sample_file import check_sample_files
from concurrent.futures import ThreadPoolExecutor
import multiprocessing as mp

@cli.command(short_help="inDrop/Drop-Seq/10X")
@click.option('--genome', '-g', help="human/mouse/mixed")
@click.option('--protocol', '-p', type=click.Choice(['10X', 'indrop', 'dropseq']))
@click.option('--cells', default=5000, help='Expected number of cells')
@click.option('--minreads', default=10000, help='Minimum reads for a barcode')
@click.option('--name', '-n', default='sample', help="sample name")
@click.option('--fq1', '-1', default='', help='fastq1 path')
@click.option('--fq2', '-2', default='', help='fastq2 path (optional)')
@click.option('--top_million_reads', default=10, help='Use the top X M reads')
@click.option('--parallel', default='4', help='Process N splitted at the same time...')
@click.option('--dir', '-d', default='./', help='Folder (./)')
@click.option('--step', default='all', help='all/star/')

def drops_pipe(genome, protocol, cells, minreads, name, fq1, fq2, dir, top_million_reads, step, parallel):
    print('Start Processing ...')
    from baseq.rna.dropseq.barcode_count import count_barcodes
    from baseq.rna.dropseq.barcode_stats import barcode_correct_filter
    from baseq.rna.dropseq.barcode_split import barcode_split
    dir = os.path.abspath(os.path.join(dir, name))
    bc_counts = os.path.join(dir, "barcode_count_{}.csv".format(name))
    bc_stats = os.path.join(dir, "barcode_stats_{}.csv".format(name))
    bc_splits_dir = os.path.join(dir, "barcode_splits")
    align_dir = os.path.join(dir, "star_align")
    tagging_dir = os.path.join(dir, "read_tagging")
    tpm_table = os.path.join(dir, "TPM.{}.txt".format(name))

    from itertools import product
    barcode_prefix = [x[0] + x[1] for x in list(product('ATCG', repeat=2))]

    dirs = [dir, align_dir, tagging_dir, bc_splits_dir]
    for dir in dirs:
        if not os.path.exists(dir):
            os.mkdir(dir)

    #count barcode
    if step in ["all", "count"]:
        print("[info] Counting the barcodes ...")
        count_barcodes(fq1, bc_counts, protocol, 30, int(top_million_reads))

    #aggregate
    if step in ["all", "stats"]:
        print("[info] Aggregating the barcodes errors ...")
        barcode_correct_filter(protocol, bc_counts, max_cell=cells, min_reads=minreads, output=bc_stats)

    #barcode split
    if step in ["all", "split"]:
        print("[info] Split the barcode ...")
        barcode_split(name, protocol, bc_stats, fq1, fq2, minreads, cells, bc_splits_dir, int(top_million_reads))

    #run alignment
    if step in ["all", "star"]:
        from baseq.rna.dropseq.run_star import run_star_multiple
        run_star_multiple(bc_splits_dir, align_dir, name, genome, parallel=4)

    #run reads tagging
    if step in ["all", "tagging"]:
        from baseq.rna.dropseq.tag_gene import tagging_reads
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
        from baseq.rna.dropseq.aggregate import read_barcode_gene_file, write_to_table
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

@cli.command(short_help="Barcode Counting")
@click.option('--name', '-n', default='sample', help="sample name")
@click.option('--protocol', '-p', type=click.Choice(['10X', 'indrop', 'dropseq']))
@click.option('--fq1', '-1', default='', help='fastq1, reads with barcode/UMI')
@click.option('--dir', '-d', default='./', help='Folder of output (./)')
def drops_barcode_counting(name, protocol, fq1, dir):
    print('Start Processing inDrop Results')
    samples = check_sample_files("", name, fq1)
    if samples == []:
        sys.exit("[error] No valid sample, Exit.")
    from baseq.rna.dropseq.barcode_count import count_barcodes
    outpath = os.path.join(dir, "barcode_count.{}.csv".format(name))
    count_barcodes(fq1, outpath, protocol, 20)

@cli.command(short_help="Barcode Stats")
@click.option('--name', '-n', default='sample', help="sample name")
@click.option('--protocol', '-p', type=click.Choice(['10X', 'indrop', 'dropseq']))
@click.option('--bcfile', '-b', default='', help='Barcode count file')
@click.option('--minreads', default=2000, help='Minimum reads for a barcode')
@click.option('--dir', '-d', default='./', help='Folder of output (./)')
def drops_barcode_stats(name, protocol, bcfile, minreads, dir):
    print('Start Processing inDrop Results')
    from baseq.rna.dropseq.barcode_stats import barcode_correct_filter
    outpath = os.path.join(dir, "barcode_stats.{}.csv".format(name))
    barcode_correct_filter(protocol = protocol, barcode_count=bcfile, min_reads=int(minreads), output=outpath)

@cli.command(short_help="Split Barcode")
@click.option('--name', '-n', default='sample', help="sample name")
@click.option('--protocol', '-p', type=click.Choice(['10X', 'indrop', 'dropseq']))
@click.option('--bcstats', '-b', help='Files generated by barcode stats step')
@click.option('--fq1', '-1', default='', help='fastq1, reads with barcode/UMI')
@click.option('--fq2', '-2', default='', help='fastq2')
@click.option('--minreads', default=2000, help='Minimum reads for a barcode')
@click.option('--maxcell', default=10000, help='Max cell number')
@click.option('--dir', '-d', default='./', help='Folder of output (./)')
def drops_barcode_split(name, protocol, bcstats, fq1, fq2, minreads, maxcell, dir):
    from baseq.rna.dropseq.barcode_split import barcode_split
    print('[info] Split the cell barcodes ...')
    barcode_split(name, protocol, bcstats, fq1, fq2, minreads, maxcell, dir)

@cli.command(short_help="Split Barcode")
@click.option('--sample', '-n', default='sample', help="sample name")
@click.option('--genome', '-g', default='hg38', help="hg38/hg19/")
@click.option('--bcdir', '-d', default='./', help='Folder of barcodes')
@click.option('--aligndir', '-o', default='./', help='Folder of output (./)')
def drops_star_align_script(sample, bcdir, genome, aligndir):
    from baseq.rna.dropseq.run_star import genrate_star_script
    print('[info] Split the cell barcodes ...')
    genrate_star_script(bcdir, genome, sample, aligndir, "qsub")

@cli.command(short_help="Split Barcode")
@click.option('--genome', '-g', default='hg38', help="Genome:hg38...")
@click.option('--bam', '-b', help='Path to the bam file')
@click.option('--out', '-o', default='./', help='Output file path')
def drops_reads_tagging(bam, genome, out):
    from baseq.rna.dropseq.tag_gene import tagging_reads
    print('[info] Tagging the reads ...')
    tagging_reads(genome, bam, out)

@cli.command(short_help="Making The Expression Table")
@click.option('--bcgenefile', '-b', help='Path to the barcode gene file...')
@click.option('--out', '-o', default='./', help='Output file path')
def drops_barcode_gene_table(bcgenefile, out):
    from baseq.rna.dropseq.aggregate import read_barcode_gene_file
    print('[info] Tagging the reads ...')
    read_barcode_gene_file(bcgenefile)