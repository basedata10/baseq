import click, os, sys
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass

#Print the Doc...
@cli.command(short_help="Print the document...")
def doc():
    from .docs import print_doc
    print_doc()

#Run the CNV pipeline for Lists of Samples
@cli.command(short_help="Run the whole pipeline for a sample list")
@click.option('--sample_file', '-m', default='', help="Tab seprated file: name, fq1, fq2")
@click.option('--fq1', '-1', default='', help="Fastq 1")
@click.option('--fq2', '-2', default='', help="Fastq 2")
@click.option('--genome', '-s', default='', help="Species hg19 or mm10")
@click.option('--outdir', '-d', default='./CNV_process', help="The scripts and output path")

def run_pipeline(sample_file, fq1, fq2, genome):
    samples = []
    with open(sample_file, 'r') as infile:
        samples = infile.readlines()
    pass

#Run BOWTIE2 Alignment
@cli.command(short_help="Alignment with bowtie2, Sorted and Index")
@click.option('--fq1', '-1', default='', help="Fastq 1")
@click.option('--fq2', '-2', default='', help="Fastq 2")
@click.option('--bamfile', '-o', default='Bowtie2.sort.bam', help="Output bamfile name (Bowtie2.sort.bam)")
@click.option('--reads', '-r', default='5000000', help="Numbers of reads")
@click.option('--thread', '-t', default=8, help="Numbers of Thread")
@click.option('--genome', '-g', default='hg19', help="Genome version : hg19/mm38")
def align(fq1, fq2, bamfile, genome, reads, thread):
    from baseq.align.bowtie2 import bowtie2_sort
    from baseq.mgt.config import get_config
    genome = get_config("CNV_ref_"+genome, "bowtie2_index")
    bowtie2_sort(fq1, fq2, bamfile, genome, reads=reads, thread=thread)

#Bin Counting
@cli.command(short_help="Couting reads in the bam file according to dynamic bins")
@click.option('--genome', '-g', default='hg19', help="Genome ID, like: hg19")
@click.option('--bamfile', '-i', default='', help="Bamfile path")
@click.option('--out', '-o', default='./sample.bincounts.txt', help="bin counts file path")
def bincount(genome, bamfile, out):
    from baseq.cnv.bincount import bin_counting
    bin_counting(genome, bamfile, out)

#Run Normalize
@cli.command(short_help="Normalize the reads to Depth and GC content")
@click.option('--genome', '-g', default='', help="genome")
@click.option('--bincount', '-i', default='', help="bin count file")
@click.option('--out', '-o', default="./sample.norm_GC.txt", help="normalized bin counts")
def normalize(genome, bincount, out):
    from baseq.cnv.normalize import Normalize_GC_py
    Normalize_GC_py(genome, bincount, out)

#Run QC
@cli.command(short_help="Normalize the reads to Depth and GC content")
@click.option('--genome', '-g', default='', help="genome")
@click.option('--bincount', '-i', default='', help="bin count file")
@click.option('--out', '-o', default="./sample.norm_GC.txt", help="normalized bin counts")
def run_quality_control(genome, bincount, out):
    from .normalize import Normalize_GC_py
    Normalize_GC_py(genome, bincount, out)

#Run CBS
@cli.command(short_help="Run CBS ...")
@click.option('--count_file', '-i', default='', help="Normalized bin count file")
@click.option('--out', '-o', default='./', help="CBS file...")
def CBS(count_file, out):
    from .segment import bin_segmentation
    bin_segmentation(count_file, out)

#Run Plot Genomes...
@cli.command(short_help="Run CBS ...")
@click.option('--counts', '-i', default='', help="Normalized bin count file")
@click.option('--countlists', '-l', default='', help="Normalized bin count file")
@click.option('--cbs', '-c', default='', help="CBS File, Genrated from run_CBS command")
@click.option('--out', '-o', default='./', help="CBS file...")
def plotgenome(counts, countlists, cbs, out):
    from baseq.cnv.plots.genome import plot_genome_py, plot_genome_multiple
    if not countlists:
        plot_genome_py(counts, cbs, out)
    else:
        plot_genome_multiple(countlists, cbs, out)

#Prepare Web Datas
@cli.command(short_help="Generating the datas for web view: baseq.io/viewcnv")
@click.option('--path', '-p', default='', help="The path to the Process Folder")
@click.option('--outfile', '-o', default='CNV_result.json', help="The path to the packed file")
def web_data(path, outfile):
    from .results import build_result
    build_result(path, outfile)