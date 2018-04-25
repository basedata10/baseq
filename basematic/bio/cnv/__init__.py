import click, sys, os
from basematic.mgt.resource import ResourceModel

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    click.echo("Welcome to Basematic-CNV")

#Run the CNV pipeline for Lists of Samples
@cli.command(short_help="Run the whole pipeline for a sample list")
@click.option('--sample_file', '-m', default='', help="Tab seprated file: name, fq1, fq2")
@click.option('--species', '-s', default='', help="Species hg19 or mm10")
@click.option('--outdir', '-d', default='./CNV_process', help="The scripts and output path")
def run_pipeline(sample_file):
    samples = []
    with open(sample_file, 'r') as infile:
        samples = infile.readlines()
    pass

#Run BOWTIE2 Alignment
@cli.command(short_help="Alignment with bowtie2, get sorted bam")
@click.argument("fastq")
@click.option('--genome', '-g', default='hg19', help="genome version name, like: hg19")
def run_alignment(fastq, genome):
    from .align import bowtie2_sort_alignment
    bowtie2_sort_alignment(fastq, genome)

#Run Bin Counting
@cli.command(short_help="Couting reads in the bam file according to dynamic bins")
@click.option('--genome', '-g', default='', help="genome")
@click.option('--bamfile', '-i', default='', help="bam file path")
@click.option('--out', '-o', default='./sample.bincounts.txt', help="bin counts file path")
def run_bincounting(genome, bamfile, out):
    from .bincounting import bin_counting
    bin_counting(genome, bamfile, out)

#Run Normalize
@cli.command(short_help="Normalize the reads to Depth and GC content")
@click.option('--genome', '-g', default='', help="genome")
@click.option('--bincount', '-i', default='', help="bin count file")
@click.option('--out', '-o', default="./sample.norm_GC.txt", help="normalized bin counts")
def run_normalize(genome, bincount, out):
    from .normalize import Normalize_GC
    Normalize_GC(genome, bincount, out)

#Run QC
@cli.command(short_help="Normalize the reads to Depth and GC content")
@click.option('--genome', '-g', default='', help="genome")
@click.option('--bincount', '-i', default='', help="bin count file")
@click.option('--out', '-o', default="./sample.norm_GC.txt", help="normalized bin counts")
def run_quality_control(genome, bincount, out):
    from .normalize import Normalize_GC
    Normalize_GC(genome, bincount, out)

#Run CBS
@cli.command(short_help="Run CBS ...")
@click.option('--count_file', '-i', default='', help="Normalized bin count file")
@click.option('--out', '-o', default='./', help="CBS file...")
def run_CBS(count_file, out):
    from .segmentation import bin_segmentation
    bin_segmentation(count_file, out)

#Run Plot Genomes...
@cli.command(short_help="Run CBS ...")
@click.option('--count_file', '-i', default='', help="Normalized bin count file")
@click.option('--cbs_file', '-c', default='', help="CBS File, Genrated from run_CBS command")
@click.option('--out', '-o', default='./', help="CBS file...")
def run_plotgenome(count_file, cbs_file, out):
    from .plot_genome import PlotGenomes
    PlotGenomes(count_file, cbs_file, out)

#Prepare Web Datas
@cli.command(short_help="Generating the datas for web view: basematic.io/viewcnv")
@click.option('--path', '-p', default='', help="The path to the Process Folder")
@click.option('--outfile', '-o', default='CNV_result.json', help="The path to the packed file")
def web_data(path, outfile):
    from .results import build_result
    build_result(path, outfile)