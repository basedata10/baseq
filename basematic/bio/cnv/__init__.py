import click
from basematic.mgt.resource import ResourceModel

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    click.echo("Welcome to Basematic-CNV")

#Run the CNV pipeline for Lists of Samples
@cli.command(short_help="Run")
@click.option('--sample_file', '-m', default='', help="Tab seprated file: name, fq1, fq2")
@click.option('--species', '-s', default='', help="Species hg19 or mm10")
def run_pipeline(sample_file):
    samples = []
    with open(sample_file, 'r') as infile:
        samples = infile.readlines()
    pass

#Run BOWTIE2 Alignment
@cli.command(short_help="Align the reads to Genome using bowtie2")
@click.argument("fastq")
@click.option('--genome', '-g', default='hg19', help="genome short name, hg19/hg38/mm10")
def run_alignment(fastq, genome):
    from .align import bowtie2_sort_alignment
    bowtie2_sort_alignment(fastq, genome)

#Run Bin Counting
@cli.command(short_help="Couting reads in the bam file according to dynamic bins")
@click.option('--genome', '-g', default='', help="genome")
@click.option('--bamfile', '-i', default='', help="bam file path")
def run_bincounting(genome, bamfile):
    from .bin_count import bin_counting
    bin_counting(genome, bamfile)

#Run Normalize
@cli.command(short_help="Couting reads in the bam file according to dynamic bins")
@click.option('--genome', '-g', default='', help="genome")
@click.option('--count_file', '-i', default='', help="bin count file")
def run_normalize(genome, count_file):
    from .normalize import Normalize_GC
    Normalize_GC(genome, count_file)