import click, re, os
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
from baseq.fastq.sample_file import check_sample_files

@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass

@cli.command(short_help="Install Salmon and related references")
@click.option('--path', '-d', default='./', help='File Path')
@click.option('--species', '-s', default='human', help='Species: human, mouse, zebrafish')
def stats_match_lenght(path, species):
    print('RNA-Seq')

@cli.command(short_help="split reads not in region...")
@click.option('--path', '-d', default='./', help='File Path')
@click.option('--species', '-s', default='human', help='Species: human, mouse, zebrafish')
def split_reads_not_in_region(path, species):
    print('RNA-Seq')

@cli.command(short_help="check_bam")
@click.option('--path', '-b', default='./', help='File Path')
def check_bam(path):
    from .stats import check_bam
    check_bam(path)