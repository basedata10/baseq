import click, sys, os
from basematic.mgt.resource import ResourceModel

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    click.echo("Welcome to Basematic-RNA")

#Run the CNV pipeline for Lists of Samples
@cli.command(short_help="RNA-Seq")
@click.option('--sample_file', '-m', default='', help="Tab seprated file: name, fq1, fq2")
@click.option('--species', '-s', default='', help="Species hg19 or mm10")
@click.option('--outdir', '-d', default='./CNV_process', help="The scripts and output path")
def run_TPM(sample_file):
    samples = []
    print("SALMON!!!!")
    with open(sample_file, 'r') as infile:
        samples = infile.readlines()
    pass