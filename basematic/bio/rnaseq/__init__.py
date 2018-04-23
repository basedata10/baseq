import click, sys
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

from basematic.bio.fastq.Files import check_infiles
message = """
Basematic-RNA : Make Sense of RNA-Seq datas fastly.
"""

@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    print(message)

@cli.command(short_help="inDrop/Drop-Seq/10X")
@click.option('--genome', help="human/mouse/mixed")
@click.option('--name', '-n', default='sample', help="sample name")
@click.option('--fq1', '-1', default='', help='fastq1 path')
@click.option('--fq2', '-2', default='', help='fastq2 path (optional)')
@click.option('--dir', '-d', default='./', help='Folder of output (./)')
def run_drops(name, fq1, fq2, genome, dir):
    print('Start Processing inDrop Results')
    samples = check_infiles("", name, fq1, fq2)
    if samples == []:
        sys.exit("[error] No valid sample, Exit.")
    from basematic.bio.rnaseq.barcode import getBarcode
    from basematic.bio.rnaseq.barcode_stats import barcode_aggregate
    getBarcode(fq1, "./barcode_count.csv", "10X", 20)
    barcode_aggregate(barcode_count="./barcode_count.csv")

@cli.command(short_help="Run")
@click.option('-species', help="human mouse zebrafish")
def run_bulk(species):
    click.echo('Start Using RNA')