from baseq.drops.cmd import cli
from baseq.mgt.config import get_config
import click, os

@cli.command(short_help="Get Barcodes for each celltype")
@click.option('--markers', '-m', default='', help="Genome:hg38...")
@click.option('--cluster', '-c', default='', help="Genome:hg38...")
@click.option('--express', '-e', default='', help="Genome:hg38...")
def get_celltype_barcodes(markers, cluster, express):
    pass