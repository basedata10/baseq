import click
from baseq.snv import cli

@cli.command(short_help="Check the enrichment quality, input: bam, interval and outpath")
@click.argument("bampath")
@click.argument("interval")
@click.argument("outpath")
def QC_enrich(bampath, interval, outpath):
    from .enrich import enrich_quality
    QC_enrich(bampath, interval, outpath)