from baseq.rna import cli
import click, os, sys

@cli.command(short_help = "Show the status of jobs")
@click.option('--namefile', '-g', default='', help='GTF file')
@click.option('--table', '-t', default='', help='GTF file')
@click.option('--out', '-o', default='', help='Outfile Name')
def table_ensg_to_name(namefile, table, out):
    from .replace import ensg_to_genename
    ensg_to_genename(namefile, table, out)

@cli.command(short_help = "Mean Expression of Table...")
@click.option('--table', '-t', default='', help='GTF file')
@click.option('--out', '-o', default='', help='Outfile Name')
def table_mean_expression(table, out):
    pass