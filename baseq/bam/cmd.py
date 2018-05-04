import click, os, sys
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.group(context_settings=CONTEXT_SETTINGS)

def cli():
    pass

@cli.command(short_help="Bam Depth...")
@click.option('--thread', '-t', default='1', help='')
@click.option('--mem', '-m', default='2G', help='Memory (2G)')
def QSUB(script_path, thread, mem):
    print(script_path, thread, mem)

