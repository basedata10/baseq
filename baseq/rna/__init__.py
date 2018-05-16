import click
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    click.echo("Welcome to baseq-RNA")

from .cmd import *