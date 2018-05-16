import click, re, os
from jinja2 import Template
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass

from .cmd import *
from .bamtype import BAMTYPE