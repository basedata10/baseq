import click, sys

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass

from baseq.snv.cmd import *
from baseq.snv.qc.cmd import *
from baseq.snv.vcf.cmd import *