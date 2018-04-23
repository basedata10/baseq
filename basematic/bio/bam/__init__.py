import click
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

from basematic.bio.fastq.Files import check_infiles
message = """
Basematic-BAM: Stats on Bam File...
"""

@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    print(message)

@cli.command(short_help="Install Salmon and related references")
@click.option('--path', '-d', default='./', help='File Path')
@click.option('--species', '-s', default='human', help='Species: human, mouse, zebrafish')
def stats_match_lenght(path, species):
    print('RNA-Seq')
