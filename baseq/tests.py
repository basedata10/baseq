from baseq.cmd import cli
import click

@cli.command(short_help = "Tests for development")
@click.argument("config",nargs=-1)
def dev(config):
    if config[0] == "bed":
        print("[info] BED FILE TESTER...")
        from baseq.bed import BEDFILE
        BEDFILE(config[1])