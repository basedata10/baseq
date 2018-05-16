from baseq.drops.cmd import cli
from baseq.mgt.config import get_config
import click, os

@cli.command(short_help="Genome position to gene name")
@click.option('--dir', '-d', default='', help="Genome:hg38...")
@click.option('--genome', '-g', default='hg38', help="Genome:hg38...")
def apa_pipeline(dir):
    pass

@cli.command(short_help="Genome position to gene name")
@click.option('--bam', '-b', default="aligned.bam", help='Path to the bam file')
@click.option('--chr', '-c', default="chr", help='Chrome')
@click.option('--start', '-s', help='Start')
@click.option('--end', '-e', help='End')
def apa_scan_region(bam, chr, start, end):
    """
    Get the peaks for chrN:start-end and print the peaks
    """
    from baseq.bam import BAMTYPE
    from baseq.drops.apa.scaner import scan
    import pandas as pd
    bam = BAMTYPE(bam)
    infos = scan(bam, "sample", chr, int(start), int(end))
    print(pd.DataFrame(infos))

@cli.command(short_help="Genome position to gene name")
@click.option('--genome', '-g', default='hg38', help="Genome:hg38...")
@click.option('--bam', '-b', default="aligned.bam", help='Path to the bam file')
@click.option('--name', '-n', default="", help='Name of analysis (ScanAPA)')
def apa_scan_genome(genome, bam, name):
    """

    """
    from baseq.drops.apa.UTR import scan_utr
    scan_utr(genome, bam, name)

@cli.command(short_help="APA Usage For Barcode or CellType")
@click.option('--bamfile', '-b', default='', help="Bam File")
@click.option('--apalist', default='', help="APA sites recovered")
@click.option('--celltype', '-t', default='', help="Cellbarcode and the its type")
@click.option('--gene', '-g', default='', help="APA sites")
def apa_sample_usage(bamfile, apalist, celltype, gene):
    from baseq.drops.apa.samples import APA_usage
    APA_usage(bamfile, apalist, celltype, gene)

@cli.command(short_help="Genome position to gene name")
@click.option('--bamfile', '-b', default='', help="Bam File")
@click.option('--apalist', default='', help="APA sites recovered")
@click.option('--gene', '-g', default='', help="APA sites")
def apa_reverse_strand_events(bamfile, apalist, gene):
    from baseq.drops.apa.samples import APA_usage
    APA_usage(bamfile, apalist, gene)