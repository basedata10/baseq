import os, sys
from baseq.bam import cli
from baseq.bam.bamtype import BAMTYPE
import click

@cli.command(short_help="check_bam")
@click.option('--bam', '-i', default='./', help='File Path')
@click.option('--bed', '-b', default='', help='Intervals of Interests, bed file')
@click.option('--regions', '-c', default=100, help='Number or regions for stats...')
def infos(bam, bed, regions):
    if not os.path.exists(bam):
        sys.exit("[info] Bam file not exists: {}".format(bam))
    bam = BAMTYPE(bam, bed)
    #stats bases and reads
    bam.stats_bases()
    bam.stats_duplicates()
    #stats on bedfile
    if bed:
        bam.stats_regions()
        bam.stats_region_coverage(regions)

@cli.command(short_help="split reads not in region...")
@click.option('--path', '-d', default='./', help='File Path')
@click.option('--species', '-s', default='human', help='Species: human, mouse, zebrafish')
def split_reads_not_in_region(path, species):
    print('RNA-Seq')

@cli.command(short_help="Samtools Sort")
@click.option('--thread', '-t', default='8', help='Thread')
@click.argument('in_bam')
@click.argument('out_bam')

def sort_bam(in_bam, thread, out_bam):
    print("[info] The sorted bam will be: {}".format(out_bam))
    cmd = "samtools sort -@ {} {} -o{}; samtools index {}".format(thread, in_bam, out_bam, out_bam)
    from baseq.mgt.command import run_cmd
    run_cmd("SortBam", cmd)