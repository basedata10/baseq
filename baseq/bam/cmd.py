import os, sys, re, click
from baseq.bam import cli
from baseq.bam.bamtype import BAMTYPE
from baseq.mgt.command import run_cmd

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
def unmapped(path, species):
    print('RNA-Seq')

@cli.command(short_help="Get reads in regions")
@click.option('--bamfile', '-b', default='', help='Bam File Path')
@click.option('--pos', '-p', default='', help='Position (chrN:XX-XX)')
@click.option('--chr', '-c', default='', help='Chrome')
@click.option('--start', '-s', default='', help='Start')
@click.option('--end', '-e', default='', help='End')
def region_reads(bamfile, pos, chr, start, end):
    if pos:
        pos = re.sub(",", "", pos)
        pos = re.split(":|-", pos)
        chr = pos[0]
        start = pos[1]
        end = pos[2]
    else:
        start = int(start)
        end = int(end)
    datas = BAMTYPE(bamfile).get_reads(chr, start, end)
    print(datas)

@cli.command(short_help="Samtools Sort")
@click.option('--thread', '-t', default='8', help='Thread')
@click.argument('in_bam')
@click.argument('out_bam')

def sort(in_bam, thread, out_bam):
    if out_bam[-4:] == ".bam":
        out_bam = out_bam[0:-4]
    print("[info] The sorted bam will be: {}".format(out_bam+".bam"))
    cmd = "samtools sort -@ {} {} {}; samtools index {}".format(thread, in_bam, out_bam, out_bam+".bam")
    run_cmd("SortBam", cmd)

@cli.command(short_help="Samtools Merge")
@click.option('--thread', '-t', default='8', help='Thread')
@click.argument('in_bam')
@click.argument('out_bam')

def merge(in_bam, thread, out_bam):
    if out_bam[-4:] == ".bam":
        out_bam = out_bam[0:-4]
    print("[info] The sorted bam will be: {}".format(out_bam+".bam"))
    cmd = "samtools sort -@ {} {} {}; samtools index {}".format(thread, in_bam, out_bam, out_bam+".bam")
    run_cmd("SortBam", cmd)