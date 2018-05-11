from baseq.drops.cmd import cli
from baseq.mgt.config import get_config
import click, os

@cli.command(short_help="Genome position to gene name")
@click.option('--dir', '-d', default='', help="Genome:hg38...")
@click.option('--genome', '-g', default='hg38', help="Genome:hg38...")
def apa_pipeline(dir):
    pass

@cli.command(short_help="Genome position to gene name")
@click.option('--genome', '-g', default='hg38', help="Genome:hg38...")
@click.option('--fastq', '-f', help='Path to the Fastq File')
@click.option('--bam', '-b', default="aligned.bam", help='Path to the bam file')
def apa_align(fastq, genome, bam):
    from baseq.align.bowtie2 import bowtie2_sort
    print('[info] Aligning the reads ...')
    index = get_config("RNA_ref_"+genome, "gencode_bowtie2")
    bowtie2_sort(fastq, "", bam, index, 1*1000*1000)

@cli.command(short_help="Genome position to gene name")
@click.option('--genome', '-g', default='hg38', help="Genome:hg38...")
@click.option('--bam', '-b', default="aligned.bam", help='Path to the bam file')
@click.option('--chr', '-c', default="chr", help='Chrome')
@click.option('--start', '-s', help='Start')
@click.option('--end', '-e', help='End')
def apa_scan_region(genome, bam, chr, start, end):
    start = int(start)
    end = int(end)
    from baseq.bam import BAMTYPE
    from baseq.drops.apa.scaner import scan
    print('[info] Aligning the reads ...')
    depth = BAMTYPE(bam).region_depth(chr, start, end, all=True)
    scan(depth, start, end)

@cli.command(short_help="Genome position to gene name")
@click.option('--genome', '-g', default='hg38', help="Genome:hg38...")
@click.option('--bam', '-b', default="aligned.bam", help='Path to the bam file')
def apa_scan_genome(genome, bam):
    from baseq.drops.apa.scaner import scan_genome
    scan_genome(genome, bam)

@cli.command(short_help="Genome position to gene name")
@click.option('--file', '-f', default='', help="APA sites")
def apa_filter_genes_with_apa(file):
    from baseq.drops.apa.genes import genes_with_APA
    genes_with_APA(file)

@cli.command(short_help="Genome position to gene name")
@click.option('--bamfile', '-b', default='', help="Bam File")
@click.option('--apalist', default='', help="APA sites recovered")
@click.option('--gene', '-g', default='', help="APA sites")
def apa_sample_usage(bamfile, apalist, gene):
    from baseq.drops.apa.samples import APA_usage
    APA_usage(bamfile, apalist, gene)

@cli.command(short_help="Genome position to gene name")
@click.option('--bamfile', '-b', default='', help="Bam File")
@click.option('--apalist', default='', help="APA sites recovered")
@click.option('--gene', '-g', default='', help="APA sites")
def apa_reverse_strand_events(bamfile, apalist, gene):
    from baseq.drops.apa.samples import APA_usage
    APA_usage(bamfile, apalist, gene)