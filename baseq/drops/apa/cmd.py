#besides: barcode-->gene-->counts
#we should generate: gene-->counts [read positions]

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

"""
baseq-Drop apa_scan_region -b star.CC_Aligned.sort.bam -c chr1 -s 12141952 -e 12145429
baseq-Drop apa_scan_region -b star.CC_Aligned.sort.bam -c chr1 -s 12204917 -e 12211873
baseq-Drop apa_scan_region -b Merged.sort.A.bam -c chr8 -s 102648432 -e 102650458
#UBR5
baseq-Drop apa_scan_region -b Merged.sort.A.bam -c chr8 -s 102251429 -e 102255483
"""
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
    index = get_config("RNA_ref_"+genome, "gencode_bowtie2")
    depth = BAMTYPE(bam).region_depth(chr, start, end, all=True)
    scan(depth, start, end)


@cli.command(short_help="Genome position to gene name")
@click.option('--genome', '-g', default='hg38', help="Genome:hg38...")
@click.option('--bam', '-b', default="aligned.bam", help='Path to the bam file')
def apa_scan_genome(genome, bam):
    """
    Report ALL THE APA CONDICATE...
    GENE/START/END/PA_START/PA_END/STRAND/MEANDEPTH/
    """
    from baseq.bam import BAMTYPE
    from baseq.drops.apa.scaner import scan
    from baseq.rna.gtf.gencode import read_gencode
    df = read_gencode(genome, "UTR")
    df['length'] = df.end-df.start
    df = df.groupby("transc").last()
    df = df.sort_values(by=['length'])
    df = df.groupby("gene").last()
    df = df.loc[df.length>1000]
    bam = BAMTYPE(bam)

    for index, row in df.iterrows():
        depth = bam.region_depth(row['chr'], row['start'], row['end'], all=True)
        if depth == []:
            continue
        if max(depth)<100:
            continue
        else:
            peaks = scan(depth, row['start'], row['end'], mindepth=100)
            if len(peaks)>=3:
                print(index, sum(depth)/len(depth), len(depth), peaks)
                print("baseq-Drop apa_scan_region -b Merged.sort.A.bam -c {} -s {} -e {}".format(row['chr'], row['start'], row['end']))


