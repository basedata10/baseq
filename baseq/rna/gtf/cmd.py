from baseq.rna import cli
import click, os, sys


docs = """
[get UTR]
baseq-RNA gtf_read -g test -t UTR

[Scan Genome...]

"""

@cli.command(short_help="Run DESeq2")
@click.option('--genome', '-g', default='hg38', help="Genome")
@click.option('--type', '-t', default='gene', help="Type of Structure")
def gtf_read(genome, type):
    from baseq.rna.gtf.gencode import read_gencode
    df = read_gencode(genome, type)
    #df = df.loc[df['']]
    print(df)

@cli.command(short_help="Run DESeq2")
@click.option('--genome', '-g', default='hg38', help="Genome")
def gtf_UTR3(genome):
    from baseq.rna.gtf.gencode import read_gencode
    df = read_gencode(genome, "UTR")
    df = df.groupby("transc").last()
    df['length'] = df.end-df.start
    df = df.loc[df.length>1000]
    print(df)
    return df