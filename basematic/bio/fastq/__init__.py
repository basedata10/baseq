import click, os, sys

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(context_settings=CONTEXT_SETTINGS)

def cli():
    click.echo("Basematic Start...")

@cli.command(short_help="Quality Control")
@click.option('--fq1', default='', help='fastq')
@click.option('--fq2', default='', help='fastq2')
@click.option('--outtype', '-t', default='table', help='(table)/csv/excel')
@click.option('--outfile', '-o', default='./qc.txt', help='file path')
def Quality(fq1, fq2, outtype, outfile):
    import pandas as pd
    print("Qualtity control of the Fastq Files ...")
    print("fq1 is ", fq1)
    print("fq2 is ", fq2)
    from basematic.bio.fastq.Stats import Stats_Fastq
    res = []
    if fq1 and os.path.exists(fq1):
        res.append(Stats_Fastq(fq1))
    else:
        sys.exit("Exit, fq1 is not Valid")
    if fq2:
        if os.path.exists(fq2):
            res.append(Stats_Fastq(fq2))
        else:
            sys.exit("Exit, fq2 is not Valid")

    print("Qualtity control results are write to :", outfile)
    pd.concat(res).to_csv(outfile, sep="\t")

@cli.command(short_help="Quality Control")
@click.option('--path', '-d', default='./', help='File Path')
@click.option('--ratio', '-d', default='1', help='ratio of reads')
@click.option('--reads', '-d', default='1', help='read counts')
@click.option('--tofasta', '-d', default='1', help='read counts')
@click.option('--path', '-d', default='./', help='File Path')
def Sampling(path):
    from basematic.bio.fastq.Stats import Stats_Fastq
    click.echo('Start Using CNV Analysis')
    print(Stats_Fastq(path))

@cli.command(short_help="List the Files")
@click.option('--path', '-d', default='./', help='File Path')
def Lists(path):
    print("List The Files... {}".format(path))