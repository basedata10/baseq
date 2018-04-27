import click, sys, os
from basematic.bio.fastq.sample_file import check_sample_files

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    click.echo("Welcome to Basematic-RNA")

@cli.command(short_help="Show the docs and examples ...")
def doc():
    from .docs import print_doc
    print_doc()

#Run the CNV pipeline for Lists of Samples
@cli.command(short_help="RNA-Seq")
@click.option('--sample_file', '-m', default='', help="Tab seprated file: name, fq1, fq2")
@click.option('--genome', '-g', default='hg19', help="Species hg19 or mm10/mm38")
@click.option('--fq1', '-1', default='', help="Fastq 1")
@click.option('--fq2', '-2', default='', help="Fastq 2")
@click.option('--outdir', '-d', default='./process_salmon', help="The scripts and output path")
def run_salmon(sample_file, fq1, fq2, genome, outdir):
    from .salmon import run_salmon, run_multiple_salmons
    if sample_file:
        run_multiple_salmons(sample_file, genome, outdir)
    else:
        run_salmon(fq1, fq2, genome, outdir)

@cli.command(short_help="Aggregate TPM and Counts and QC Tables for multiple samples")
@click.option('--sample_file', '-m', default='', help="Tab seprated file: name, fq1, fq2")
@click.option('--processdir', '-d', default='./', help="Combine all the TPMs under the folder")
@click.option('--outpath', '-o', default='./tpm.txt', help="Prefix of the TPM and Count file")
def aggr_TPM_QC(sample_file, processdir, outpath):
    from .salmon import build_tpm_table
    print("[info] Aggregate TPM into {}".format(outpath))
    build_tpm_table(processdir, sample_file, outpath)


@cli.command(short_help="Run DESeq2")
@click.option('--groupfile', '-g', default='', help="Tab seprated file: sample, group")
@click.option('--tpmfile', '-t', default='./', help="TPM filepath")
@click.option('--outpath', '-o', default='./tpm.txt', help="Prefix of the TPM and Count file")
def deseq2(tpmfile, groupfile, outpath):
    from .diff_express import deseq2
    print("[info] DESeq2 {}".format(outpath))
    deseq2(tpmfile, groupfile, outpath)

@cli.command(short_help="inDrop/Drop-Seq/10X")
@click.option('--genome', help="human/mouse/mixed")
@click.option('--name', '-n', default='sample', help="sample name")
@click.option('--fq1', '-1', default='', help='fastq1 path')
@click.option('--fq2', '-2', default='', help='fastq2 path (optional)')
@click.option('--dir', '-d', default='./', help='Folder of output (./)')
def run_drops(name, fq1, fq2, genome, dir):
    print('Start Processing inDrop Results')
    samples = check_sample_files("", name, fq1, fq2)
    if samples == []:
        sys.exit("[error] No valid sample, Exit.")
    from basematic.bio.rna.barcode import getBarcode
    from basematic.bio.rna.barcode_stats import barcode_aggregate
    getBarcode(fq1, "./barcode_count.csv", "10X", 20)
    barcode_aggregate(barcode_count="./barcode_count.csv")