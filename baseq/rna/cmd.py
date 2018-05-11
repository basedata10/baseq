import click
from baseq.rna import cli

@cli.command(short_help="Show the docs and examples ...")
def doc():
    from .docs import print_doc
    print_doc()

@cli.command(short_help="salmon quatification")
@click.option('--sample_file', '-m', default='', help="Tab seprated file: name, fq1, fq2")
@click.option('--genome', '-g', default='hg38', help="Species hg19 or mm10/mm38")
@click.option('--name', '-n', default='RNA', help="Name of The Process (RNA)")
@click.option('--parallel', '-p', default='6', help="How many runs at the same time...")
@click.option('--fq1', '-1', default='', help="Fastq 1")
@click.option('--fq2', '-2', default='', help="Fastq 2")
def run_salmon(sample_file, fq1, fq2, genome, name, parallel):
    from .salmon import run_salmon, run_multiple_salmons
    if sample_file:
        run_multiple_salmons(sample_file, genome, name, parallel)
    else:
        run_salmon(fq1, fq2, genome, name)

#run hisat, samtools and cufflinks
@cli.command(short_help="hisat alignment and cufflinks quatification")
@click.option('--sample_file', '-m', default='', help="Tab seprated file: name, fq1, fq2")
@click.option('--genome', '-g', default='hg38', help="Species hg19 or mm10/mm38")
@click.option('--fq1', '-1', default='', help="Fastq 1")
@click.option('--fq2', '-2', default='', help="Fastq 2")
@click.option('--outdir', '-d', default='', help="The scripts and output path")
def run_hisat(sample_file, fq1, fq2, genome, outdir):
    from baseq.rna.hisat import run_hisat, run_multiple_hisat
    if sample_file:
        run_multiple_hisat(sample_file, genome, outdir)
    else:
        run_hisat(fq1, fq2, genome, outdir)

@cli.command(short_help="STAR alignment and cufflinks quatification")
@click.option('--sample_file', '-m', default='', help="Tab seprated file: name, fq1, fq2")
@click.option('--genome', '-g', default='hg38', help="Species hg19 or mm10/mm38")
@click.option('--fq1', '-1', default='', help="Fastq 1")
@click.option('--fq2', '-2', default='', help="Fastq 2")
@click.option('--outdir', '-d', default='', help="The scripts and output path")
def run_star(sample_file, fq1, fq2, genome, outdir):
    from baseq.rna.star import run_star, run_multiple_star
    if sample_file:
        run_multiple_star(sample_file, genome, outdir)
    else:
        run_star(fq1, fq2, genome, outdir)

@cli.command(short_help="Aggregate TPM and Counts and QC Tables for multiple samples")
@click.option('--sample_file', '-m', default='', help="Tab seprated file: name, fq1, fq2")
@click.option('--processdir', '-d', default='./', help="Combine all the TPMs under the folder")
@click.option('--name', '-n', default='sample', help="Name of Analysis")
def aggr_TPM_QC(processdir, sample_file, name):
    from .salmon import build_tpm_table
    build_tpm_table(processdir, sample_file, name)

@cli.command(short_help="Aggregate TPM and Counts and QC Tables for multiple samples")
@click.option('--sample_file', '-m', default='', help="Tab seprated file: name, fq1, fq2")
@click.option('--processdir', '-d', default='./', help="Combine all the TPMs under the folder")
@click.option('--outpath', '-o', default='./', help="Prefix of the TPM and Count file")
def aggr_FPKM_QC(processdir, sample_file, outpath):
    from .hisat import build_FPKM_table
    print("[info] Aggregate FPKM into {}".format(outpath))
    build_FPKM_table(processdir, sample_file, outpath)

@cli.command(short_help="plot tpm correlation figure between samples")
@click.option('--name1', '-1', default='', help="sample name")
@click.option('--name2', '-2', default='', help="sample name")
@click.option('--table', '-t', default='', help="Table path")
@click.option('--figname', '-n', default='scatter', help="The scripts and output path")
def plot_corelation_scatter(name1, name2, table, figname):
    from baseq.rna.qc.scatterplot import plot_corelation_fig
    plot_corelation_fig(name1, name2, table, figname)

@cli.command(short_help="Correlation Heatmap")
@click.option('--table', '-t', default='', help="Table path")
@click.option('--name', '-n', default='heatmap', help="The scripts and output path")
def corr_heatmap(table, name):
    from baseq.rna.qc.corrheatmap import correlation_heatmap
    correlation_heatmap(table, name)

@cli.command(short_help="Correlation Heatmap")
@click.option('--table', '-t', default='', help="Table path")
@click.option('--group', '-g', default='', help="Group file")
@click.option('--name', '-n', default='pca_analysis', help="The scripts and output path")
def pca_analysis(table, group, name):
    from baseq.rna.qc.pca import pca_analysis
    pca_analysis(table, group)

@cli.command(short_help="Correlation Heatmap")
@click.option('--table', '-t', default='', help="Table path")
@click.option('--group', '-g', default='', help="Group file")
@click.option('--comparefile', '-p', default='', help="Tab seprated file: group1, group2")
@click.option('--qcfile', '-q', default='', help="QC files")
@click.option('--name', '-n', default='pca_analysis', help="The scripts and output path")
def diff_power_analysis(table, group, comparefile, qcfile, name):
    from baseq.rna.qc.pca import pca_score_power
    pca_score_power(table, group, comparefile, qcfile)

from .table.cmd import *
from .deseq.cmd import *
from .gtf.cmd import *