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
@click.option('--outpath', '-o', default='./', help="Prefix of the TPM and Count file")
def aggr_TPM_QC(processdir, sample_file, outpath):
    from .salmon import build_tpm_table
    print("[info] Aggregate TPM into {}".format(outpath))
    build_tpm_table(processdir, sample_file, outpath)


@cli.command(short_help="Run DESeq2")
@click.option('--groupfile', '-g', default='', help="Tab seprated file: samplename, groups")
@click.option('--comparefile', '-p', default='', help="Tab seprated file: group1, group2")
@click.option('--tpmfile', '-t', default='./', help="TPM filepath")
@click.option('--countfile', '-c', default='./', help="Read Count filepath")
@click.option('--outpath', '-o', default='./diff_exp', help="Folder for exportation...")
def deseq2(tpmfile, countfile, groupfile, comparefile, outpath):
    from .diff_express import deseq2
    print("[info] DESeq2 {}".format(outpath))
    deseq2(tpmfile, countfile, groupfile, comparefile, outpath)

@cli.command(short_help="Run DESeq2")
@click.option('--groupfile', '-g', default='', help="Tab seprated file: samplename, groups")
@click.option('--comparefile', '-p', default='', help="Tab seprated file: group1, group2")
@click.option('--outdir', '-o', default='./diff_exp', help="Folder for exportation...")
def deseq2_result(groupfile, comparefile, outdir):
    from .diff_express import pack_DeSeq2_result
    print("[info] DESeq2 {}".format(outdir))
    pack_DeSeq2_result(groupfile, comparefile, outdir)


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

@cli.command(short_help="plot tpm correlation figure between samples")
@click.option('--name1', '-1', default='', help="sample name")
@click.option('--name2', '-2', default='', help="sample name")
@click.option('--table', '-t', default='', help="Table path")
@click.option('--outdir', '-d', default='./process_salmon', help="The scripts and output path")
def Plot_corelation_fig(name1, name2, table, outdir):
    from .salmon import Plot_corelation_fig
    Plot_corelation_fig(name1,name2, table)

@cli.command(short_help="Correlation Heatmap")
@click.option('--table', '-t', default='', help="Table path")
@click.option('--name', '-n', default='heatmap', help="The scripts and output path")
def corr_heatmap(table, name):
    from .quality_control import correlation_heatmap
    import pandas as pd
    import numpy as np
    df = pd.read_table(table, index_col=0)
    list(df)
    correlations = df.corr()
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cax = ax.matshow(correlations, vmin=0, vmax=1)
    fig.colorbar(cax)
    ticks = np.arange(0, 48, 4)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels([list(df)[x] for x in ticks])
    ax.set_yticklabels([list(df)[x] for x in ticks])
    plt.xticks(rotation=90)
    plt.savefig("./"+name+".png")

# @cli.command(short_help="slice|")
# def tools(name1, name2, table, outdir):
#     from .salmon import Plot_corelation_fig
#     Plot_corelation_fig(name1,name2, table)

@cli.command(short_help="Genrate Excel....")
def Excel():
    import xlsxwriter

    # Create an new Excel file and add a worksheet.
    workbook = xlsxwriter.Workbook('demo.xlsx')

    #Set default...
    workbook.formats[0].set_font_size(12)
    workbook.formats[0].set_font_name('arial')

    worksheet = workbook.add_worksheet("Sheet1")

    format_main = workbook.add_format({'bold': False, 'font_size':12, 'font_name':'arial'})
    format_header = workbook.add_format({'bold': True, 'font_size':15, 'font_name':'arial'})

    # Widen the first column to make the text clearer.
    worksheet.set_column('A:A', 20)

    # Write some simple text.
    worksheet.write('A1', 'Hello')

    # Text with formatting.
    worksheet.write('A2', 'World', format_header)

    # Write some numbers, with row/column notation.
    worksheet.write(2, 0, "XXXX", format_main)
    worksheet.write(3, 0, 123.456, format_main)
    worksheet.write(3, 0, 12341234, format_main)

    # Sheet 2...
    worksheet2 = workbook.add_worksheet("Sheet2")
    worksheet2.write('A2', 'World', format_header)
    worksheet2.write(10, 10, 'World 10 10')

    # Insert an image.
    worksheet.insert_image('B5', 'logo.png')

    workbook.close()