from baseq.rna import cli
import click, os, sys

@cli.command(short_help="Run DESeq2")
@click.option('--groupfile', '-g', default='', help="Tab seprated file: samplename, groups")
@click.option('--comparefile', '-p', default='', help="Tab seprated file: group1, group2")
@click.option('--tpmfile', '-t', default='./', help="TPM filepath")
@click.option('--countfile', '-c', default='./', help="Read Count filepath")
@click.option('--outpath', '-o', default='./diff_exp', help="Folder for exportation...")
def diff_deseq2(tpmfile, countfile, groupfile, comparefile, outpath):
    from baseq.rna.deseq.deseq2 import deseq2
    print("[info] DESeq2 {}".format(outpath))
    deseq2(tpmfile, countfile, groupfile, comparefile, outpath)

@cli.command(short_help="Run DESeq2")
@click.option('--groupfile', '-g', default='', help="Tab seprated file: samplename, groups")
@click.option('--comparefile', '-p', default='', help="Tab seprated file: group1, group2")
@click.option('--outdir', '-o', default='./diff_exp', help="Folder for exportation...")
def diff_deseq2_result(groupfile, comparefile, outdir):
    from baseq.rna.deseq.deseq2 import pack_DeSeq2_result
    print("[info] DESeq2 {}".format(outdir))
    pack_DeSeq2_result(groupfile, comparefile, outdir)