import click, re, os
from baseq.snv import cli
import pandas as pd
import numpy as np

@cli.command(short_help="Stats The VCF File")
@click.option('--vcf', '-i', default = '', help = 'VCF path')
@click.option('--vcf_lists', '-l', default = '', help = 'VCF path List, sample name and path')
@click.option('--depth', '-d', default = 100, help = 'Filter Depth')
@click.option('--name', '-n', default = 'sample', help = 'Name of process')

def vcf_stats(vcf, vcf_lists, depth, name):
    from baseq.snv.vcf.GATK import vcf_stats

    #build VCF lists
    if vcf and os.path.exists(vcf):
        vcfs = [name, vcf]
    elif vcf_lists and os.path.exists(vcf_lists):
        with open(vcf_lists, 'r') as infile:
            lines = infile.readlines()
        infos = [re.split("\s+", x) for x in lines]
        vcfs = [x for x in infos if x[0] and os.path.exists(x[1])]
    else:
        pass

    results = []
    import multiprocessing as mp
    pool = mp.Pool(processes=20)

    for vcf in vcfs:
        results.append(pool.apply_async(vcf_stats, (vcf[0], vcf[1], int(depth),)))
    pool.close()
    pool.join()
    results = [x.get() for x in results]
    MAF = [[x['sample']] + x["MAF"].tolist() for x in results]

    writer = pd.ExcelWriter('VCF_stats.xlsx', engine='xlsxwriter', options={'font_name':'arial'})

    # pd.DataFrame(results, columns=["sample", "counts", "mean_depth", "GT_01", "GT_02"]).to_excel("VCF.xls")
    # pd.DataFrame(MAF, columns=["sample"]+[str(round(x/50, 2)) for x in range(50)]).to_excel("MAF.xls")

    pd.DataFrame(results, columns=["sample", "counts", "mean_depth", "GT_01", "GT_11"]).to_excel(writer, sheet_name='VCF')
    pd.DataFrame(MAF, columns=["sample"] + [str(round(x/100, 2)) for x in range(50)]).to_excel(writer, sheet_name='MAF')
    writer.book.formats[0].set_font_name('arial')
    writer.save()