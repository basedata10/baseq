import click
from baseq.snv import cli
import pandas as pd

@cli.command(short_help="Check the enrichment quality, input: bam, interval and filename")
@click.option('--bamlist', '-l', default = '', help = 'List of Bam file')
@click.option('--bampath', '-b', default = '', help = 'Bampath')
@click.option('--intervals', '-r', default = '', help = '')
@click.option('--name', '-n', default = 'Sample', help = 'Outfile (Excel/Table)')
def qc_enrich(bamlist, bampath, intervals, name):
    from baseq.snv.qc.enrich import quality_of_enrich_sample
    if bamlist:
        with open(bamlist, 'r') as file:
            lines = file.readlines()
        bams = [line.strip().split() for line in lines]
    else:
        bams = [[name, bampath]]

    print("[info] {} Bam files".format(len(bams)))
    results = []
    import multiprocessing as mp
    pool = mp.Pool(processes = 5)
    for bam in bams:
        results.append(pool.apply_async(quality_of_enrich_sample, (bam[0], bam[1], intervals)))
    pool.close()
    pool.join()

    results = [x.get() for x in results]
    df = pd.DataFrame(results, columns=["Sample", "Total", "Mapped", "Map_Ratio", "Dup_ratio", "Mean_Depth", "PCT_10X", "PCT_30X", "PCT_50X", "PCT_100X"])
    print("[info] Write QC Infos to {}".format("QC.xls"))
    df.to_excel("QC.xls")

