import os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import peakutils
import pandas as pd

from baseq.drops.apa.scaner import scan

def scan_utr(genome, bam, name):
    """
    For a genome, read the gencode annotationm get the logest UTR for each gene (>=1000bp)\
    Apply the 'scan' function for each UTR (default 20 threads...)\
    Call the peaks for each UTR.\
    Build and Write the APA Peaks for all the genes.


    """
    from baseq.bam import BAMTYPE
    from baseq.rna.gtf.gencode import read_gencode

    # Read Gencode And Get UTR Region (>1000bp)
    # Sort by UTR length and get the logest by gene...
    df = read_gencode(genome, "UTR")
    df['length'] = df.end - df.start
    df = df.sort_values(by=['length'])
    df = df.groupby("gene").last()
    df = df.loc[df.length>1000]
    bam = BAMTYPE(bam)

    import multiprocessing as mp
    pool = mp.Pool(processes=20)
    results = []
    for index, row in df.iterrows():
        results.append(pool.apply_async(scan, (bam, index, row['chr'], row['start'], row['end'], 50)))
    pool.close()
    pool.join()

    results = [x.get() for x in results]
    results = [y for x in results for y in x]
    results = [df.loc[x[0],:].tolist() + x  for x in results]
    df_peaks = pd.DataFrame(results, columns=["chr", "start", 'end', 'strand', 'transc',
                                              'exon', 'length', 'gene', 'pos', 'mean_depth',
                                              "mid", "left", "right", "counts"])
    df_peaks = df_peaks.drop(columns=["transc", 'exon'])
    file_tsv = "peaks.{}.txt".format(name)
    file_xls = "peaks.{}.xls".format(name)
    df_peaks.to_csv("peaks.{}.txt".format(name), sep="\t")
    df_peaks.to_excel("peaks.{}.xls".format(name))
    print("[info] The peaks files are write to: {}/{}".format(file_tsv, file_xls))