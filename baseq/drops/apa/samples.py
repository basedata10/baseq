import os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import matplotlib.style as style

def APA_usage(bamfile, APA_sitefile, celltype, gene):
    """ Get the abundance for each cell barcode for each APA in the gene.

    :param bamfile: method for the new :class:`Request` object.
    :param APA_sitefile: URL for the new :class:`Request` object.
    :param celltype: (optional) The celltype file genreated from cellranger
    :type gene: string

    Usage:
      >>> import requests
      >>> req = requests.request('GET', 'http://httpbin.org/get')
      <Response [200]>

    Returnsass:
        Generate a heatmap;
        Print the Read count;
    """
    from baseq.bam import BAMTYPE
    bam = BAMTYPE(bamfile)

    #Read CellType Table...
    if celltype:
        df_type = pd.read_csv(celltype)
        df_type["cell"] = [x.split("-")[0] for x in df_type.Barcode.tolist()]
        df_type = df_type.drop("Barcode", axis=1)
        df_type = df_type.set_index('cell')

    #Read APA Site Table...
    df_apa = pd.read_table(APA_sitefile)
    df_gene = df_apa[df_apa.gene == gene]
    sample_usage = []

    #Get The Mapped Read Infos For Each Peak...
    for idx, row in df_gene.iterrows():
        chr = row['chr']
        start = row['pos']-100
        end = row['pos']+100
        reads = bam.get_reads(chr, start, end)
        for read in reads:
            read_header = read[0].split("_")
            sample_usage.append([read_header[1], read_header[2], str(idx)])

    #Build a Table
    df_counts = pd.DataFrame(sample_usage, columns=["sample", "UMI", "APA"])
    df_counts['reads'] = 1
    df_counts = df_counts.groupby(by=["sample", "UMI", "APA"]).sum().reset_index()

    df_counts = df_counts.drop(["UMI"], axis=1)
    df_counts = df_counts.groupby(by=["sample", "APA"]).count().reset_index()
    df_counts = df_counts.pivot(index='sample', columns='APA', values='reads').fillna(0)
    df_counts["total"] = df_counts.sum(axis=1)
    df_counts = df_counts[df_counts.total>=1]
    df_counts = df_counts.sort_values("total", ascending=False)

    #Aggregate By Cell Type...
    if celltype:
        df = df_counts.join(df_type)
        df = df.groupby("Cluster").sum()
        print(df)
        df = df.div(df.total/100, axis=0)
        print(df)

    #plot heatmap....
    style.use('seaborn-poster')
    plt.figure()
    df_counts = df_counts.drop(["total"], axis=1)
    sns.heatmap(df_counts.iloc[1:40, :], cmap="YlGnBu_r")
    plt.savefig("hehe.png")
    print("[info] Figure Export To {}".format("hehe.png"))