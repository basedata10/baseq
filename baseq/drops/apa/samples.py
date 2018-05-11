import os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import matplotlib.style as style

def assign_reads_to_APA(file):
    pass

def APA_usage(bamfile, APA_sitefile, gene):
    from baseq.bam import BAMTYPE
    bam = BAMTYPE(bamfile)
    df_apa = pd.read_table(APA_sitefile)
    df_gene = df_apa[df_apa.gene == gene]
    sample_usage = []
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
    print(df_counts)
    df_counts = df_counts.groupby(by=["sample", "APA"]).count().reset_index()
    print(df_counts)
    df_counts = df_counts.pivot(index='sample', columns='APA', values='reads').fillna(0)
    df_counts["total"] = df_counts.sum(axis=1)
    df_counts = df_counts[df_counts.total>=3]
    df_counts = df_counts.sort_values("total", ascending=False)
    print(df_counts)

    #plot heatmap....
    style.use('seaborn-poster')
    plt.figure()
    df_counts = df_counts.drop(["total"], axis=1)
    sns.heatmap(df_counts.iloc[1:40,:], cmap="YlGnBu_r")
    plt.savefig("hehe.png")