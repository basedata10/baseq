import os, sys
from baseq.mgt import get_config, run_cmd
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

def plot_genome_py(bincount, cbs_path, path_out):
    df = pd.read_table(bincount)
    #Plot The Dots...
    plt.figure(figsize=(10, 1.4))
    plt.margins(x=0, y=0)
    plt.gcf().subplots_adjust(bottom=0.2)
    plt.gcf().subplots_adjust(left=0.03, right=0.97)
    plt.scatter(df.absstart, df.norm_by_GC_Ploidy, edgecolors='dodgerblue', s=1)

    #Plot Genomes
    df_chr_pos = df.groupby(by=["chr"]).agg({"absstart" : "first"})
    for idx, row in df_chr_pos.iterrows():
        plt.axvline(x=row['absstart'], color="grey")

    #Chr Name
    df_chr_pos = df.groupby(by=["chr"]).agg({"absstart": "mean", "start": "mean", "chr":"first"}).sort_values(by="absstart")
    df_chr_pos['offset'] = df_chr_pos.absstart - df_chr_pos.start
    labels = df_chr_pos.iloc[::2].chr.tolist()
    labels = [x.split("chr")[-1] for x in labels]
    plt.xticks(df_chr_pos.iloc[::2].absstart, labels)

    #Add Segs
    df_cbs = pd.read_table(cbs_path)
    for idx, row in df_cbs.iterrows():
        offset = df_chr_pos.loc[row['chrom'], 'offset']
        xmin = row['loc.start']+offset
        xmax = row['loc.end']+offset
        plt.plot([xmin, xmax], [row['CN'], row['CN']], color="red")

    print("[info] The fig path is {}".format("Genome12.png"))
    plt.savefig("Genome12.png")

def plot_genome_multiple(bincount, cbs_path, path_out):
    df = pd.read_table(bincount)

    #Plot The Dots...
    fig, axes = plt.subplots(10, 1, figsize=(12, 14), sharex='col')
    plt.margins(x=0, y=0)
    plt.gcf().subplots_adjust(bottom=0.2)
    plt.gcf().subplots_adjust(left=0.03, right=0.97)

    #Plot Genomes
    for x in range(10):
        #scatter
        axes[x].scatter(df.absstart, df.norm_by_GC_Ploidy, edgecolors='dodgerblue', s=1)
        #chr line
        df_chr_pos = df.groupby(by=["chr"]).agg({"absstart" : "first"})
        for idx, row in df_chr_pos.iterrows():
            axes[x].axvline(x=row['absstart'], color="grey")
        axes[x].xaxis.set_ticks_position('none')
        #Segs
        df_chr_pos = df.groupby(by=["chr"]).agg({"absstart": "mean", "start": "mean", "chr": "first"}).sort_values(
            by="absstart")
        df_chr_pos['offset'] = df_chr_pos.absstart - df_chr_pos.start
        df_cbs = pd.read_table(cbs_path)
        for idx, row in df_cbs.iterrows():
            offset = df_chr_pos.loc[row['chrom'], 'offset']
            xmin = row['loc.start'] + offset
            xmax = row['loc.end'] + offset
            axes[x].plot([xmin, xmax], [row['CN'], row['CN']], color="red")

    #Chr Name
    df_chr_pos = df.groupby(by=["chr"]).agg({"absstart": "mean", "chr":"first"}).sort_values(by="absstart").iloc[::2]
    labels = df_chr_pos.chr.tolist()
    labels = [x.split("chr")[-1] for x in labels]
    axes[9].set_xticks(df_chr_pos.absstart)
    axes[9].set_xticklabels(labels)

    print("[info] The fig path is {}".format("Genome.png"))
    plt.savefig("Genome_20.png")