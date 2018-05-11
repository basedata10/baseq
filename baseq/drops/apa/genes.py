import os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

def genes_with_APA(file):
    df = pd.read_table(file)
    df["left_ratio"] = df.left/df.mid
    df["right_ratio"] = df.right/df.mid
    print(df.iloc[1:10,])
    #Filter Exon Events
    df_raw = df.loc[(df.left_ratio>=0.5) & (df.left_ratio>=0.5), :]

    #Aggregate By Gene
    df = df_raw.groupby(["gene"]).agg({"pos":[max, min, len], "mid":["mean"], "strand":["first"]})
    df.add_suffix('_apa').reset_index()
    df.columns = ["max_apa", "min_apa", "counts", "mean_depth", "strand"]
    df['distance'] = df.max_apa - df.min_apa
    df = df.loc[(df.distance>=1000) & (df.mean_depth>=100)]

    genes_pos = df.loc[df.strand=="+"].index.tolist()
    genes_neg = df.loc[df.strand=="-"].index.tolist()
    print(genes_pos, genes_neg)

    #Get_F_Strand
    df_gene = df_raw[df_raw.gene.isin(genes_pos)]
    g1_l = df_gene.sort_values("pos").groupby("gene").first()
    g1_r = df_gene.sort_values("pos").groupby("gene").last()
    df_gene = df_raw[df_raw.gene.isin(genes_neg)]
    g2_l = df_gene.sort_values("pos").groupby("gene").last()
    g2_r = df_gene.sort_values("pos").groupby("gene").first()

    df_1 = pd.DataFrame([g1_l.pos, g1_r.pos, g1_l.mid, g1_r.mid], index=["left", "right", "depth_l", "depth_r"]).transpose()
    df_2 = pd.DataFrame([g2_l.pos, g2_r.pos, g2_l.mid, g2_r.mid], index=["left", "right", "depth_l", "depth_r"]).transpose()
    df_2 = df_2.append(df_1)

    plt.figure(figsize=(10, 10))
    plt.subplot(2, 2, 1)
    plt.scatter(df_2.depth_l, df_2.depth_r)
    plt.loglog()

    plt.subplot(2, 2, 2)
    sns.distplot(df_2.depth_l, hist=False, rug=True, hist_kws={'log':True}).set(xlim=(0, 1000))
    sns.distplot(df_2.depth_r, hist=False, rug=True, hist_kws={'log':True}).set(xlim=(0, 1000))

    plt.subplot(2, 2, 3)
    plt.hist(df_2.depth_l, density=True, color="r", bins=50)
    plt.hist(df_2.depth_r, density=True, color="b", bins=50)
    plt.loglog()

    plt.subplot(2, 2, 4)
    sns.distplot(df_2.depth_r/df_2.depth_l, hist=False, rug=True).set(xlim=(0, 10))

    print(df_2.depth_l.median(), df_2.depth_r.median())
    df_2.to_csv("Biomodule.genes.txt", sep="\t")
    plt.savefig("Plots.png")