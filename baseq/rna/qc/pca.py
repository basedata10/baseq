import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import seaborn as sns
from baseq.rna.filetype import *
from .exp_table import log2_transform

def pca_analysis(expTable, groupfile):
    """
    PCA for QC
    """
    #read groupfile
    groups = pd.read_table(groupfile, sep=" ")
    print("[info] Groups lists ...")
    print(groups)
    samples = groups["sample"]

    df = pd.read_table(expTable, index_col=0)
    df = df[samples]
    df = df.loc[df.mean(axis=1)>=10, :]
    df = log2_transform(df)
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(df.transpose())
    pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'])

    pca_df['group'] = groups['groups']
    sns_plot = sns.pairplot(data=pca_df, x_vars=["PC1"], y_vars=["PC2"],  hue="group", size=5)
    sns_plot.savefig("PCA.png")

def pca_score_power(exp_table, samplegroup, group_pairs, qcfile):
    """
    For each group pair, do the PCA analysis between them...
    """
    #prepare Expression files ...
    df_exp = pd.read_table(exp_table, index_col=0)
    df_sample = pd.read_table(samplegroup, sep='\s+')
    df_compare = pd.read_table(group_pairs, sep='\s+')
    df_qc = pd.read_table(qcfile, sep="\s+", index_col="sample")

    for idx, row in df_compare.iterrows():
        samples1 = df_sample.loc[df_sample["groups"]== row.tolist()[0], "sample"].tolist()
        samples2 = df_sample.loc[df_sample["groups"] == row.tolist()[1], "sample"].tolist()
        df = df_exp[samples1+samples2]
        df = df.loc[df.mean(axis=1) >= 10, :]
        df = log2_transform(df)
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(df.transpose())
        pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'], index=samples1+samples2)

        #QC
        qcs = df_qc.loc[samples1+samples2]
        corr = pca_df['PC1'].corr(qcs['ratio'])
        corr2 = pca_df['PC2'].corr(qcs['ratio'])
        print("PC1, PC2 Vs Ratio Correlation", corr, corr2)

def pca_score_exploration():
    pass