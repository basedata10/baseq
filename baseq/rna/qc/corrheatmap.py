import os, sys

def get_house_keeping_genes_expression():
    pass

def correlation_heatmap(table, name):
    """Spearman Correlation Heatmap Plot"""
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