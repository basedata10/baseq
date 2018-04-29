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

def plot_corelation_fig(name1, name2, counttable, figname):
    import pandas as pd
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt

    tpm_file = pd.read_table(counttable)
    sample_names = list(tpm_file.columns.values)

    print("[info] Figure axis log scale, min 1, max 5000;")

    if name1 and name2 in sample_names:
        plot_data = tpm_file[(tpm_file[name1] > 0) | (tpm_file[name2] > 0)]
        plt.figure(figsize=(5, 5))
        plt.scatter(plot_data[name1], plot_data[name2], s=3)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlim(1, 5000)
        plt.ylim(1, 5000)
        plt.savefig(os.path.join(figname + '.png'))

    elif not name1 in sample_names:
        print("[error] undefined sample name {}".format(name1))

    elif not name2 in sample_names:
        print("[error] undefined sample name {}".format(name2))