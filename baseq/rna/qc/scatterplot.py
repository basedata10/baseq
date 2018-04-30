import os, sys

def plot_corelation_fig(name1, name2, exp_table, figname):
    """
    plot a figure with name:
    """
    import pandas as pd
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt

    tpm_file = pd.read_table(exp_table)
    sample_names = list(tpm_file.columns.values)
    path = os.path.join(figname + '.png')

    print("[info] Figure axis log scale, min 1, max 5000;")
    print("[info] Figure save to {}".format(path))

    if name1 and name2 in sample_names:
        plot_data = tpm_file[(tpm_file[name1] > 0) | (tpm_file[name2] > 0)]
        plt.figure(figsize=(5, 5))
        plt.scatter(plot_data[name1], plot_data[name2], s=3)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlim(1, 5000)
        plt.ylim(1, 5000)
        plt.savefig(path)

    elif not name1 in sample_names:
        print("[error] undefined sample name {}".format(name1))

    elif not name2 in sample_names:
        print("[error] undefined sample name {}".format(name2))