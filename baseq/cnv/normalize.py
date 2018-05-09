from subprocess import call
from baseq.setting import r_script_dir
import os, sys
from baseq.mgt import get_config, run_cmd
import pandas as pd
import numpy as np
import statsmodels.api as sm
from scipy.interpolate import interp1d
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

def Normalize_GC(genome, bincount, path_out):
    """Normalize the Raw bin counts with GC contents...
    Run a Rscript named as 'Lowess.R'.
    Args:
        genome : Raw Read Counts
        bincount : Genome GC content file
    Process:
        Read Counts are first normalized by mean;
        Then normalized by GC using Lowess Regression;
    Output:
        GC_content_image: images
        Normalized bin counts (1M)
    """
    script = os.path.join(r_script_dir, "Lowess.R")
    dynamic_bin = get_config("CNV_ref_"+genome, "dynamic_bin")
    cmd = "Rscript {} {} {} {}".format(script, bincount, dynamic_bin, path_out)
    try:
        run_cmd("Normalize ", cmd)
    except:
        sys.exit("[error] Failed to run the Normalize Rscript ...")

def Normalize_GC_py(genome, bincount, outpath):
    """
    Using a python function ...
    """
    dynamic_bin = get_config("CNV_ref_"+genome, "dynamic_bin")
    df_counts = pd.read_table(bincount)
    df = pd.read_table(dynamic_bin)
    df['counts'] = df_counts['counts']
    #First, Aggregate the Bins into 500kb...
    print("[info] Aggregate the bins into 500kb...")
    df = df.groupby(df.index // 10).agg({"chr":"first", "start": "mean", "absstart": "mean", "GC": "mean", "length" : "sum", "counts" : "sum"})

    #Do the normalization on Length
    df['norm_counts'] = df['counts']/df['length']
    df['norm_counts'] = df['norm_counts']/df['norm_counts'].mean()
    df['norm_counts_log'] = np.log(df.norm_counts + 0.01)

    #Do Norm on GC...
    lowess = sm.nonparametric.lowess
    z = lowess(df.norm_counts_log, df.GC)
    f = interp1d(list(zip(*z))[0], list(zip(*z))[1], bounds_error = False)
    df['norm_by_GC'] = np.exp(df.norm_counts_log - f(df['GC']))-0.01
    df['norm_by_GC'] = df['norm_by_GC']/df['norm_by_GC'].mean()

    #plot GC correction...
    plt.figure(figsize=(10, 10))
    plt.subplot(2, 2, 1)
    plt.title('Raw Normalized Reads (500Kb)')
    plt.scatter(df.GC, df.norm_counts, facecolors='none', edgecolors='r')

    plt.subplot(2, 2, 2)
    plt.title('GC corrected (500Kb)')
    plt.scatter(df.GC, df.norm_by_GC, facecolors='none', edgecolors='b')

    #Peaks Detection
    plt.subplot(2, 2, 3)
    plt.title('Peaks')
    plt.hist(df.norm_by_GC, bins=300, density=True)

    #Peaks Detection..
    Ploidy_Lists = [x/40 for x in range(60, 240, 1)]
    SoS = []
    for ploidy in Ploidy_Lists:
        errors = round(df['norm_by_GC']*ploidy)-df['norm_by_GC']*ploidy
        SoS.append(sum([x*x for x in errors]))
    estimate_ploidy = Ploidy_Lists[SoS.index(min(SoS))]

    print("[info] The estimated plodity is {}".format(estimate_ploidy))
    plt.subplot(2, 2, 4)
    plt.title('Plodity Estimate: {}'.format(estimate_ploidy))
    plt.ylabel("Errors (Should be Minimized)")
    plt.xlabel("Ploidy")
    plt.plot(Ploidy_Lists, SoS)
    df['norm_by_GC_Ploidy'] =  df['norm_by_GC']*estimate_ploidy

    df_export = df[['chr', 'start', 'absstart', 'norm_by_GC', 'norm_by_GC_Ploidy']]
    print("[info] Write to {}".format(outpath))
    df_export.to_csv(outpath, sep="\t", float_format='%.2f')

    plt.savefig("Normalize.png")
    print(df)