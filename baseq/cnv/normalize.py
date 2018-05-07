from subprocess import call
from baseq.setting import r_script_dir
import os, sys
from baseq.mgt import get_config, run_cmd
import pandas as pd

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

def Normalize_GC_py(genome, bincount, path_out):
    """
    Using a python function ...
    """
    dynamic_bin = get_config("CNV_ref_"+genome, "dynamic_bin")
    df_counts = pd.read_table(bincount)
    df = pd.read_table(dynamic_bin)
    df['counts'] = df_counts['counts']
    #First, Aggregate the Bins into 500kb...
    print("[info] Aggregate the bins into 500kb...")
    df = df.groupby(df.index // 10).agg({"chr":"first", "absstart": "mean", "GC": "mean", "length" : "sum", "counts" : "sum"})
    #Do the normalization on Length
    df['norm_counts'] = df['counts']/df['length']
    print(df)
