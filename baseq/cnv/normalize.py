from subprocess import call
from baseq.setting import r_script_dir
import os, sys
from baseq.mgt import get_config, run_cmd

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