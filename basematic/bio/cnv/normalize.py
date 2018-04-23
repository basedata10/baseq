from subprocess import call
from basematic.setting import r_script_dir
import os, sys
from basematic.mgt import get_config, run_cmd

def Normalize_GC(genome, bincount):
    """Lowess call a Rscript named as 'Lowess.R'.
    Input:
        Raw Read Counts
        Genome GC content file
    Process:
        Read Counts are first normalized by mean;
        Then normalized by GC using Lowess Regression;
    Output:

    """
    script = os.path.join(r_script_dir, "Lowess.R")
    dynamic_bin = get_config("CNV_ref_"+genome, "dynamic_bin")
    cmd = "Rscript {} {} {}".format(script, bincount, dynamic_bin)
    run_cmd("Normalize ", cmd)