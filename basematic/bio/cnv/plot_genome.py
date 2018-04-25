from subprocess import call
from basematic.setting import r_script_dir
import os, sys
from basematic.mgt import get_config, run_cmd

def PlotGenomes(bincount, cbs_path, path_out):
    """Plot the genome
    """
    script = os.path.join(r_script_dir, "Genome_Plot.R")
    #dynamic_bin = get_config("CNV_ref_"+genome, "dynamic_bin")
    cmd = "Rscript {} {} {} {}".format(script, bincount, cbs_path, path_out)
    try:
        run_cmd("Plot Genome", cmd)
        #build the json...
    except:
        sys.exit("[error] Failed to run the Normalize Rscript ...")