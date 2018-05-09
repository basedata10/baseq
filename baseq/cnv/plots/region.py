from baseq.setting import r_script_dir
import os, sys
from baseq.mgt import get_config, run_cmd

def plot_genome(bincount, cbs_path, path_out):
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