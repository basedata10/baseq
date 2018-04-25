from subprocess import call
from basematic.setting import r_script_dir
import os, sys
from basematic.mgt import get_config, run_cmd

def bin_segmentation(infile, path_out):
    """ Run DNACopy.R file ...
    """
    script = os.path.join(r_script_dir, "DNACopy.R")
    cmd = "Rscript {} {} {}".format(script, infile, path_out)
    try:
        run_cmd("Normalize ", cmd)
        print("[info] Segment file write to {}".format(path_out))
    except:
        sys.exit("[error] Failed to run the Normalize Rscript ...")

