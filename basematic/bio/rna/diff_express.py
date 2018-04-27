from subprocess import call
from basematic.setting import r_script_dir
import os, sys
from basematic.mgt import get_config, run_cmd

def deseq2(tpmfile, groupfile, outpath):
    """ Run DNACopy.R file ...
    """
    Rscript = get_config("RNA", "deseq")
    script = os.path.join(r_script_dir, "DESeq2.R")
    cmd = "{} {} {} {} {}".format(Rscript, script, tpmfile, groupfile, outpath)
    try:
        run_cmd("Normalize ", cmd)
        print("[info] Segment file write to {}".format(path_out))
    except:
        sys.exit("[error] Failed to run the Normalize Rscript ...")