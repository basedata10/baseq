from subprocess import call
from basematic.setting import r_script_dir
import os, sys

def CBS_GC():
    """Lowess call a Rscript named as 'Lowess.R'.
    Input:
        Raw Read Counts
        Genome GC content file
    Process:
        Read Counts are first normalized by mean;
        Then normalized by GC using Lowess Regression;
    Output:

    """
    script = os.path.join(r_script_dir, "DNACopy.R")

