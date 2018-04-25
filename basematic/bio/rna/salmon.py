import sys, os, time, re
from subprocess import Popen, PIPE, call
from basematic.mgt import get_config, run_cmd

def run_build_salmon_index():
    pass

def run_salmon(fq1, fq2, genome, path="./"):

    salmon = get_config("RNA", "salmon")
    salmon_ref = get_config("RNA_ref_"+genome, "salmon_index")
    gene_map = ""
    outdir = "./outdir"

    # Run salmon
    if os.path.exists(fq1) and os.path.exists(fq2):
        salmon_cmd = [salmon, 'quant', '-i', salmon_ref, '-l A', '-1', fq1, '-2', fq2, '-p 8', '-g', gene_map,
                      '-o', outdir]
    else:
        salmon_cmd = [salmon, 'quant', '-i', salmon_ref, '-l A', '-U', fq1, '-p 8', '-g', gene_map,
                      '-o', outdir]
    run_cmd("bowtie alignment", " ".join(salmon_cmd))

#Get SALMON STATS...
def get_salmon_stats(outdir):
    pass

def build_tpm_table():
    pass