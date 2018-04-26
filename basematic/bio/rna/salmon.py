import sys, os, time, re
from subprocess import Popen, PIPE, call
from basematic.mgt import get_config, run_cmd

def run_build_salmon_index():
    pass

def run_salmon(fq1, fq2, genome, outdir):
    print(fq1, fq2, genome, outdir)
    salmon = get_config("RNA", "salmon")
    salmon_ref = get_config("RNA_ref_"+genome, "salmon_index")
    gene_map = get_config("RNA_ref_"+genome, "gene_map")

    # Run salmon
    if fq1 and fq2 and os.path.exists(fq1) and os.path.exists(fq2):
        salmon_cmd = [salmon, 'quant', '-i', salmon_ref, '-l A', '-1', fq1, '-2', fq2, '-p 8', '-g', gene_map,
                      '-o', outdir]
    elif fq1 and os.path.exists(fq1):
        salmon_cmd = [salmon, 'quant', '-i', salmon_ref, '-l A', '-r', fq1, '-p 8', '-g', gene_map, '-o', outdir]
    else:
        pass

    run_cmd("Salmon Quantification", " ".join(salmon_cmd))

def run_multiple_salmons(path, genome, outdir):
    from basematic.bio.fastq.sample_file import check_sample_files
    samples = check_sample_files(path)
    print(samples)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        print("[info] Create outdir in: {}".format(outdir))
    script_lists = []
    script_main_path = os.path.join(outdir, "work.sh")
    for sample in samples:
        name = sample[0]
        path = os.path.join(outdir, name)
        script_path = os.path.join(path, "work.sh")
        if not os.path.exists(path):
            os.mkdir(path)
            print("[info] Create outdir for sample {} in: {}".format(name, path))
        else:
            print("[info] Outdir for sample {} exists in: {}".format(name, path))
        if sample[2]:
            script = "basematic-RNA run_salmon -1 {} -2 {} -g {} -d {}".format(sample[1], sample[2], genome, path)
        else:
            script = "basematic-RNA run_salmon -1 {} -g {} -d {}".format(sample[1], genome, path)
        with open(script_path, "w") as file:
            file.writelines(script+"\n")
            print("[info] work script written in {}".format(script_path))
        script_lists.append(script_path)
    #write the main script
    with open(script_main_path, "w") as file:
        file.writelines("\n".join(["bash " + x for x in script_lists]))
    print("[info] Main script written in {}".format(script_main_path))

#Get SALMON STATS...
def get_salmon_stats(outdir):
    pass

def build_tpm_table():
    pass