import sys, os, time, re, json
from subprocess import Popen, PIPE, call
from baseq.mgt import get_config, run_cmd
from baseq.fastq.sample_file import check_sample_files
from baseq.rna.hisat import run_cufflinks
import pandas as pd


script ="""
cd {} 
{} --runThreadN 20 --genomeDir {} --readFilesCommand zcat --readFilesIn {} {} --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif
{} flagstat Aligned.sortedByCoord.out.bam > flagstat.txt
"""
script1="""
cd {} 
{} --runThreadN 20 --genomeDir {} --readFilesCommand zcat --readFilesIn {} --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif
{} flagstat Aligned.sortedByCoord.out.bam > flagstat.txt
"""


def run_star(fq1, fq2, genome, outdir, run=True):
    star = get_config("RNA", "star")
    star_index = get_config("RNA_ref_" + genome, "star_index")
    samtools = get_config("RNA", "samtools")
    # Run hisat, samtools and cufflinks
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        print("[info] Create outdir in: {}".format(outdir))
    if fq1 and fq2 and os.path.exists(fq1) and os.path.exists(fq2):
        star_cmd = script.format(outdir, star, star_index, fq1, fq2,samtools)
    elif fq1 and os.path.exists(fq1):
        star_cmd = script1.format(outdir, star, star_index, fq1,samtools)
    else:
        pass
    cufflinks_cmd = run_cufflinks(genome, method="star")
    if run:
        run_cmd("star analysis", "".join(star_cmd))
        run_cmd("cufflinks analysis", "".join(cufflinks_cmd))
    return star_cmd + "\n" + cufflinks_cmd

def run_multiple_star(path, genome, outdir):
    samples = check_sample_files(path)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        print("[info] Create outdir in: {}".format(outdir))
    script_lists = []
    script_main_path = os.path.join(outdir, "work.sh")
    qsub_main_path = os.path.join(outdir, "qsub_work.sh")
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
            script_cmd = run_star(sample[1], sample[2], genome, path, False)
        else:
            script_cmd = run_star(sample[1],"",genome,path, False)
        with open(script_path, "w") as file:
            file.writelines("#!/bin/bash"+"\n"+script_cmd+"\n")
            print("[info] work script written in {}".format(script_path))
        script_lists.append(script_path)
     # write the main script
    with open(script_main_path, "w") as file:
        file.writelines("\n".join(["bash " + x for x in script_lists]))
    with open(qsub_main_path, "w") as file:
        file.writelines("\n".join(["qsub -cwd -l vf=8g,p=8 " + x for x in script_lists]))
    print("[info] Main script written in {}".format(script_main_path))