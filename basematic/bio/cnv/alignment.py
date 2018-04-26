import sys, os, time, re
from subprocess import Popen, PIPE, call
from basematic.mgt import get_config, run_cmd

def script(samples):
    from basematic.bio.fastq import sample_file as FQfiles
    samples = FQfiles.check_sample_files(samples)
    if len(samples) == 0:
        sys.exit("[error] No valid samples")
    else:
        print("[info] {} samples detected".format(len(samples)))

def bowtie2_sort_alignment(fastq, genome, path="./"):
    bowtie2 = get_config("CNV", "bowtie2")
    bowtie2_ref = get_config("CNV_ref_"+genome, "bowtie2_index")
    samtools = get_config("CNV", "samtools")

    samfile = os.path.join(path, "sample.sam")
    bamfile = os.path.join(path, "sample.bam")
    statsfile = os.path.join(path, "sample.stats.txt")

    print("[info] The bamfile will be write to: {}".format(bamfile))

    #Run Bowtie
    bowtie_cmd = [bowtie2, '-p', '10', '-x', bowtie2_ref, '-u', str(5*1000*1000), '-U', fastq, '>', samfile]
    run_cmd("bowtie alignment", " ".join(bowtie_cmd))

    #run Samtools
    samtools_sort_cmd = [samtools, 'sort -@ 8 -o', bamfile, samfile, ";", samtools, "index", bamfile, "; rm", samfile]
    run_cmd("samtools sort", " ".join(samtools_sort_cmd))

    #run flagstats
    cmd_stats = [samtools, "flagstat", bamfile, ">", statsfile]
    run_cmd("samtools stats", " ".join(cmd_stats))
    return bamfile

