import sys, os, time, re
from subprocess import Popen, PIPE, call
from basematic.mgt import get_config, run_cmd

def bin_counting(genome, bamfile, out):

    import pandas as pd
    import bisect

    dynamic_bin = get_config("CNV_ref_"+genome, "dynamic_bin")
    df_raw = pd.read_table(dynamic_bin, names=["chr", "start", 'absstart', 'end', 'range', 'length', 'GC'],sep=" ")

    chrs = df_raw["chr"].tolist()
    start = df_raw["start"].tolist()
    end = df_raw["end"].tolist()

    #Build region position for each chromosome...
    chr_bin_starts = {}
    chr_bin_ends = {}
    chr_bin_idx = {}
    for idx, chr in enumerate(chrs):
        if chr not in chr_bin_starts:
            chr_bin_starts[chr] = [start[idx]]
            chr_bin_ends[chr] = [end[idx]]
            chr_bin_idx[chr] = [idx]
        else:
            chr_bin_starts[chr].append(start[idx])
            chr_bin_ends[chr].append(end[idx])
            chr_bin_idx[chr].append(idx)

    lines = len(chrs)
    counts = [0]*lines

    #Read bam file
    reading = Popen(["samtools", "view", bamfile], stdout=PIPE, bufsize=1000000)
    infile = reading.stdout

    while True:
        data = infile.readline().decode('utf8')
        if data == "":
            break
        infos = data.split()

        #Filter on quality
        quality = int(infos[4])
        if quality<40:
            continue
        # if re.search("XS:i", data):
        #     continue

        #Counting reads to bins
        chr = infos[2]
        pos = int(infos[3])
        if not chr in chr_bin_starts:
            continue

        idx_chr = bisect.bisect_right(chr_bin_starts[chr], pos)

        if idx_chr<len(chr_bin_ends[chr]) and chr_bin_ends[chr][idx_chr]>=pos:
            idx = chr_bin_idx[chr][idx_chr]
            counts[idx] += 1

    #print
    df_raw["counts"] = counts
    print("[info] Save the bin count file to {}".format(out))
    df_raw['counts'].to_csv(out, header=True, sep="\t")