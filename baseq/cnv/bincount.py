from subprocess import Popen, PIPE, call
from baseq.mgt import get_config, run_cmd
import pandas as pd

def bin_counting(genome, bamfile, out):
    import bisect

    dynamic_bin = get_config("CNV_ref_"+genome, "dynamic_bin")
    df_bin = pd.read_table(dynamic_bin, sep="\t")

    chrs = df_bin["chr"].tolist()
    start = df_bin["start"].tolist()
    end = df_bin["end"].tolist()

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

        #Counting reads to bins
        chr = infos[2]
        pos = int(infos[3])
        if not chr in chr_bin_starts:
            continue

        idx_chr = bisect.bisect_right(chr_bin_starts[chr], pos)

        if idx_chr<len(chr_bin_ends[chr]) and chr_bin_ends[chr][idx_chr]>=pos:
            idx = chr_bin_idx[chr][idx_chr]
            counts[idx] += 1

    df_bin["counts"] = counts
    print("[info] Bin Count File: {}".format(out))
    df_bin['counts'].to_csv(out, header=True, sep="\t")