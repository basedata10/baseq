import time, subprocess
import pandas as pd
try:
    import cPickle as pickle
except:
    import _pickle as pickle

def HammingDistance(seq1, seq2):
    return sum([1 for x in zip(seq1, seq2) if x[0] != x[1]])

def get_barcode(protocol, seq):
    if protocol == "10X":
        return seq[0:16]
    if protocol == "dropseq":
        return seq[0:12]
    if protocol == "indrop":
        w1 = "GAGTGATTGCTTGTGACGCCTT"
        if w1 in seq:
            w1_pos = seq.find(w1)
            if 7 < w1_pos < 12:
                return seq[0:w1_pos] + seq[w1_pos + 22:w1_pos + 22 + 8]
        else:
            for i in range(8, 12):
                w1_mutate = seq[i:i + 22]
                if HammingDistance(w1_mutate, w1) < 2:
                    return seq[0:i] + seq[i + 22 : i + 22 + 8]
                    break
        return ""

from basematic.file.Read import Path_Iterator

def getBarcode(path, output, protocol, min_reads):
    bc_counts = {}
    index = 0
    print("[info] Process the top 10M reads in {}".format(path))
    print("[info] The protocol is {}".format(protocol))
    lines = Path_Iterator(path, 10*1000*1000, 4)
    for line in lines:
        index += 1
        bc = get_barcode(protocol, line[1])
        if bc == "":
            continue
        if bc in bc_counts:
            bc_counts[bc] += 1
        else:
            bc_counts[bc] = 1
        if index % 1000000 == 1:
            print("[info] Processed {} lines".format(index-1))

    bc_counts_filter = []
    for k, v in bc_counts.items():
        if v >= min_reads:
            bc_counts_filter.append([k, v])

    print("[info] Write the Barcode depth file to {}".format(output))
    df = pd.DataFrame(bc_counts_filter, columns=["barcode", "counts"])
    df.to_csv(output, sep=",", index=False)

def splitBarcode(outfile, prefix=2):
    pass