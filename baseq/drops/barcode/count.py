import pandas as pd
from time import time
from baseq.utils.file_reader import read_file_by_lines

def HammingDistance(seq1, seq2):
    return sum([1 for x in zip(seq1, seq2) if x[0] != x[1]])

def extract_barcode(protocol, seq):
    """Extract cell barcode from reads

    - 10X: seq[0:16]
    - indrop: seq[0:i] + seq[i + 22 : i + 22 + 8] (i is length of barcode 1)
    - dropseq: seq[0:12]

    :param protocol: 10X/indrop/drop-seq.
    :param seq: The sequence containing cellbarcode.

    Return:
        barcode: barcode, if no valid barcode, return ""
    """

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

def count_barcodes(path, output, protocol, min_reads, topreads=100):
    """Count thre number of Each barcode

    :param path: fastq file.
    :param output: The stats will write to ...
    :param protocol: Protocol
    :param min_reads: minimum reads
    :param topreads: process max N million reads

    Return:
        A barcode_count file will be generated.
        cellbarcode/counts
    """

    bc_counts = {}
    index = 0
    start = time()
    print("[info] Process the top {}M reads in {}".format(topreads, path))
    print("[info] Barcode with less than {} reads is discard".format(min_reads))
    lines = read_file_by_lines(path, topreads * 1000 * 1000, 4)
    for line in lines:
        index += 1
        bc = extract_barcode(protocol, line[1])
        if index % 1000000 == 0:
            print("[info] Processed {}M lines in {}s".format(index/1000000, round(time()-start, 2)))
            start = time()
        if bc == "":
            continue
        if bc in bc_counts:
            bc_counts[bc] += 1
        else:
            bc_counts[bc] = 1

    bc_counts_filter = []
    for k, v in bc_counts.items():
        if v >= min_reads:
            bc_counts_filter.append([k, v])

    print("[info] Barcode depth file: {}".format(output))
    df = pd.DataFrame(bc_counts_filter, columns=["barcode", "counts"])
    df.to_csv(output, sep=",", index=False)