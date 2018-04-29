import pandas as pd
from baseq.utils.file_reader import read_file_by_lines

def HammingDistance(seq1, seq2):
    return sum([1 for x in zip(seq1, seq2) if x[0] != x[1]])

def cut_seq_barcode(protocol, seq):
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

def count_barcodes(path, output, protocol, min_reads):
    bc_counts = {}
    index = 0
    print("[info] Process the top 10M reads in {}".format(path))
    print("[info] The Method is : {}".format(protocol))
    print("[info] Barcode with less {} reads is discard".format(min_reads))
    lines = read_file_by_lines(path, 10 * 1000 * 1000, 4)
    for line in lines:
        index += 1
        bc = cut_seq_barcode(protocol, line[1])
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

    print("[info] Write the Barcode depth CSV to {}".format(output))
    df = pd.DataFrame(bc_counts_filter, columns=["barcode", "counts"])
    df.to_csv(output, sep=",", index=False)