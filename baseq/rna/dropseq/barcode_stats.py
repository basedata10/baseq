import pandas as pd
from baseq.utils.file_reader import read_file_by_lines

def mutate_single_base(seq):
    mutated = []
    bases = ['A', 'T', 'C', 'G']
    for index in range(0, len(seq)):
        temp = list(seq)
        base_raw = temp[index]
        for base in bases:
            if base != base_raw:
                temp[index] = base
                mutated.append(''.join(temp))
    return mutated

def rev_comp(seq):
    ___tbl = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(___tbl[s] for s in seq[::-1])

class WhiteListCheck:
    def __init__(self, protocol, white_list):
        if protocol == "10X":
            bc_white = [set()]
            with open(white_list, 'r') as infile:
                for line in infile:
                    bc_white[0].add(line.strip())

        if protocol == "indrop":
            bc_white = [set(), set()]
            with open(white_list, 'r') as infile:
                for line in infile:
                    bc_rev = rev_comp(line.strip())
                    bc_white[0].add(bc_rev)

            with open(white_list, 'r') as infile:
                for line in infile:
                    bc_rev = rev_comp(line.strip())
                    bc_white[1].add(bc_rev)

        if protocol == "dropseq":
            bc_white = []
        self.protocol = protocol
        self.bc_white = bc_white

    def check(self, barcode):
        if self.protocol == "10X":
            if barcode in self.bc_white[0]:
                return 1
            else:
                return 0
        if self.protocol == "dropseq":
            if barcode in ['TCAAAAGCAGTG']:
                return 0
            else:
                return 1
        if self.protocol == "indrop":
            if barcode[0:-8] in self.bc_white[0] and barcode[-8:] in self.bc_white[1]:
                return 1
            else:
                return 0

def barcode_aggregate(white_list="10X", protocol="", barcode_count="", max_cell=20000, min_reads=2000, output="./bc_stats.txt"):
    print("[info] Stats the barcodes counts in {}".format(barcode_count))
    df = pd.read_csv(barcode_count).sort_values("counts", ascending=False)
    df["mismatch_reads"] = 0
    df["mismatch_bc"] = ""
    df["mutate_last_base"] = 0
    df = df.reset_index(drop=True)
    data = df.values.tolist() #matrix: barcode/reads/mismatch_reads
    bcs = [x[0] for x in data] #barcode names
    #Aggregate by 1 Mismatch
    bc_counts = len(df.index)
    for id in range(0, bc_counts):
        bc = df.iloc[id, 0]
        count = df.iloc[id, 1]
        print(bc, df.iloc[id, 1])

        if df.iloc[id, 1] == 0 or count <= 0.25 * min_reads:
            continue

        bc_mis = mutate_single_base(bc)

        #index for these mismatches
        index = df['barcode'].isin(bc_mis)
        df_mis = df.loc[index].sort_values("counts", ascending=False)
        barcode_mis = df_mis['barcode'].tolist()

        #determine if mutate in the last base
        if len(barcode_mis)>=3 and sum([1 for x in barcode_mis[0:3] if x[0:-1]==bc[0:-1]])==3:
            df.iloc[id, 4] = 1

        df.iloc[id, 2] = sum(df_mis['counts'])
        df.iloc[id, 3] = "_".join(df_mis['barcode'])
        df.loc[index, "counts"] = 0

    print("[info] Filtering the barcodes exceeds number {}".format(min_reads))
    df = df.loc[df['counts'] >= min_reads]
    df.to_csv(output, index=False)

def getUMI(protocol, seq1, UMI_pos):
    if protocol == "10X":
        UMI = seq1[16:26]
    if protocol == "dropseq":
        if 1:
            UMI = seq1[11:19]
        else:
            UMI = seq1[12:20]
    if protocol == "indrop":
        UMI = seq1[len(barcode) + 22:len(barcode) + 22 + 6]
    return UMI

# from baseq.rna.dropseq.barcode_counting import cut_seq_barcode
# def barcode_split(fq1, fq2, barcode_info, protocol, outdir):
#     read1 = read_file_by_lines(fq1, 10 * 1000 * 1000, 4)
#     read2 = read_file_by_lines(fq2, 10 * 1000 * 1000, 4)
#     bc = cut_seq_barcode(read1[1])