import json
import pandas as pd
from basematic.file.Read import Path_Iterator

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
    print("[error] Starts stating the barcodes counts from {}".format(barcode_count))
    df = pd.read_csv(barcode_count).sort_values("counts", ascending=False)
    df["mismatch_reads"] = 0
    df.set_index('barcode')

    # Aggregate by 1 Mismatch
    bc_counts = len(df.index)
    for id in range(0, bc_counts):
        bc = df.loc[id, 'barcode']
        bc_mis = mutate_single_base(bc)
        df_mis = df.loc[df['barcode'].isin(bc_mis)]
        df.loc[id, 'mismatch_reads'] = sum(df_mis['counts'])
        df.loc[df['barcode'].isin(bc_mis), 'counts'] = 0
    print("[info] Filtering the barcodes exceeds number {}".format(min_reads))
    df = df.loc[df['counts']>=min_reads]

    # Detemine UMI Start Pos

    df["UMI_pos"] = 0
    df.to_csv(output)

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

from basematic.bio.rnaseq.barcode import get_barcode
def barcode_split(fq1, fq2, barcode_info, protocol, outdir):
    read1 = Path_Iterator(fq1, 10 * 1000 * 1000, 4)
    read2 = Path_Iterator(fq2, 10 * 1000 * 1000, 4)
    bc = get_barcode(read1[1])