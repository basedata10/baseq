import pandas as pd
from baseq.utils.file_reader import read_file_by_lines

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

from baseq.rna.dropseq.barcode_counting import cut_seq_barcode
def barcode_split(fq1, fq2, bcstats, protocol, outdir):
    #read

    fq1 = read_file_by_lines(fq1, 10 * 1000 * 1000, 4)
    fq2 = read_file_by_lines(fq2, 10 * 1000 * 1000, 4)

    for read1 in fq1:
        read2 = next(fq2)
        bc = cut_seq_barcode(read1[1])

