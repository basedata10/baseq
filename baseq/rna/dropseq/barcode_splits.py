import os, sys
import pandas as pd
import numpy as np
from baseq.utils.file_reader import read_file_by_lines
from baseq.rna.dropseq.barcode_counting import cut_seq_barcode

from itertools import product
barcode_prefix = [x[0]+x[1] for x in list(product('ATCG', repeat=2))]

#build 16 files
def open_splited_files(dir, name):
    files = {}
    buffers = {}
    if not os.path.exists(dir):
        try:
            os.mkdir(dir)
        except:
            sys.exit("[error] Failed to make output dir...")

    for barcode in barcode_prefix:
        path = os.path.join(dir, "split.{}.{}.fq".format(name, barcode))
        files[barcode] = open(path, 'wt')
        buffers[barcode] = []
    return (files, buffers)

def write_buffers():
    pass

def getUMI(protocol, barcode, seq1, mutate_last_base):
    if protocol == "10X":
        UMI = seq1[16:26]
    if protocol == "dropseq":
        if mutate_last_base:
            UMI = seq1[11:19]
        else:
            UMI = seq1[12:20]
    if protocol == "indrop":
        UMI = seq1[len(barcode) + 22:len(barcode) + 22 + 6]
    return UMI

def barcode_split(name, protocol, bcstats, fq1, fq2, minreads, maxcell, dir):
    #barcode infos
    barcode_corrected = {}
    barcode_mutate_last = []
    #read barcode stats...
    bc_stats = pd.read_csv(bcstats)
    bc_stats = bc_stats.replace(np.nan, '', regex=True)
    for index, row in bc_stats.iterrows():
        barcode = row['barcode']
        barcode_corrected[barcode] = barcode
        mutate_last = row['mutate_last_base']
        if mutate_last:
            barcode_mutate_last.append(barcode)
        if row['mismatch_bc']:
            barcode_mis = row['mismatch_bc'].split("_")
            for bc in barcode_mis:
                barcode_corrected[bc] = barcode
    print(barcode_corrected)
    files, buffers = open_splited_files(dir, name)

    fq1 = read_file_by_lines(fq1, 1 * 1000 * 1000, 4)
    fq2 = read_file_by_lines(fq2, 1 * 1000 * 1000, 4)
    counter = 0
    for read1 in fq1:
        read2 = next(fq2)
        seq1 = read1[1]

        counter += 1
        bc = cut_seq_barcode(protocol, read1[1])
        if bc in barcode_corrected:
            bc_corrected = barcode_corrected[bc]
        else:
            continue

        bc_prefix = bc_corrected[0:2]
        if bc_prefix not in barcode_prefix:
            continue

        #get UMI
        if barcode in barcode_mutate_last:
            UMI = getUMI(protocol, barcode, seq1, 1)
        else:
            UMI = getUMI(protocol, barcode, seq1, 0)

        #build new header
        header = "_".join(['@', bc_corrected, UMI, str(counter)])
        seq = read2[1].strip()
        quality = read2[3].strip()

        buffers[bc_prefix].append("\n".join([header, seq, "+", quality]))

        if counter % 100000 ==1:
            print("[info] Processed 10K reads")
            for key in buffers:
                files[key].writelines("\n".join(buffers[key]))
                buffers[key] = []

    for key in buffers:
        files[key].writelines("\n".join(buffers[key]))
        files[key].close()