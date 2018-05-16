import os, sys
import pandas as pd
import numpy as np
from time import time
from baseq.utils.file_reader import read_filelines
from baseq.drops.barcode.count import extract_barcode

from itertools import product
barcode_prefix = [x[0]+x[1] for x in list(product('ATCG', repeat=2))]

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
        files[barcode] = open(path, 'w')
        buffers[barcode] = []
    return (files, buffers)

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

def barcode_split(name, protocol, bcstats, fq1, fq2, dir, size, skip):
    print("[info] SPLIT!!! {}".format(name))
    #barcode infos
    barcode_corrected = {}
    barcode_mutate_last = []
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

    #make buffers and file...
    files, buffers = open_splited_files(dir, name)
    fq1 = read_filelines(fq1, size, 4, skip)
    fq2 = read_filelines(fq2, size, 4, skip)

    counter = 0
    start = time()
    for read1 in fq1:
        read2 = next(fq2)
        seq1 = read1[1]
        counter += 1
        if counter % 1000000 == 0:
            print("[info] Processed {}M lines in {}s".format(counter/1000000, round(time()-start, 2)))
            start = time()

        if counter % 100000 == 0:
            for key in buffers:
                files[key].writelines("\n".join(buffers[key])+"\n")
                buffers[key] = []

        bc = extract_barcode(protocol, read1[1])
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
        header = "_".join(['@', bc_corrected, UMI, str(counter+skip)])
        seq = read2[1].strip()
        quality = read2[3].strip()
        buffers[bc_prefix].append("\n".join([header, seq, "+", quality]))

    for key in buffers:
        files[key].writelines("\n".join(buffers[key])+"\n")
        files[key].close()
    return 1

def barcode_splits_all(name, protocol, bcstats, fq1, fq2, minreads, maxcell, dir):
    import multiprocessing as mp
    pool = mp.Pool(processes=10)
    results = []
    size = 10 * 1000 * 1000
    names = []
    for x in range(18):
        skip = x * size
        name2 = "{}_{}".format(name, x)
        names.append(name2)
        results.append(pool.apply_async(barcode_split, [name2, protocol, bcstats, fq1, fq2, dir, size, skip]))
    pool.close()
    pool.join()
    [x.get() for x in results]
    #combine file...
    bases = [x+y for x in "ATCG" for y in "ATCG"]
    for base in bases:
        files = [os.path.join(dir, "split.{}.{}.fq".format(name, base)) for name in names]
        files = " ".join(files)
        merged_file = os.path.join(dir, "split.{}.{}.fq".format(name, base))
        cmd = "cat {} > {}".format(files, merged_file)
        print(cmd)