# -*- coding: utf-8 -*-
import itertools, json, re, sys, os
import random

def HammingDistance(seq1, seq2):
    return sum([1 for x in zip(seq1, seq2) if x[0] != x[1]])

def sampling_UMIs(UMIs, ratio):
    UMIs_sample = {}
    for UMI in UMIs:
        UMIs_sample[UMI] = sum([1 for x in range(UMIs[UMI]) if random.random()<=ratio])
    return {x:UMIs_sample[x] for x in UMIs_sample if UMIs_sample[x]>0}

def calculate_UMI_with_mismatch(UMIs):
    """
    Corrected the mismatches in UMIs
    input: UMI sequences and their counts;
    return: Corrected unique UMI sequences
    """
    if len(UMIs.keys()) == 1:
        return [x for x in UMIs if UMIs[x]>0]

    UMIs = sorted(UMIs.items(), key=lambda k: k[1], reverse=True)
    UMI_info = {x[0]:x[1] for x in UMIs}
    umi_num = len(UMIs)
    if umi_num <= 10:
        for idx1 in range(0, umi_num-1):
            for idx2 in range(idx1+1, umi_num):
                umi_1 = UMIs[idx1][0]
                umi_2 = UMIs[idx2][0]
                if HammingDistance(umi_1, umi_2) <= 1:
                    UMI_info[umi_1] += UMI_info[umi_2]
                    UMI_info[umi_2] = 0
    return [x for x in UMI_info if UMI_info[x]>0]

def read_barcode_gene_file(filepath):
    """
    Return the UMIs and Counts for each aggragate file...
    [barcode in the file, Reads in the file and UMIs in the file]
    """
    sampling_factor = 1
    UMIs = {}
    Reads = {}
    barcodes = []

    with open(filepath, 'r') as file:
        for line in file:
            infos = re.split('\t', line)
            barcode = infos[0]
            barcodes.append(barcode)
            umi_counts = json.loads(infos[1])
            for gene in umi_counts:
                gene_UMIs = sampling_UMIs(umi_counts[gene], sampling_factor)
                bcs = calculate_UMI_with_mismatch(gene_UMIs)
                reads = sum(gene_UMIs.values())
                if gene in UMIs:
                    UMIs[gene][barcode] = len(bcs)
                    Reads[gene][barcode] = Reads
                else:
                    UMIs[gene] = {}
                    UMIs[gene][barcode] = len(bcs)
                    Reads[gene] = {}
                    Reads[gene][barcode] = reads
    return [barcodes, UMIs, Reads]

#Export Genes Which express in at least 10 samples
def write_to_table(barcodes, UMIs, outfile):
    header = "\t".join(["genes"]+barcodes)
    print(len(barcodes), barcodes)
    with open(outfile, 'w') as file:
        file.write(header+"\n")
        for gene in UMIs:
            umis = [gene]
            if len(UMIs[gene].keys())<= 0.02*len(barcodes) or sum(UMIs[gene].values())<=0.1*len(barcodes):
                continue
            for bc in barcodes:
                if bc in UMIs[gene]:
                    umis.append(UMIs[gene][bc])
                else:
                    umis.append(0)
            umis_str = "\t".join([str(x) for x in umis])
            file.write(umis_str+"\n")
    print("[info] Success Write to file {}".format(outfile))