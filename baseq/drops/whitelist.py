import os, sys

def rev_comp(seq):
    ___tbl = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(___tbl[s] for s in seq[::-1])

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

from baseq.mgt.config import get_config

def read_whitelist(protocol):
    whilelistdir = get_config("Drops", "whitelistDir")
    if protocol == "10X":
        bc_white = [set()]
        white_path = os.path.join(whilelistdir, "whitelist.10X.txt")
        with open(white_path, 'r') as infile:
            for line in infile:
               bc_white[0].add(line.strip())

    if protocol == "indrop":
        bc_white = [set(), set()]
        white_path = os.path.join(whilelistdir, "whitelist.indrop_1.txt")
        with open(white_path, 'r') as infile:
            for line in infile:
                bc_rev = rev_comp(line.strip())
                bc_white[0].add(bc_rev)
        white_path = os.path.join(whilelistdir, "whitelist.indrop_2.txt")
        with open(white_path, 'r') as infile:
            for line in infile:
                bc_rev = rev_comp(line.strip())
                bc_white[1].add(bc_rev)

    if protocol == "dropseq":
        bc_white = []

    return bc_white

def whitelist_check(bc_white, protocol, barcode):
    if protocol == "10X":
        if barcode in bc_white[0]:
            return 1
        else:
            return 0
    if protocol == "dropseq":
        if barcode in ['TCAAAAGCAGTG']:
            return 0
        else:
            return 1
    if protocol == "indrop":
        if barcode[0:-8] in bc_white[0] and barcode[-8:] in bc_white[1]:
            return 1
        else:
            return 0