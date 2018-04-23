seqs = {{}}
filepath = {{}}
linecount = {{}}

seqcounts = max(len(seqs), 1)

import subprocess

def read_file():
    if filepath.endswith("gz") or filepath.endswith("gzip") or filepath.endswith("gz2"):
        reading = subprocess.Popen(["gunzip", "-c", filepath], stdout=subprocess.PIPE, bufsize=1000000)
        infile = reading.stdout
        while True:
            data = [infile.readline().decode('utf8') for i in range(linecount)]
            if data[0] == "":
                break
            yield data
    else:
        infile = open(filepath, 'r')
        while True:
            data = [infile.readline() for i in range(linecount)]
            if data[0] == "":
                break
            yield data

def logic():
    seqlen = len(seqs[0])
    content = {}

    for base in ['A', 'T', 'C', 'G']:
        content[base] = [0]*seqlen

    for seq in seqs:
        for idx, base in enumerate(seq):
            if base in ['A', 'T', 'C', 'G'] and idx<seqlen:
                content[base][idx] += 1

    for base in ['A', 'T', 'C', 'G']:
        content[base] = [float(x)/seqcounts for x in content[base]]

    return [content[base] for base in ['A', 'T', 'C', 'G']]