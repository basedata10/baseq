from basematic.file.Read import Path_Iterator
import pandas as pd
pd.set_option('precision', 3)

def fastq_basecontent_quality(fastq_path, maxLines = 100000):
    """
    Generate the basic quality stats of the fastq file
    return:
        dataframe: A/T/C/G/quality;
        base content figure in base64;
        base quality figure in base64;
    """
    inlines = Path_Iterator(fastq_path, maxLines, 4)
    seqs = []
    quals = []
    content = {}
    for line in inlines:
        seqs.append(line[1].strip())
        quals.append(line[3].strip())

    seqlen = len(seqs[0])
    quality = [0] * seqlen
    bases = ['A', 'T', 'C', 'G']

    for base in bases:
        content[base] = [0] * seqlen

    #count bases
    for seq in seqs:
        for idx, base in enumerate(seq):
            if base in ['A', 'T', 'C', 'G'] and idx < seqlen:
                content[base][idx] += 1
    #count quality
    for qual in quals:
        for idx, base in enumerate(qual):
            if idx < seqlen:
                quality[idx] += ord(base) - 33

    content['qual'] = [int(q / len(seqs)) for q in quality]
    for base in ['A', 'T', 'C', 'G']:
        content[base] = [float(x) / len(seqs) for x in content[base]]

    #plot figures

    return pd.DataFrame(content, columns=['A', 'T', 'C', 'G', 'quality'])