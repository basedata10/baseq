from baseq.utils.file_reader import read_file_by_lines
import pandas as pd
pd.set_option('precision', 3)

def fastq_basecontent_quality(sample, fastq_path, maxLines = 10000):
    """
    Generate the basic quality stats of the fastq file
    return:
        dataframe: A/T/C/G/quality;
        base content figure in base64;
        base quality figure in base64;
    """
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt

    inlines = read_file_by_lines(fastq_path, maxLines, 4)
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
                if not (ord(base) - 33):
                    print(ord(base) - 33)
                quality[idx] += ord(base) - 33
    #high_bias_pos: position where one base exceed 50% of all coverage;
    high_bias_pos = 0
    for base in ['A', 'T', 'C', 'G']:
        content[base] = [float(x) / len(seqs) for x in content[base]]
        high_bias_pos += sum([1 for x in content[base] if x>=0.5])

    #quality and mean quality of all bases...
    content['quality'] = [q / len(seqs) for q in quality]
    mean_quality= round(sum(content['quality'])/len(content['quality']), 2)

    #plot basecontent...
    plt.figure(figsize=(4, 2))
    plt.plot(range(1, seqlen+1), content['A'])
    plt.ylim((0, 1))
    plt.savefig("./{}_basecontent.png".format(sample))

    #plot quality
    plt.figure(figsize=(4, 2))
    plt.plot(range(1, seqlen+1), content['quality'])
    plt.savefig("./{}_basequality.png".format(sample))

    pd.DataFrame(content, columns=['A', 'T', 'C', 'G', 'quality'])
    return ("./{}_basecontent.png".format(sample), "./{}_basequality.png".format(sample), mean_quality, high_bias_pos)