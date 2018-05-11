import os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import peakutils
import pandas as pd

def scan(bam, gene, chr, start, end, min_depth=10):
    """
    The returned result:

    """
    print("[info] {}".format(gene))
    start = start-100
    end = end + 100
    #Get Depth
    depth = bam.region_depth(chr, start, end, all=True)

    if depth == []:
        return []

    if max(depth) < 100:
        return []

    #Average Every 20 bases
    kernel_0 = [1] * 20
    depth_window_slice = np.convolve(depth, kernel_0,  'same')/20

    #Left Right Difference
    kernel_1 = [-1]*25 + [0] + [1] * 25
    depth_diffs = np.convolve(depth_window_slice, kernel_1,  'same')/50

    #Flanking Difference Patterns
    kernel_2 = [1]*150 + [0] + [-1] * 150
    res2 = np.convolve(depth_diffs, kernel_2,  'same')/300

    indexes = peakutils.indexes(res2, thres=0.02 / max(res2), min_dist=100)

    x = range(start, end+1)

    plt.figure(figsize=(10, 7))
    plt.subplot(3, 1, 1)
    plt.title('Read Depth')
    plt.plot(x, depth)

    plt.subplot(3, 1, 2)
    plt.title('First Order Deriviate')
    plt.plot(range(start, end+1), depth_diffs)
    plt.hlines(0, start, end)

    plt.subplot(3, 1, 3)
    plt.title('Second Order Deriviate')
    plt.plot(range(start, end+1), res2)
    plt.hlines(0, start, end)

    peaks = []
    for idx in indexes:
        if idx>=100 and idx+100<=len(depth) and depth_window_slice[idx]>=min_depth and res2[idx]>0:
            position = range(start, end + 1)[idx]
            mean_left_50_100 = np.mean(depth[max(0, idx-100):max(0, idx-50)])
            mean_right_50_100 = np.mean(depth[min(len(depth), idx+50):min(len(depth), idx+100)])
            mean_middle_100 = np.mean(depth[idx-50:idx+50])
            plt.vlines(position, min(res2), max(res2))
            peaks.append([gene, position, depth_window_slice[idx], mean_middle_100, mean_left_50_100, mean_right_50_100])

    peaks = [x + [len(peaks)] for x in peaks]
    if len(peaks)>=2:
        plt.savefig("Depth_{}.png".format(gene))
        return peaks

    else:
        return []

def scan_genome(genome, bam):
    """
    Report ALL THE APA CONDICATE...
    GENE/START/END/PA_START/PA_END/STRAND/MEANDEPTH/
    """
    from baseq.bam import BAMTYPE
    from baseq.drops.apa.scaner import scan
    from baseq.rna.gtf.gencode import read_gencode
    # Read Gencode And Get UTR Region (>1000bp)
    # Sort by UTR length and get the logest by gene...
    df = read_gencode(genome, "UTR")
    df['length'] = df.end - df.start
    df = df.sort_values(by=['length'])
    df = df.groupby("gene").last()
    df = df.loc[df.length>1000]
    bam = BAMTYPE(bam)

    import multiprocessing as mp
    pool = mp.Pool(processes=20)
    results = []
    for index, row in df.iterrows():
        results.append(pool.apply_async(scan, (bam, index, row['chr'], row['start'], row['end'], 50)))
    pool.close()
    pool.join()

    results = [x.get() for x in results]
    results = [y for x in results for y in x]
    results = [df.loc[x[0],:].tolist() + x  for x in results]
    df_peaks = pd.DataFrame(results, columns=["chr", "start", 'end', 'strand', 'transc', 'exon', 'length', 'gene', 'pos', 'mean_depth', "mid", "left", "right", "counts"])
    df_peaks = df_peaks.drop(columns=["transc", 'exon'])

    df_peaks.to_csv("APA.txt", sep="\t")
    df_peaks.to_excel("APA.xls")