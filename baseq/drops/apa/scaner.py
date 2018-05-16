import os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import peakutils
import pandas as pd

def scan(bam, name, chr, start, end, min_depth=10):
    """ SCAN THE GENOME...

    """
    start = start-100
    end = end + 100
    #Get Depth
    depth = bam.region_depth(chr, start, end, all=True)
    print(depth)
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

    plt.figure(figsize=(100, 7))
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
            peaks.append([
                name, position, depth_window_slice[idx],
                mean_middle_100, mean_left_50_100, mean_right_50_100]
            )

    #Get Strand Information...
    #For each Peaks Get the reads in the 50upstream and 50 downstream
    peaks = [x + [len(peaks)] for x in peaks]
    print("[info] {}:Depth_{}.png".format(name, name))
    if len(peaks)>=2:
        plt.savefig("Depth_{}.png".format(name))
        return peaks
    else:
        return []