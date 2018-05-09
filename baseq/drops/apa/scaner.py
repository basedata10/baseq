import os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import peakutils

def scan(depth, start, end, mindepth=10):
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

    plt.figure(figsize=(10, 5))
    plt.subplot(3, 1, 1)
    plt.title('Read Depth')
    plt.plot(x, depth)

    plt.subplot(3, 1, 2)
    plt.title('Read Depth')
    plt.plot(range(start, end+1), depth_diffs)
    plt.hlines(0, start, end)

    plt.subplot(3, 1, 3)
    plt.title('Read Depth')
    plt.plot(range(start, end+1), res2)
    plt.hlines(0, start, end)

    peaks = []
    for idx in indexes:
        if depth_window_slice[idx]>=mindepth and res2[idx]>0:
            plt.vlines(range(start, end + 1)[idx], min(res2), max(res2))
            peaks.append(range(start, end + 1)[idx])
    plt.savefig("Depth.png")
    return peaks