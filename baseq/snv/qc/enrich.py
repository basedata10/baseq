import os, sys, time
from baseq.bam import BAMFILE
from baseq.bed import BEDFILE
from multiprocessing.pool import ThreadPool

def mean_depth_in_bedfiles(bam, bed, number=500):
    intervals = bed.sampling(number)
    pool = ThreadPool(processes=10)
    results = []
    for interval in intervals:
        results.append(pool.apply_async(bam.bin_mean_depth, (interval[0][3:], interval[1], interval[2])))
    pool.close()
    pool.join()
    region_depth = [r.get() for r in results]
    return region_depth

def enrich_quality(bampath, bedfile, outpath):
    bam = BAMFILE(bampath)
    bed = BEDFILE(bedfile)

    #get bases in the region
    region_mean_depth = mean_depth_in_bedfiles(bam, bed)
    region_size = bed.length
    region_bases = region_mean_depth * region_size

    #get bases in the whole genome
    bam.bam_stats()
    bam.match_length()
    total_bases = bam.mapped_reads * bam.match_length

    #10X 30X 50X 100X ratio
    
