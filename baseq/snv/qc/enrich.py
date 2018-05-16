
def enrich_qc(samplename, bampath, intervals):
    """
    Check the coverage depth and enrichment quality.

    Usage:
    ::
        enrich_qc("sample01", "xx.bam", "panel.bed")

    Return:
        Sample/Total/Mapped/Map_Ratio/Dup_ratio/PCT_10X/PCT_30X/...
    """

    from baseq.bam.bamtype import BAMTYPE

    bam = BAMTYPE(bampath, bedfile=intervals)
    bam.stats_bases()
    bam.stats_duplicates()
    bam.stats_regions()
    bam.stats_region_coverage(1000)

    stats = {
        "Sample" : samplename,
        "Total" : bam.reads_total,
        "Mapped" : bam.reads_mapped,
        "Map_Ratio" : bam.mapping_ratio,
        "Dup_ratio" : bam.dup_ratio,
        "Mean_Depth": bam.mean_depth,
        "PCT_10X" :  bam.pct_10X,
        "PCT_30X": bam.pct_30X,
        "PCT_50X": bam.pct_50X,
        "PCT_100X": bam.pct_100X,
    }

    return stats