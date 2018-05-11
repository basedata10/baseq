docs = """
#stats multiple stuffs about a bam file...
baseq-BAM infos -i S6.final.bam -b ../refs/merged.all.target.bed -c 1000
baseq-BAM infos -i /mnt/gpfs/Users/wufan/p12_HEC/GATK/12wes/P055_M1/P055_M1_marked_bqsr.bam -b /mnt/gpfs/Users/wufan/p12_HEC/GATK/resources/ref_b37/Broad.human.exome.b37.interval_list -c 1000

[sort]
baseq-BAM sort_bam AA.sort.bam sorted.AA.bam

[region=>reads] the first 5 columns of a bamfile
baseq-BAM region_reads -b Merged.sort.A.bam -p chr2:74,155,931-74,156,000


"""