docs ="""

[Scan Region]
baseq-Drop apa_scan_region -b star.CC_Aligned.sort.bam -c chr1 -s 12141952 -e 12145429
baseq-Drop apa_scan_region -b star.CC_Aligned.sort.bam -c chr1 -s 12204917 -e 12211873
baseq-Drop apa_scan_region -b Merged.sort.A.bam -c chr8 -s 102648432 -e 102650458
baseq-Drop apa_scan_region -b Merged.sort.A.bam -c chr8 -s 102251429 -e 102255483

[Scan Whole Genome]
# Generate a APA.txt with detected PA sites...
baseq-Drop apa_scan_genome -g hg38 -b Merged.sort.A.bam

[Filter Genes With Multilpe APA Site]
#generate Biomodule.genes.txt which 
#gene	left	right	depth_l	depth_r
/mnt/gpfs/Users/zhangxiannian/basematic/dropseq/APA

[Sample Genes APA usage]: generate a table of UMIs for each APA...
baseq-Drop apa_sample_usage -b Merged.sort.A.bam --apalist APA.txt -g POU2F2
baseq-Drop apa_sample_usage -b Merged.A.bam --apalist APA.txt -g SAMHD1


"""