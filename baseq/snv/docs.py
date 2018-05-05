"""
#Enrich Quality
baseq-SNV qc_enrich ./bam ./bedfile ./out


#Alignment
baseq-SNV run_bwa -1 Read1M.P457.1.fq.gz -2 Read1M.P457.2.fq.gz -g hg38 -n Test -o Test.bam

#MarkDuplicate
baseq-SNV run_markdup -b Test.bam -m Test.marked.bam



"""