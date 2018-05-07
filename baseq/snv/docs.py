doc="""
#Enrich Quality
baseq-SNV qc_enrich ./bam ./bedfile ./out


#Alignment
baseq-SNV run_bwa -1 Read1M.P457.1.fq.gz -2 Read1M.P457.2.fq.gz -g hg38 -n Test -o Test.bam -t 10

#MarkDuplicate
baseq-SNV run_markdup -b Test.bam -m Test.marked.bam

#bqsr
baseq-SNV run_bqsr -m Test.marked.bam -g hg38 -q Test.marked.bqsr.bam

#call variants
baseq-SNV run_callvar -q Test.marked.bqsr.bam -r Test.raw.indel.snp.vcf -g hg38

#select variants
baseq-SNV run_selectvar -r Test.raw.indel.snp.vcf -s Test.raw.snp.vcf -f Test.filtered.snp.vcf -g hg38

#annovar annotation
baseq-SNV run_annovar -g hg38 -n Test -f Test.filtered.snp.vcf -a Test.snps.avinput

#run gatk pipeline
baseq-SNV run_gatkpipe -1 Read1M.P457.1.fq.gz -2 Read1M.P457.2.fq.gz -n Test -g hg38 -d ./ 

#run gatk pipeline from bam file
baseq-SNV run_gatkpipe -m 
"""