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

#create PoN files
baseq-SNV create_pon -p /mnt/gpfs/Users/wufan/p12_HEC/GATK/baseq_mutect_test/ -l /mnt/gpfs/Users/wufan/p12_HEC/GATK/28wuchen/mutect_call.txt -L /mnt/gpfs/Users/wufan/p12_HEC/GATK/28wuchen/merge.all.target.1.list -g hg37

#single mutect:/mnt/gpfs/Users/wufan/p12_HEC/GATK/28wuchen/N506/
baseq-SNV run_mutect2 -g hg37 -n N506 -N N506_marked_bqsr.bam -t T506 -T T506_marked_bqsr.bam -o ./

#filter mutect call
baseq-SNV filter_mutect_call -r /mnt/gpfs/Users/wufan/p12_HEC/GATK/resources/ref_b37/small_exac_common_3_b37.vcf.gz -s /mnt/gpfs/Users/wufan/p12_HEC/GATK/mutect_test/T506.vcf.gz -T /mnt/gpfs/Users/wufan/p12_HEC/GATK/28wuchen/T506/T506_marked_bqsr.bam -o ./ -t T506
"""