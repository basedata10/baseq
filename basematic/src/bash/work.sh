#!/usr/bin/env bash

samtools=/mnt/gpfs/Database/softs/anaconda2/bin/samtools
java=/mnt/gpfs/Database/softs/anaconda2/bin/java
annovar=/mnt/gpfs/Users/wufan/p12_HEC/GATK/annovar_new/annovar
picard=/mnt/gpfs/Users/wufan/p12_HEC/GATK/picard.jar
hg38_bwa_index=/mnt/gpfs/Database/ref/hg38/hg38.fa
bwa=/mnt/gpfs/Database/softs/anaconda2/bin/bwa
genome=/mnt/gpfs/Database/ref/hg38/hg38.fa
temp_dir=/tmp
gatk=/mnt/gpfs/Users/wufan/p12_HEC/GATK/GATK_code/gatk-4.0.3.0/gatk
dbsnp=/mnt/gpfs/Users/wufan/p12_HEC/GATK/resources/dbsnp_138.hg38.vcf.gz
annovar_db_hg38=/mnt/gpfs/Users/wufan/p12_HEC/GATK/annovar/human38db

sample=S001
fq1=/mnt/gpfs8/Users/zhangxiannian/projects/p15_panel/01_vcf_test/Read1M.P457.1.fq.gz
fq2=/mnt/gpfs8/Users/zhangxiannian/projects/p15_panel/01_vcf_test/Read1M.P457.2.fq.gz
interval=pad_probe.intervals
outdir=/mnt/gpfs8/Users/zhangxiannian/projects/p15_panel/tests2/S001

cd ${outdir}

#Map to Reference
$bwa mem -t 20 -M -R "@RG\tID:$sample\tSM:$sample\tLB:WES\tPL:Illumina" $hg38_bwa_index $fq1 $fq2 1>${sample}.sam
$java -Djava.io.tmpdir=$temp_dir -Xmx40g -jar $picard SortSam SORT_ORDER=coordinate INPUT=${sample}.sam OUTPUT=${sample}.bam
$samtools index ${sample}.bam

#Mark Duplicates
$java -Djava.io.tmpdir=$temp_dir -Xmx40g -jar $picard MarkDuplicates INPUT=${sample}.bam OUTPUT=${sample}_marked.bam METRICS_FILE=${sample}.metrics
$samtools index ${sample}_marked.bam

#BQSR
$GATK BaseRecalibrator -R $genome -L $intervals -I ${sample}_marked.bam --known-sites $dbSNP --known-sites $SNP --known-sites $INDEL -O ${sample}_temp.table
$GATK ApplyBQSR -R $genome -L $intervals -I ${sample}_marked.bam -bqsr ${sample}_temp.table -O ${sample}_marked_bqsr.bam

#validate bam file
$java -Djava.io.tmpdir=$temp_dir -Xmx40g -jar $picard Validat

#Haplotype Caller
$GATK --java-options "-Xmx4g" HaplotypeCaller -R $genome -I ${sample}_marked_bqsr.bam -O ${sample}_raw.snps.indels.vcf -bamout bamout.bam -L $intervals --native-pair-hmm-threads 20

#SelectVariants
$GATK SelectVariants -R $genome -V ${sample}_raw.snps.indels.vcf --select-type-to-include SNP -O ${sample}_raw_snps.vcf
$GATK VariantFiltration -R $genome -V ${sample}_raw_snps.vcf -O ${sample}_filtered_snps.vcf --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" --filter-name "my_snp_filter"

#annotation
$path_annovar/convert2annovar.pl -format vcf4 ${sample}_raw_snps.vcf > ${sample}_snps.avinput
$path_annovar/table_annovar.pl ${sample}_snps.avinput $ref_annovar -buildver hg38 -out GATK_test -remove -protocol \
refGene,esp6500siv2_all,1000G2015aug_ALL,1000G2015aug_EAS,exac03,avsnp147,dbnsfp30a,clinvar_20170130,cosmic70,dbscsnv11,cytoBand -operation g,f,f,f,f,f,f,f,f,f,r -nastring . -csvout