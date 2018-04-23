#!/bin/bash
export PATH=/mnt/gpfs/Database/softs/anaconda2/bin:$PATH
TMPDIR=/tmp
GENOME=/mnt/gpfs/Database/ref/hg38/hg38.fa
GATK=/mnt/gpfs/Users/wufan/p12_HEC/GATK/GATK_code/gatk-4.0.3.0/gatk
PICARD=/mnt/gpfs/Users/wufan/p12_HEC/GATK/picard.jar
DBSNP=/mnt/gpfs/Users/wufan/p12_HEC/GATK/resources/dbsnp_138.hg38.vcf.gz
path_annovar=/mnt/gpfs/Users/wufan/p12_HEC/GATK/annovar_new/annovar
ref_annovar=/mnt/gpfs/Users/wufan/p12_HEC/GATK/annovar/human38db

sample={{sample}}
intervals={{interval}}
fq1={{fq1}}
fq2={{fq2}}
logfile=./log.txt

##Make Work Dir
mkdir ./$sample
cd ./$sample

##Check File1 Exists
#basematic checkfq $fq1 $fq2

##Bwa alignment
bwa mem -t 20 -M -R "@RG\tID:$sample\tSM:$sample\tLB:WES\tPL:Illumina" $GENOME $fq1 $fq2 1>$sample.sam
samtools flagstat $sample.sam >flagstat.txt

##Picard Sortsam
java -Djava.io.tmpdir=$TMPDIR -Xmx40g -jar $PICARD SortSam SORT_ORDER=coordinate INPUT=$sample.sam OUTPUT=$sample.bam
samtools index $sample.bam

##Picard Markduplicates
java -Djava.io.tmpdir=$TMPDIR -Xmx40g -jar $PICARD MarkDuplicates INPUT=$sample.bam OUTPUT=${sample}_marked.bam METRICS_FILE=$sample.metrics

##Picard FixMateInformation
java -Djava.io.tmpdir=$TMPDIR -Xmx40g -jar $PICARD FixMateInformation INPUT=${sample}_marked.bam OUTPUT=${sample}_marked_fixed.bam SO=coordinate
samtools index ${sample}_marked_fixed.bam

##GATK SplitNCigarReads
$GATK SplitNCigarReads -R $GENOME -L $intervals -I ${sample}_marked_fixed.bam -O ${sample}_marked_fixed_split.bam

##GATK BaseRecalibrator
$GATK BaseRecalibrator -R $GENOME -L $intervals -I ${sample}_marked_fixed_split.bam -O ${sample}_temp.table --known-sites $DBSNP

##GATK HaplotypeCaller
$GATK --java-options "-Xmx4g" HaplotypeCaller -R $GENOME -I ${sample}_recal.bam -O ${sample}_raw.snps.indels.vcf -bamout ${sample}_Haplotype.bam -L $intervals --native-pair-hmm-threads 20
$GATK SelectVariants -R $GENOME -V ${sample}_raw.snps.indels.vcf --select-type-to-include SNP -O ${sample}_raw_snps.vcf

##Annovar
$path_annovar/convert2annovar.pl -format vcf4 ${sample}_raw_snps.vcf > ${sample}_snps.avinput
$path_annovar/table_annovar.pl ${sample}_snps.avinput $ref_annovar -buildver hg38 -out ${sample} -remove -protocol refGene,genomicSuperDups,esp6500siv2_all,1000G2015aug_ALL,1000G2015aug_EAS,exac03,avsnp147,dbnsfp30a,clinvar_20170130,cosmic70,dbscsnv11,cytoBand -operation g,r,f,f,f,f,f,f,f,f,f,r -nastring . -csvout
#basematic-SNV annovar --intype gatk --out ......

##Filter VCFs
#basematic-snp vcffilter --from annovar --type csv --file ${sample}.hg38_multianno.csv --out annovarfiles.txt

##Stats Result
#basematic-bam stats --file ${sample}.bam --out annovarfiles.txt

##Build Result
#basematic visualize --logfile ./logs.txt --out result.txt

##Remove Temp Files
rm ${sample}.sam
rm ${sample}.bam*
rm ${sample}_marked.bam*
rm ${sample}_marked_fixed.bam*
rm ${sample}_marked_fixed_split*
rm ${sample}_Haplotype.bam
rm ${sample}_recal.b*