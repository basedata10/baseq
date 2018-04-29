
sample={{name}}
fq1={{fq1}}
fq2={{fq2}}
intervals={{interval}}
outdir={{outdir}}

cd $outdir

##Map to Reference
$bwa mem -t 20 -M -R "@RG\tID:${sample}\tSM:${sample}\tLB:WES\tPL:Illumina" $bwa_index $fq1 $fq2 1>${sample}.sam
$java -Djava.io.tmpdir=$temp_dir -Xmx40g -jar $picard SortSam SORT_ORDER=coordinate INPUT=${sample}.sam OUTPUT=${sample}.bam
$samtools index ${sample}.bam

##Mark Duplicates
$java -Djava.io.tmpdir=$temp_dir -Xmx40g -jar $picard MarkDuplicates INPUT=${sample}.bam OUTPUT=${sample}_marked.bam METRICS_FILE=${sample}.metrics
$samtools index ${sample}_marked.bam

##GATK BQSR
$GATK BaseRecalibrator -R $genome -L $intervals -I ${sample}_marked.bam --known-sites $dbSNP --known-sites $SNP --known-sites $INDEL -O ${sample}_temp.table
$GATK ApplyBQSR -R $genome -L $intervals -I ${sample}_marked.bam -bqsr ${sample}_temp.table -O ${sample}_marked_bqsr.bam

##Haplotype Caller
$GATK --java-options "-Xmx4g" HaplotypeCaller -R $genome -I ${sample}_marked_bqsr.bam -O ${sample}_raw.snps.indels.vcf -bamout bamout.bam -L $intervals --native-pair-hmm-threads 20

##SelectVariants
$GATK SelectVariants -R $genome -V ${sample}_raw.snps.indels.vcf --select-type-to-include SNP -O ${sample}_raw_snps.vcf
$GATK VariantFiltration -R $genome -V ${sample}_raw_snps.vcf -O ${sample}_filtered_snps.vcf --filter-expression "QD < 2 || FS > 60 || MQ < 40" --filter-name "snp_filter"

##Remove Temp Files
rm ${sample}.sam
rm ${sample}.bam*
rm ${sample}_marked.bam*
rm ${sample}_recal.ba*