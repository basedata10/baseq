import os
from baseq.setting import bash_script_dir
from baseq.mgt import get_config, run_cmd


PICARD = get_config("SNV", "picard")
GATK = get_config("SNV", "GATK")



bwa_cmd_script_p = r"""{0} mem -t 8 -M -R "@RG\tID:{1}\tSM:{1}\tLB:WES\tPL:Illumina" {2} {3} {4} 1>{1}.sam"""
bwa_cmd_script_s = r"""{0} mem -t 8 -M -R "@RG\tID:{1}\tSM:{1}\tLB:WES\tPL:Illumina" {2} {3} 1>{1}.sam"""
sort_index_cmd_script = """
{3} -jar {0} SortSam SORT_ORDER=coordinate INPUT={1}.sam OUTPUT={1}.bam
{2} index {1}.bam
"""
def run_alignment(fq1,fq2,sample,genome,run=True):
    bwa = get_config("SNV", "bwa")
    PICARD = get_config("SNV", "picard")
    java = get_config("SNV", "java")
    samtools = get_config("SNV", "samtools")
    genome = get_config("SNV_ref_"+genome, "bwa_index")
    if fq1 and fq2 and os.path.exists(fq1) and os.path.exists(fq2):
        bwa_cmd = bwa_cmd_script_p.format(bwa,sample,genome,fq1,fq2)
    elif fq1 and os.path.exists(fq1):
        bwa_cmd = bwa_cmd_script_s.format(bwa,sample,genome,fq1)
    sort_index_cmd=sort_index_cmd_script.format(PICARD,sample,samtools,java)
    if run:
        run_cmd("bwa alignment","".join(bwa_cmd))
        run_cmd("PICARD SortSam","".join(sort_index_cmd))
    return bwa_cmd+"\n"+sort_index_cmd



markdup_cmd_script ="""
{0} -jar {1} MarkDuplicates INPUT={2}.bam OUTPUT={2}_marked.bam METRICS_FILE={2}.metrics
{3} index {2}_marked.bam
"""
def run_markdup(sample,run=True):
    java = get_config("SNV", "java")
    PICARD = get_config("SNV", "picard")
    samtools = get_config("RNA", "samtools")
    markdup_cmd = markdup_cmd_script.format(java,PICARD,sample,samtools)
    if run:
        run_cmd("Mark duplicates","".join(markdup_cmd))
    return markdup_cmd

bqsr_cmd_script = """
{0} BaseRecalibrator -R {1} -L {2} -I {3}_marked.bam --known-sites {4} --known-sites {5} --known-sites {6} -O {3}_temp.table
{0} ApplyBQSR -R {1} -L {2} -I {3}_marked.bam -bqsr {3}_temp.table -O {3}_marked_bqsr.bam
"""
def bqsr(sample,genome,run=True):
    GATK = get_config("SNV", "GATK")
    index = get_config("SNV_ref_"+genome,"bwa_index")
    DBSNP = get_config("SNV_ref_"+genome,"DBSNP")
    SNP = get_config("SNV_ref_"+genome,"SNP")
    INDEL = get_config("SNV_ref_"+genome,"INDEL")
    interval = get_config("SNV_ref_"+genome,"interval")
    bqsr_cmd = bqsr_cmd_script.format(GATK,index,interval,sample,DBSNP,SNP,INDEL)
    if run:
        run_cmd("BaseRecalibrator","".join(bqsr_cmd))
    return bqsr_cmd

callvar_cmd_script = """
{0} --java-options "-Xmx4g" HaplotypeCaller -R {1} -L {2} -I {3}_marked_bqsr.bam -O {3}_raw.snps.indels.vcf -bamout bamout.bam --native-pair-hmm-threads 20
"""
def run_callvar(sample,genome,run=True):
    GATK = get_config("SNV", "GATK")
    index = get_config("SNV_ref_" + genome, "bwa_index")
    interval = get_config("SNV_ref_" + genome, "interval")
    callvar_cmd = callvar_cmd_script.format(GATK, index, interval, sample)
    if run:
        run_cmd("call variants","".join(callvar_cmd))
    return callvar_cmd

selectvar_cmd_script = """
{0} SelectVariants -R {1} -V {2}_raw.snps.indels.vcf --select-type-to-include SNP -O {2}_raw_snps.vcf
{0} VariantFiltration -R {1} -V {2}_raw_snps.vcf -O {2}_filtered_snps.vcf --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" --filter-name "my_snp_filter"
"""
def selectvar(sample,genome,run=True):
    GATK = get_config("SNV", "GATK")
    index = get_config("SNV_ref_" + genome, "bwa_index")
    selectvar_cmd = selectvar_cmd_script.format(GATK,index,sample)
    if run:
        run_cmd("SelectVariants","".join(selectvar_cmd))
    return selectvar_cmd


