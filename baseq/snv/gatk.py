import os
from baseq.setting import bash_script_dir
from baseq.mgt import get_config, run_cmd

PICARD = get_config("SNV", "picard")
GATK = get_config("SNV", "GATK")

bwa_cmd_script_p = r"""{bwa} mem -t 8 -M -R "@RG\tID:{sample}\tSM:{sample}\tLB:WES\tPL:Illumina" {genome} {fq1} {fq2}  1>{samfile}"""
bwa_cmd_script_s = r"""{bwa} mem -t 8 -M -R "@RG\tID:{sample}\tSM:{sample}\tLB:WES\tPL:Illumina" {genome} {fq1} 1>{samfile}"""
sort_index_cmd_script = """
{java} -jar {picard} SortSam SORT_ORDER=coordinate INPUT={samfile} OUTPUT={outfile}
{samtools} index {outfile}
rm {samfile}
"""
def run_alignment(fq1, fq2, sample, genome, outfile):
    bwa = get_config("SNV", "bwa")
    picard = get_config("SNV", "picard")
    java = get_config("SNV", "java")
    samtools = get_config("SNV", "samtools")
    genome = get_config("SNV_ref_"+genome, "bwa_index")
    samfile = outfile+".sam"

    if fq1 and fq2 and os.path.exists(fq1) and os.path.exists(fq2):
        bwa_cmd = bwa_cmd_script_p.format(bwa=bwa, sample=sample, genome=genome, fq1=fq1, fq2=fq2, samfile=samfile)

    elif fq1 and os.path.exists(fq1):
        bwa_cmd = bwa_cmd_script_s.format(bwa=bwa, sample=sample, genome=genome, fq1=fq1, samfile=samfile)

    sort_index_cmd=sort_index_cmd_script.format(picard=picard, sample=sample, samtools=samtools, java=java, outfile=outfile, samfile=samfile)

    run_cmd("bwa alignment", "".join(bwa_cmd))
    run_cmd("PICARD SortSam", "".join(sort_index_cmd))

    return bwa_cmd+"\n"+sort_index_cmd

markdup_cmd_script ="""
{java} -jar {picard} MarkDuplicates INPUT={bamfile} OUTPUT={markedbam} METRICS_FILE={markedbam}.metrics
{samtools} index {markedbam}
"""

def run_markdup(bamfile, markedbam):
    java = get_config("SNV", "java")
    picard = get_config("SNV", "picard")
    samtools = get_config("RNA", "samtools")
    cmd = markdup_cmd_script.format(java=java, picard=picard, samtools=samtools, markedbam=markedbam, bamfile=bamfile)
    run_cmd("Mark duplicates","".join(cmd))
    return cmd

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