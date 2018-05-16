import os
from baseq.mgt import get_config, run_cmd

# PICARD = get_config("SNV", "picard")
# GATK = get_config("SNV", "GATK")

bwa_cmd_script_p = r"""{bwa} mem -t {thread} -M -R "@RG\tID:{sample}\tSM:{sample}\tLB:WES\tPL:Illumina" {genome} {fq1} {fq2}  1>{samfile}"""
bwa_cmd_script_s = r"""{bwa} mem -t {thread} -M -R "@RG\tID:{sample}\tSM:{sample}\tLB:WES\tPL:Illumina" {genome} {fq1} 1>{samfile}"""
sort_index_cmd_script = """
{samtools} view -b -u -S {samfile}>{viewedbam}
{samtools} sort -@ 8 {viewedbam} -o {outfile}
{samtools} index {outfile}
rm {samfile}
"""

def alignment(fq1, fq2, sample, genome, outfile, thread=8):
    """
    Map fastq1/2 files into genome using BWA. Add tags to bamfile using the input sample name. The bamfile is named as outfile.
    """
    bwa = get_config("SNV", "bwa")
    samtools = get_config("SNV", "samtools")
    genome = get_config("SNV_ref_"+genome, "bwa_index")
    viewedbam = outfile + ".view.bam"
    samfile = outfile+".sam"
    if fq1 and fq2 and os.path.exists(fq1) and os.path.exists(fq2):
        bwa_cmd = bwa_cmd_script_p.format(bwa=bwa, sample=sample, genome=genome, fq1=fq1, fq2=fq2, samfile=samfile, thread=thread)
    elif fq1 and os.path.exists(fq1):
        bwa_cmd = bwa_cmd_script_s.format(bwa=bwa, sample=sample, genome=genome, fq1=fq1, samfile=samfile, thread=thread)
    sort_index_cmd=sort_index_cmd_script.format(samtools=samtools, outfile=outfile, samfile=samfile,viewedbam=viewedbam)
    run_cmd("bwa alignment", "".join(bwa_cmd))
    run_cmd("samtools sort", "".join(sort_index_cmd))
    return bwa_cmd+"\n"+sort_index_cmd

markdup_cmd_script ="""
{java} -jar {picard} MarkDuplicates INPUT={bamfile} OUTPUT={markedbam} METRICS_FILE={markedbam}.metrics
{samtools} index {markedbam}
"""

def run_markdup(bamfile, markedbam):
    """
    Run MarkDuplicate of Picard. Generate the bai for the marked bamfile.
    ::
        run_markdup("in.bam", "in.marked.bam")
    """
    java = get_config("SNV", "java")
    picard = get_config("SNV", "picard")
    samtools = get_config("RNA", "samtools")
    cmd = markdup_cmd_script.format(java=java, picard=picard, samtools=samtools, markedbam=markedbam, bamfile=bamfile)
    run_cmd("Mark duplicates","".join(cmd))
    return cmd

bqsr_cmd_script = """
{gatk} BaseRecalibrator -R {index} -L {interval} -I {markedbam} --known-sites {dbsnp} --known-sites {snp} --known-sites {indel} -O {markedbam}.table
{gatk} ApplyBQSR -R {index} -L {interval} -I {markedbam} -bqsr {markedbam}.table -O {bqsrbam}
"""
bqsr_cmd_script_DRF = """
{gatk} BaseRecalibrator -R {index} -L {interval} -I {markedbam} --known-sites {dbsnp} --known-sites {snp} --known-sites {indel} --disable-read-filter NotDuplicateReadFilter -O {markedbam}.table
{gatk} ApplyBQSR -R {index} -L {interval} -I {markedbam} -bqsr {markedbam}.table --disable-read-filter NotDuplicateReadFilter -O {bqsrbam}
"""

def bqsr(markedbam, bqsrbam, genome, disable_dup_filter=False):
    """
    Run BQSR_.
    ::
        bqsr()
        This will generate a XXXX...

    .. _BQSR: https://gatkforums.broadinstitute.org/gatk/discussion/44/base-quality-score-recalibration-bqsr
    """
    gatk = get_config("SNV", "GATK")
    index = get_config("SNV_ref_"+genome,"bwa_index")
    DBSNP = get_config("SNV_ref_"+genome,"DBSNP")
    SNP = get_config("SNV_ref_"+genome,"SNP")
    INDEL = get_config("SNV_ref_"+genome,"INDEL")
    interval = get_config("SNV_ref_"+genome,"interval")

    if not disable_dup_filter:
        bqsr_cmd = bqsr_cmd_script.format(gatk=gatk, index=index, interval=interval, markedbam=markedbam,
                                          bqsrbam=bqsrbam, dbsnp=DBSNP, snp=SNP, indel=INDEL)
    else:
        bqsr_cmd = bqsr_cmd_script_DRF.format(gatk=gatk, index=index, interval=interval, markedbam=markedbam,
                                          bqsrbam=bqsrbam, dbsnp=DBSNP, snp=SNP, indel=INDEL)
    run_cmd("BaseRecalibrator","".join(bqsr_cmd))

callvar_cmd_script = """
{gatk} --java-options "-Xmx4g" HaplotypeCaller -R {index} -L {interval} -I {bqsrbam} -O {rawvcf} -bamout bamout.bam --native-pair-hmm-threads 20 {extrainfos}
"""

def run_callvar(bqsrbam, rawvcf, genome, disable_dup_filter = False):
    GATK = get_config("SNV", "GATK")
    index = get_config("SNV_ref_" + genome, "bwa_index")
    interval = get_config("SNV_ref_" + genome, "interval")
    extra = ""
    if disable_dup_filter:
        extra = "--disable-read-filter NotDuplicateReadFilter"
    callvar_cmd = callvar_cmd_script.format(gatk=GATK, index=index, interval=interval, bqsrbam=bqsrbam, rawvcf=rawvcf, extrainfos=extra)
    run_cmd("call variants","".join(callvar_cmd))
    return callvar_cmd

selectvar_cmd_script = """
{gatk} SelectVariants -R {index} -V {rawvcf} --select-type-to-include SNP -O {selectvcf}
{gatk} VariantFiltration -R {index} -V {selectvcf} -O {filtervcf} --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" --filter-name "my_snp_filter"
"""

def selectvar(rawvcf,selectvcf,filtervcf,genome,run=True):
    """
    Select Variants, it process the XX from XX, genrate XX for XX...
    """
    GATK = get_config("SNV", "GATK")
    index = get_config("SNV_ref_" + genome, "bwa_index")
    selectvar_cmd = selectvar_cmd_script.format(gatk=GATK,index=index,rawvcf=rawvcf,selectvcf=selectvcf,filtervcf=filtervcf)
    if run:
        run_cmd("SelectVariants","".join(selectvar_cmd))
    return selectvar_cmd

mutect2_cmd_simplified_script="""
{gatk} Mutect2 -R {index} -I {normalbam} -normal {normalname} -I {tumorbam} -tumor {tumorname} -O {vcffile}
"""
mutect2_cmd_stardand_script="""
{gatk} Mutect2 -R {index} -I {normalbam} -normal {normalname} -I {tumorbam} -tumor {tumorname} --germline-resource {germline} --panel-of-normals {pon} -O {vcffile}
"""
def mutect2(genome, normalname, normalbam, tumorname, tumorbam, vcffile, pon, germline):
    gatk = get_config("SNV","GATK")
    index = get_config("SNV_ref_" + genome, "bwa_index")
    if pon :
        if not germline:
            germline = get_config("SNV_ref_" + genome, "germline")
            mutect2_cmd = mutect2_cmd_stardand_script.format(gatk=gatk, index=index, normalbam=normalbam,
                                                             normalname=normalname, tumorbam=tumorbam,
                                                             tumorname=tumorname, vcffile=vcffile, germline=germline,
                                                             pon=pon)
        else:
            mutect2_cmd = mutect2_cmd_stardand_script.format(gatk=gatk, index=index, normalbam=normalbam,
                                                             normalname=normalname, tumorbam=tumorbam,
                                                             tumorname=tumorname, vcffile=vcffile, germline=germline,
                                                             pon=pon)
    else:
        mutect2_cmd = mutect2_cmd_simplified_script.format(gatk=gatk, index=index, normalbam=normalbam,
                                                           normalname=normalname, tumorbam=tumorbam,
                                                           tumorname=tumorname, vcffile=vcffile)


    run_cmd("mutect annlysis", "".join(mutect2_cmd))


filtertable_cmd_script = """
{gatk} GetPileupSummaries -I {tumorbam} -V {resource} -O {gps_table}
{gatk} CalculateContamination -I {gps_table} -O {calcontam_table}
"""
def get_filter_table(tumorbam, resource, gps_table, calcontam_table):
    gatk = get_config("SNV","GATK")
    filtertable_cmd = filtertable_cmd_script.format(gatk=gatk, tumorbam=tumorbam, resource=resource, gps_table=gps_table,
                                                    calcontam_table=calcontam_table)
    run_cmd("obatin filter table for mutect calls","".join(filtertable_cmd))



filtercall_cmd_script = """
{gatk} FilterMutectCalls -V {somaticvcf} --contamination-table {calcontam_table} -O {filter_call}
"""
def filter_mutect_vcf(somaticvcf,calcontam_table,filter_call):
    gatk = get_config("SNV","GATK")
    filtercall_cmd = filtercall_cmd_script.format(gatk=gatk, somaticvcf=somaticvcf, calcontam_table=calcontam_table,
                                                    filter_call=filter_call)
    run_cmd("filter mutect calls using contamination table","".join(filtercall_cmd))


def listofvcf(path):
    normalargs=open(os.path.join(path,"normalvcf_for_pon.args"),"w")
    filelist=[]
    filename = os.listdir(path)
    for fn in filename:
        if fn[-6:] == "vcf.gz":
           fullfilename = os.path.join(path,fn)
           filelist.append(fullfilename)
    normalargs.write("\n".join([list for list in filelist]) + "\n")
    return os.path.join(path,"normalvcf_for_pon.args")

normalvcf_cmd_script="""
{gatk} Mutect2 -R {index} -I {normalbam} -tumor {samplename} --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter -L {interval} -O {normalvcf}
"""
normalvcf_cmd_script1="""
{gatk} Mutect2 -R {index} -I {normalbam} -tumor {samplename} --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter -O {normalvcf}
"""
pon_cmd_script="""
{gatk} CreateSomaticPanelOfNormals -vcfs {normalvcfs} -O {ponvcf}
"""
def create_pon(genome,list,path,interval):
    index = get_config("SNV_ref_" + genome, "bwa_index")
    gatk = get_config("SNV", "GATK")
    if not os.path.exists(path):
        print("[ERROR] No such file or directory")
    else:
        path_pon = os.path.join(path, "pon")
        if not os.path.exists(path_pon):
            os.mkdir(path_pon)

    with open(list, "r") as file:
        lines = file.readlines()
    sample_info = [line.strip().split() for line in lines]
    import multiprocessing as mp
    pool = mp.Pool(processes=6)
    results = []
    for sample in sample_info:
        normalvcf = os.path.join(path_pon, "{}_tumor-only.vcf.gz".format(sample[0]))
        if interval:
            normalvcf_cmd = normalvcf_cmd_script.format(gatk=gatk, index=index, normalbam=sample[1],
                                                        samplename=sample[0],
                                                        normalvcf=normalvcf, interval=interval)
        else:
            normalvcf_cmd = normalvcf_cmd_script1.format(gatk=gatk, index=index, normalbam=sample[1],
                                                         samplename=sample[0],
                                                         normalvcf=normalvcf)
        results.append(pool.apply_async(run_cmd, ("creat normal vcf file", "".join(normalvcf_cmd,))))
    pool.close()
    pool.join()
    [x.get() for x in results]
    normalargs = listofvcf(path_pon)
    ponvcf = os.path.join(path, "pon.vcf.gz")
    pon_cmd = pon_cmd_script.format(gatk=gatk, normalvcfs=normalargs, ponvcf=ponvcf)
    run_cmd("create panel of normals", "".join(pon_cmd))

def run_createsomatic_pon(path_pon,path):
    gatk = get_config("SNV", "GATK")
    normalargs = listofvcf(path_pon)
    ponvcf = os.path.join(path,"pon.vcf.gz")
    pon_cmd = pon_cmd_script.format(gatk=gatk, normalvcfs=normalargs, ponvcf=ponvcf)
    run_cmd("create panel of normals", "".join(pon_cmd))


