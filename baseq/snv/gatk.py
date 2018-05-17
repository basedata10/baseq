import os
from baseq.mgt import get_config, run_cmd

# PICARD = get_config("SNV", "picard")
# GATK = get_config("SNV", "GATK")

bwa_cmd_script_p = r"""{bwa} mem -t {thread} -M -R "@RG\tID:{sample}\tSM:{sample}\tLB:WES\tPL:Illumina" {genome} {fq1} {fq2}  1>{samfile}"""
bwa_cmd_script_s = r"""{bwa} mem -t {thread} -M -R "@RG\tID:{sample}\tSM:{sample}\tLB:WES\tPL:Illumina" {genome} {fq1} 1>{samfile}"""
sort_index_cmd_script = """
{samtools} view -b -u -S {samfile}>{viewedbam}
{samtools} sort -@ 8 {viewedbam} {sample}
{samtools} index {sample}.bam
rm {samfile} {viewedbam}
"""

def alignment(fq1, fq2, sample, genome, thread=8):
    """
    Map low-divergent sequences against reference genome using BWA.
    Add ReadGroup(more details about ReadGroup_ )to bamfile using the input sample name.
    Outfile is in BAM format and indexed for the downstream analysis.

    .. _ReadGroup: https://software.broadinstitute.org/gatk/documentation/article.php?id=6472

    Usage:
    ::
      baseq-SNV run_bwa -1 Reads.1.fq.gz -2 Read.2.fq.gz -g hg38 -n Test

    Return:
    ::
      Test.bam
      Test.bam.bai

    """
    bwa = get_config("SNV", "bwa")
    samtools = get_config("SNV", "samtools")
    genome = get_config("SNV_ref_"+genome, "bwa_index")
    viewedbam = sample + ".view.bam"
    samfile = sample + ".sam"
    if fq1 and fq2 and os.path.exists(fq1) and os.path.exists(fq2):
        bwa_cmd = bwa_cmd_script_p.format(bwa=bwa, sample=sample, genome=genome, fq1=fq1, fq2=fq2, samfile=samfile, thread=thread)
    elif fq1 and os.path.exists(fq1):
        bwa_cmd = bwa_cmd_script_s.format(bwa=bwa, sample=sample, genome=genome, fq1=fq1, samfile=samfile, thread=thread)
    sort_index_cmd=sort_index_cmd_script.format(samtools=samtools, sample=sample, samfile=samfile,viewedbam=viewedbam)
    run_cmd("bwa alignment", "".join(bwa_cmd))
    run_cmd("samtools sort", "".join(sort_index_cmd))
    return bwa_cmd+"\n"+sort_index_cmd

markdup_cmd_script ="""
{java} -jar {picard} MarkDuplicates INPUT={bamfile} OUTPUT={markedbam} METRICS_FILE={markedbam}.metrics
{samtools} index {markedbam}
"""

def run_markdup(bamfile, markedbam):
    """
    Run MarkDuplicates of Picard. this function tags duplicate reads with "markduplicate" in BAM file.
    See also MarkDuplicates_ in GATK.


    .. _MarkDuplicates: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_sam_markduplicates_MarkDuplicates.php

    Usage:
    ::
      baseq-SNV run_markdup -b Test.bam -m Test.marked.bam

    Return:
    metrics file indicates the numbers of duplicates for both single- and paired-end reads
    ::
      Test.marked.bam
      Test.marked.bam.bai
      Test.marked.bam.metrics
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
    Run BQSR_. This function performs the two-steps process called base quality score recalibration. the first
    ster generates a recalibration table based on various covariates which is recruited to the second step to
    correct the systematic bias in input BAM file. More details about BaseRecalibrator_ and ApplyBQSR_ .


    .. _BQSR: https://gatkforums.broadinstitute.org/gatk/discussion/44/base-quality-score-recalibration-bqsr
    .. _BaseRecalibrator: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_BaseRecalibrator.php
    .. _ApplyBQSR: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_ApplyBQSR.php

    Usage:

    * Default mode filters duplicate reads (reads with "markduplicate" tags) before applying BQSR
      ::
         baseq-SNV run_bqsr -m Test.marked.bam -g hg38 -q Test.marked.bqsr.bam

    * Disable reads filter before analysis.
      ::
        baseq-SNV run_bqsr -m Test.marked.bam -g hg38 -q Test.marked.bqsr.bam -f Yes

    Return:
    ::
      Test.marked.bam.table
      Test.marked.bqsr.bai
      Test.marked.bqsr.bam
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
    """
    Call germline SNPs and indels via local re-assembly of haplotypes. BAM file recalbrated by BQSR do recommand as
    input BAM file and this functin only run the single sample genotypeVCF calling. More details see also
    HaplotypeCaller_

    .. _HaplotypeCaller: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php

    Usage:
    ::
      baseq-SNV run_callvar -q Test.marked.bqsr.bam -r Test.raw.indel.snp.vcf -g hg38

    Return:
    ::
      Test.raw.indel.snp.vcf

    """
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
    This function selects SNPs from a VCF file which is usually the output file of
    HaplotypeCaller. Then, all SNPs are filtered by certain criteria based on INFO and/or FORMAT annotations.
    Criteria used here is "QD < 2.0 || FS > 60.0 || MQ < 40.0".
    More details about SelectVariants_ and VariantFiltration_

    .. _SelectVariants: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_variantutils_SelectVariants.php
    .. _VariantFiltration: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_filters_VariantFiltration.php

    Usage:
    ::
      baseq-SNV run_selectvar -r Test.raw.indel.snp.vcf -s Test.raw.snp.vcf -f Test.filtered.snp.vcf -g hg38

    Return:
    ::
      Test.raw.snp.vcf
      Test.filtered.snp.vcf
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
    """
    Mutect2 is aim to call somatic SNVs and indels via local assembly of haplotypes. This function requires both
    tumor BAM file and its matched normal BAM file. tumorname and normalname should be consistent with the ReadGroup(ID) of tumor
    BAM file and normal BAM file respectively. PoN is refer to panel of normals callset(more infomation about PoN and how to
    create it, please see PoN_ ). Germline resource, also in VCF format, is used to annotate variant alleles. Default germline resource is
    downloaded from here_ .

    .. _here: https://software.broadinstitute.org/gatk/download/bundle
    .. _PoN: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_mutect_CreateSomaticPanelOfNormals.php

    Usage:

    * Simplified Mutect2 command line
      ::
        # single sample
        baseq-SNV run_mutect2 -g hg37 -n normal -N normal_marked_bqsr.bam \\
                                      -t tumor -T tumor_marked_bqsr.bam -o ./
        # multiple samples
        baseq-SNV run_mutect2 -g hg37 -l sample_list.txt -o ./

    * Specify PoN(panels of normals) VCF file and germline VCF file
      Default germline VCF file comes form GATK resource bundle and is recruited if germline isn't designated.
      ::
        # single sample
        baseq-SNV run_mutect2 -g hg37 -n normal -N normal_marked_bqsr.bam \\
                                      -t tumor -T tumor_marked_bqsr.bam -o ./ \\
                                      -p pon.vcf.gz -G af-only-gnomad.raw.sites.b37.vcf.gz
        # multiple samples
        baseq-SNV run_mutect2 -g hg37 -l sample_list.txt -o ./ \\
                                      -p pon.vcf.gz -G af-only-gnomad.raw.sites.b37.vcf.gz

    """



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
    """
    Create_pon function helps you create PoN(panel of normals) file necessary for mutect2 calling. The PoN captures
    common artifactual and germline variant sites. Mutect2 then uses the PoN to filter variants at the site-level.

    Example of samples list (tab delimited):

    * Content of columns are normal sample name, normal BAM file, tumor sample name, tumor BAM file(order cannot be distruped).
      Absoulte path of all BAM files should be added if directory of BAM files and analysis directory are different.
      ::
        N504    N504_marked_bqsr.bam   T504    T504_marked_bqsr.bam
        N505    N505_marked_bqsr.bam   T505    T505_marked_bqsr.bam
        N506    N506_marked_bqsr.bam   T506    T506_marked_bqsr.bam
        N509    N509_marked_bqsr.bam   T509    T509_marked_bqsr.bam
        N510    N510_marked_bqsr.bam   T510    T510_marked_bqsr.bam

    Usage:
    
    * Interval list defines genomic regions where analysis is restricted. Introduction of interval list format and its function, please see here_.
      ::
        # designated a intervals.list
        baseq-SNV create_pon -g hg37 -l sample_list.txt -p ./ -L interval.list
        # Using the dafalut intervals.list
        baseq-SNV create_pon -g hg37 -l sample_list.txt -p ./

    .. _here: https://software.broadinstitute.org/gatk/documentation/article?id=11009
    """
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


