from baseq.mgt import get_config, run_cmd

def bowtie2_sort(fq1, fq2, bamfile, genome, reads=5*1000*1000, thread=8):
    bowtie2 = get_config("CNV", "bowtie2")
    bowtie2_ref = get_config("CNV_ref_"+genome, "bowtie2_index")
    samtools = get_config("CNV", "samtools")

    samfile = bamfile+".sam"
    bamfile = bamfile
    statsfile = bamfile+".stat"

    print("[info] Bamfile Path : {}".format(bamfile))

    #Run Bowtie
    if fq1 and fq2:
        bowtie_cmd = [bowtie2, '-p', str(thread), '-x', bowtie2_ref, '-u', str(reads), '-1', fq1, '-2', fq2, '>', samfile]
    else:
        bowtie_cmd = [bowtie2, '-p', str(thread), '-x', bowtie2_ref, '-u', str(reads), '-U', fq1, '>', samfile]
    run_cmd("bowtie alignment", " ".join(bowtie_cmd))

    #run Samtools
    samtools_sort = [samtools, 'sort -@ ', str(thread), '-o', bamfile, samfile, ";", samtools, "index", bamfile, "; rm", samfile]
    run_cmd("samtools sort", " ".join(samtools_sort))

    #run flagstats
    cmd_stats = [samtools, "flagstat", bamfile, ">", statsfile]
    run_cmd("samtools stats", " ".join(cmd_stats))

    return bamfile