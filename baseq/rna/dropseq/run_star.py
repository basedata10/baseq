
import os

star_script = """
cd ${workdir}
${STAR} --genomeLoad LoadAndRemove --genomeDir ${STAR_REF} \
    --readFilesIn $1 --runThreadN 10 --outSAMunmapped Within \
    --outSAMtype BAM Unsorted --outFileNamePrefix ${outPrefix}
#sort the bam file using name    
${samtools} sort -n -@ 10 ${name}.out.bam -o ${name}.sort.bam
"""
from baseq.mgt.config import get_config
from string import Template
def genrate_star_script(bc_dir, genome, workmode, workdir):
    STAR = get_config("RNA", "STAR")
    STAR_REF = get_config("RNA_ref_"+genome, "STAR_REF")
    samtools = get_config("RNA", "samtools")
    script = Template(star_script).substitute(
        STAR = STAR,
        STAR_REF = STAR_REF,
        samtools = samtools,
        workdir = workdir
    )

    #list fastq sampels
    files = os.listdir(bc_dir)
    for file in files:
        pass

    bases = [[x+y for y in "ATCG"] for x in "ATCG"]
    print(bases)
