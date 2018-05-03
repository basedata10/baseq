import os

star_script = """
mkdir ${workdir}
cd ${workdir}
${STAR} --genomeLoad LoadAndRemove --genomeDir ${STAR_REF} \
--readFilesIn ${file} --runThreadN 10 --outSAMunmapped Within \
--outSAMtype BAM Unsorted --outFileNamePrefix ${name}_
  ${samtools} sort -n -@ 10 ${name}_Aligned.out.bam -o ${name}.sort.bam
rm ${name}_Aligned.out.bam
"""

from baseq.mgt.config import get_config
from string import Template
from baseq.mgt.command import run_cmd

def genrate_star_script(bc_dir, genome, sample, workdir, workmode):

    if not os.path.exists(workdir):
        os.mkdir(workdir)

    STAR = get_config("RNA", "star")
    STAR_REF = get_config("RNA_ref_" + genome, "star_index")
    samtools = get_config("RNA", "samtools")

    def work_script(name, file):
        return Template(star_script).substitute(
            STAR=STAR,
            STAR_REF=STAR_REF,
            samtools=samtools,
            workdir=workdir,
            name=name,
            file=os.path.abspath(file)
        )

    from itertools import product
    barcode_prefix = [x[0] + x[1] for x in list(product('ATCG', repeat=2))]
    bash_files = []

    for base in barcode_prefix:
        fastq = os.path.join(bc_dir, "split.{}.{}.fq".format(sample, base))
        script = work_script(base, fastq)
        bash_path = os.path.join(workdir, "work_star_{}.sh".format(base))
        bash_files.append(bash_path)
        with open(bash_path, 'w') as file:
            file.writelines(script)

    from baseq.utils.workscript import write_bash, write_bash_qsub
    write_bash("./", bash_files, "work.align.sh")
    write_bash_qsub("./", bash_files, "work.qsub.align.sh")

def run_star_multiple(bc_dir, workdir, sample, genome, parallel):
    from concurrent.futures import ThreadPoolExecutor
    pool = ThreadPoolExecutor(parallel)

    STAR = get_config("RNA", "star")
    STAR_REF = get_config("RNA_ref_" + genome, "star_index")
    samtools = get_config("RNA", "samtools")

    def work_script(name, file):
        return Template(star_script).substitute(
            STAR = STAR,
            STAR_REF = STAR_REF,
            samtools = samtools,
            workdir = os.path.join(workdir, name),
            name = name,
            file = os.path.abspath(file)
        )

    from itertools import product
    barcode_prefix = [x[0] + x[1] for x in list(product('ATCG', repeat=2))]
    for bc in barcode_prefix:
        fastq = os.path.join(bc_dir, "split.{}.{}.fq".format(sample, bc))
        script = work_script(bc, fastq)
        pool.submit(run_cmd, "Star Alignment", script)