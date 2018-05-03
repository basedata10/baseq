
star_script = """
mkdir ${workdir}
cd ${workdir}
${bowtie2} -x ${bowtie2_index} -1 ${fq1} -2 ${fq2} --no-unal -p 20 -t -S bt2.sample.sam
${samtools} sort -n -@ 10 ${name}_Aligned.out.bam -o ${name}.sort.bam
rm bt2.sample.sam
"""

import os
from baseq.mgt.config import get_config
from string import Template

def genrate_bowtie_script(bowtie2, bowtie2_index, samtools, fq1, fq2, sample, workdir):

    if not os.path.exists(workdir):
        os.mkdir(workdir)

    def work_script(name, file):
        return Template(star_script).substitute(
            bowtie2 = bowtie2,
            bowtie2_index = bowtie2_index,
            samtools = samtools,
            workdir = workdir,
            name = name,
            file = os.path.abspath(file)
        )

    from itertools import product
    barcode_prefix = [x[0] + x[1] for x in list(product('ATCG', repeat=2))]
    bash_files = []

    for base in barcode_prefix:
        fastq = os.path.join("", "split.{}.{}.fq".format(sample, base))
        script = work_script(base, fastq)
        bash_path = os.path.join(workdir, "work_star_{}.sh".format(base))
        bash_files.append(bash_path)
        with open(bash_path, 'w') as file:
            file.writelines(script)

    from baseq.utils.workscript import write_bash, write_bash_qsub
    write_bash("./", bash_files, "work.align.sh")
    write_bash_qsub("./", bash_files, "work.qsub.align.sh")