

bowtie2 = """
    ${bowtie2} -x ${ref} -1 ${fq1} -2 ${fq2} --no-unal -p ${thread} -t -S ${samfile}
    """

templete = """
${align}
${samtools_sort}
"""

samtools_sort = """
samtools view -b -u -S $1.sam >$1.bam
samtools sort -@ 6 $1.bam -o $1.sort.bam
samtools index $1.sort.bam
rm $1.bam $1.sam
"""

readme = 'Nothing found... 你来到了没有知识的荒原'