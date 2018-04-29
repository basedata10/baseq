
doc = """
[The Integrated Pipeline]
    baseq-CNV run_pipeline ./Tn5_S1.fq.gz -h hg19

[The Pipeline in Steps]
    baseq-CNV run_alignment ./Tn5_S1.fq.gz -g hg19
    baseq-CNV run_bincounting -g hg19 -i ./sample.bam -o bincounts.txt
    baseq-CNV run_normalize -g hg19 -i ./bincounts.txt -o normbincounts.txt
    baseq-CNV run_cbs -i ./sample.norm.txt -o ./out.file
    baseq-CNV run_plotgenome -i ./sample.norm.txt -c ./out2.file -o ./plot.ong

[Main file formats]
dynamicBin: CSV with header, 7 columns
    chr start abstart end range length GC
bincounts: Generate from "run_bincounting"
    tsv with index/counts/
normbincounts: Generated from "run_normalize"
    
[Steps]    
alignment:
    using bowtie2, for fast mode, the reads counts is limited to 5M reads;
bincounting...
"""

def print_doc():
    print(doc)