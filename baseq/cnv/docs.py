
doc = """
[Whole Pipeline]
baseq-CNV run_pipeline ./Tn5_S1.fq.gz -g hg19

[Align]
baseq-CNV align -1 Tn5_S1.fq.gz -r 4000000 -g hg19 -t 10

[Bincount]
baseq-CNV bincount -g hg19 -i ./sample.bam -o normbincounts.txt

[Normalize]
baseq-CNV normalize -g hg19 -i ./bincounts.txt -o bincounts_norm.txt   

[CBS]
baseq-CNV cbs -i ./bincounts_norm.txt  -o ./out.file

[QC]: MAD...
baseq-CNV quality_control -i ./bincounts_norm.txt  -o ./out.file

[plot]
baseq-CNV plotgenome -i ./bincounts_norm.txt -c ./out.file

[File formats]
dynamicBin: CSV with header, 7 columns
    chr start abstart end range length GC
bincounts: Generate from "run_bincounting"
    tsv with index/counts/
normbincounts: Generated from "run_normalize"
    
"""

def print_doc():
    print(doc)