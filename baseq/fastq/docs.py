doc = """
[Functions]
baseq-fastq list_samples ./

[Sampling Fastq]


[split barcode]
baseq-fastq split_barcode ./barcode.txt R1.fq.gz ./ -s R1.fq &
baseq-fastq split_barcode ./barcode.txt R2.fq.gz ./ -s R2.fq &

[quality control]
baseq-fastq qc -1 R1.fq.gz -2 R2.fq.gz
baseq-fastq qc -m samples.txt 

#trim the fastq files...
baseq-fastq filter_polyat -m samples.txt -t 10 --seqfile ./seqs.txt
baseq-fastq filter_polyat -m samples.txt -t 5 --seqfile ./seqs.txt

##seqs.txt

"""

def print_doc():
    print(doc)