doc = """
[Functions]
basematic-fastq list_samples ./

#split barcode
basematic-fastq split_barcode ./barcode.txt R1.fq.gz ./ -s R1.fq &
basematic-fastq split_barcode ./barcode.txt R2.fq.gz ./ -s R2.fq &

#quality control
basematic-fastq qc -1 R1.fq.gz -2 R2.fq.gz
basematic-fastq qc -m samples.txt 
"""

def print_doc():
    print(doc)
