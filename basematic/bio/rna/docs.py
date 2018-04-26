doc = """
[Steps]
    #Single sample:
    basematic-RNA run_salmon -1 HMC_01.1.fq.gz -2 HMC_01.2.fq.gz -g hg38
    #Multiple samples:
    basematic-RNA run_salmon -m samples.txt -g hg38
    
    basematic-RNA aggr_tpm --folder ./CNV_process -o ./table_tpm.txt
    basematic-RNA qc_tpm -o ./table_tpm.txt
    basematic-RNA diff_genes --groupfile ./file

[Files]
    TPM/Count file: table with column and row names;
    groupfile: A file with two columns: sample name/group label
    
[Step Detail]:
    qc_tpm ...
"""

def print_doc():
    print(doc)

