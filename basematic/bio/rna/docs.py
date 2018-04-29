doc = """
[Steps]
    #Single sample:
    basematic-RNA run_salmon -1 HMC_01.1.fq.gz -2 HMC_01.2.fq.gz -g hg38
    #Multiple samples:
    basematic-RNA run_salmon -m samples.txt -g hg38
    #Aggregate: TPM Counts and QC...
    basematic-RNA aggr_tpm_qc -m samples.txt -d ./salmon_process -o ./salmon
    
    #scatter plot...
    basematic-RNA plot_corelation_fig -1 sample03 -2 sample02 -t ./exptpm.txt
    #heatmap for a tpm... (correlation_heatmap.png)
    basematic-RNA corr_heatmap -t ./tpm.txt -n correlation_heatmap  
    
    basematic-RNA qc_tpm -o ./table_tpm.txt
    basematic-RNA diff_genes --groupfile ./file

[Files]
    TPM/Count file: table with column and row names;
    groupfile: A file with two columns: sample name/group label
    
[Step Detail]:
    qc_tpm ...
    
[to do]:
#check tpm files, list their samples...

"""

def print_doc():
    print(doc)

