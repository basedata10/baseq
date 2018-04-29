doc = """
[Steps]
    #Single sample:
    baseq-RNA run_salmon -1 HMC_01.1.fq.gz -2 HMC_01.2.fq.gz -g hg38
    #Multiple samples:
    baseq-RNA run_salmon -m samples.txt -g hg38
    #Aggregate: TPM Counts and QC...
    baseq-RNA aggr_tpm_qc -m samples.txt -d ./salmon_process -o ./salmon
    
    #HISAT2 Pipeline [to do...]
    baseq-RNA run_hisat2 -1 HMC_01.1.fq.gz -2 HMC_01.2.fq.gz -g hg38
    baseq-RNA run_star -1 HMC_01.1.fq.gz -2 HMC_01.2.fq.gz -g hg38
    baseq-RNA run_featurecount -1 HMC_01.1.fq.gz -2 HMC_01.2.fq.gz -g hg38
    baseq-RNA run_cufflinks -1 HMC_01.1.fq.gz -2 HMC_01.2.fq.gz -g hg38
    
    #scatter plot...
    baseq-RNA plot_corelation_scatter -1 sample03 -2 sample02 -t ./exptpm.txt
    #heatmap for a tpm... (correlation_heatmap.png)
    baseq-RNA corr_heatmap -t ./tpm.txt -n correlation_heatmap  
    
    baseq-RNA qc_tpm -o ./table_tpm.txt
    baseq-RNA diff_genes --groupfile ./file

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

