doc = """
[Salmon]
    #Single sample:
    baseq-RNA run_salmon -1 HMC_01.1.fq.gz -2 HMC_01.2.fq.gz -g hg38
    
    #Multiple samples run 4 samples at the same time ...
    baseq-RNA run_salmon -m samples.txt -g hg38 -n RNA
    
    #Aggregate: TPM Counts and QC...
    baseq-RNA aggr_tpm_qc -m samples.txt -d ./salmon -o ANZHEN_RNA4
    
[hisat+cufflinks]
    baseq-RNA run_hisat -m samples.txt -g hg38 -d ./hisat
    baseq-RNA run_hisat -1 HMC_01.1.fq.gz -2 HMC_01.2.fq.gz -g hg38

[star+featurecounts]
    baseq-RNA run_star -1 HMC_01.1.fq.gz -2 HMC_01.2.fq.gz -g hg38
    baseq-RNA run_featurecount -1 HMC_01.1.fq.gz -2 HMC_01.2.fq.gz -g hg38
    baseq-RNA run_cufflinks -1 HMC_01.1.fq.gz -2 HMC_01.2.fq.gz -g hg38
    
[visualization]    
    #scatter plot...
    baseq-RNA plot_corelation_scatter -1 sample03 -2 sample02 -t ./exptpm.txt
    
    #heatmap for a tpm... (correlation_heatmap.png)
    baseq-RNA corr_heatmap -t ./tpm.txt -n correlation_heatmap  
    
    baseq-RNA qc_tpm -o ./table_tpm.txt

[Differential Expression]
    baseq-RNA diff_power_analysis -g group.txt -p group_compare.txt -t ./tpm.txt -q qc.txt     
    baseq-RNA deseq2 -g group.txt -p group_compare.txt -t ./tpm.txt -c ./count.txt -o ./deseq_rep2
    
[Files]
    TPM/Count file: table with column and row names;
    groupfile: A file with two columns: sample name/group label
    
[Step Detail]:
    qc_tpm ...
    
[Quality Control]:
    baseq-RNA pca_analysis -t tpm.txt -g groups2.txt (generate PCA.png)
    
[to do]:
#check tpm files, list their samples...

[dropseq]
    #download barcode whitelist
    wget https://github.com/basedata10/DropRNA/raw/master/whitelist/whitelist.10X.txt
    wget https://github.com/basedata10/DropRNA/raw/master/whitelist/whitelist.indrop_1.txt	
    wget https://github.com/basedata10/DropRNA/raw/master/whitelist/whitelist.indrop_2.txt
        
    #All
    baseq-RNA drops_pipe -g hg38 -p 10X --cells 2000 --minreads 10000 -n 10X_1 -1 10X_1.1.fq.gz -2 10X_1.2.fq.gz -d ./ --parallel 4 --step tagging

    #barcode counting: Write the Barcode depth CSV to ./barcode_count.10X1.csv
    baseq-RNA drops_barcode_counting -p 10X -n 10X1 -1 10X_1.1.fq.gz -d ./
    baseq-RNA drops_barcode_counting -p indrop -n indrop_1 -1 inDrop_1.1.fq.gz -d ./
    baseq-RNA drops_barcode_counting -p dropseq -n dropseq -1 DropSeq_1.1.fq.gz -d ./    (=>barcode_count.dropseq.csv)
        
    #barcode stats: generate barcode_stats file which:
    #barcode/counts/mismatch_reads/mismatch_bc/mutate_last_base/total_reads
    baseq-RNA drops_barcode_stats -n 10X1 -p 10X --minreads 2000 -b ./barcode_count.10X1.csv (=>barcode_stats.10X1.csv)
    baseq-RNA drops_barcode_stats -n indrop1 -p indrop --minreads 1500 -b ./barcode_count.indrop_1.csv
    
    #split barcode:
    baseq-RNA drops_barcode_split -n 10X1 -p 10X -b ./barcode_stats.10X1.csv -1 10X_1.1.fq.gz -2 10X_1.2.fq.gz -d ./bcsplit_10X1

    #Star Alignment
    baseq-RNA drops_star_align_script -n drop1 -g hg38 -d ./bcsplit_dropseq1 -o ./star_aligns
    
    #reads Tagging: write generate file: barcode/genes UMI counts
    baseq-RNA drops_reads_tagging -b ./star_10X1/AA.sort.bam  -o ./tags.txt -g hg38
    
    #aggregate the barcodes...
    python 05_levels.py reads_indrop2 0.1 reads_indrop2/Counts_indrop2_1.txt reads_indrop2/Reads_indrop2_1.txt reads_indrop2/StatsUMI_indrop2_1.txt    
    baseq-RNA drops_barcode_gene_table -b ./10X_1/read_tagging/tagging.AA.txt
"""

def print_doc():
    print(doc)

