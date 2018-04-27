#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library("ggplot2")
library("DESeq2")

tpm_file = args[1]
group_file = args[2]
out_file = args[3]

#Read Group Infos
groups = read.table(group_file)
samples = groups[,1]

df_sample_groups = do.call("rbind", apply(samples, 1,
                     function(x){data.frame(sample=x$samples, groups = x$name, stringsAsFactors = F)}))
rownames(df_sample_groups) = df_sample_groups$sample

#Read TPM File
tpm = read.table(tpm_file, header=T)
rownames(tpm) = tpm$gene
tpm = tpm[,2:dim(tpm)[2]]
storage.mode(tpm) = "integer"

print(tpm[1:10, samples])

#build DESeq2 Object

dds <- DESeqDataSetFromMatrix(
        countData = tpm[, samples],
        colData = df_sample_groups,
        design = ~ groups
    )
