#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library("ggplot2")
library("DESeq2")

tpm_file = args[1]
group_file = args[2]
out_file = args[3]

#Read Group Infos
df_sample_groups = read.table(group_file,header=TRUE,row.names = 1)
samples = rownames(df_sample_groups)

#Read TPM File
tpm = read.table(tpm_file, header = TRUE,row.names = 1)
storage.mode(tpm) = "integer"

print(tpm[1:10, samples])

#build DESeq2 Object
dds <- DESeqDataSetFromMatrix(
        countData = tpm[, samples],
        colData = df_sample_groups,
        design = ~ groups)
#plot PCA
vsd_data <- vst(dds,blind = FALSE)
pcaData <- plotPCA(vsd_data,intgroup = "group",returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw()+ theme(panel.border = element_blank(),axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank())

#select variable genes
 dds <- dds[rowSums(counts(dds)) > 10,]
 dds <- DESeq(dds)

#using results() function to calculate padj, but how to transfer parameters ?
res1 <- results(dds,contrast = c("","","") #transfer para form colData