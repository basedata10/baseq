#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library("ggplot2")
library("DESeq2")
library("Cairo")

tpm_file = args[1]
count_file = args[2]
group_file = args[3]
group_compare_file = args[4]
out_path = args[5]

#Read Group Infos
df_sample_groups = read.table(group_file,header=TRUE,row.names = 1)
samples = make.names(rownames(df_sample_groups))
rownames(df_sample_groups) = samples

#Read Group Comparation File...
df_comp = read.table(group_compare_file,header=TRUE, stringsAsFactors =FALSE)
print(df_comp)

#Read TPM File
print(paste("Read TPM File ", tpm_file))
df_tpm = read.table(tpm_file, header = TRUE, row.names = 1)

#Read count_file
print(paste("Read Count File ", count_file))
df_counts = read.table(count_file, header = TRUE, row.names = 1)

#Check sample names in TPM table
#to be added...

#build DESeq2 Object
dds <- DESeqDataSetFromMatrix(
        countData = round(df_counts[, samples], 0),
        colData = df_sample_groups,
        design = ~ groups
    )

#Plot PCA
vsd_data <- vst(dds,blind = FALSE)
pcaData <- plotPCA(vsd_data, intgroup = "groups", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

figurefile = "./PCA.png"
CairoPNG(figurefile, width=700, height=600)

p2 = ggplot(pcaData, aes(PC1, PC2, color=group)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    theme_bw()+
    theme(panel.border = element_blank(),axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank())
print(p2)
dev.off()

compare_two_groups <- function(group1, group2){

    samples_group1 = rownames(subset(df_sample_groups, groups==group1))
    samples_group2 = rownames(subset(df_sample_groups, groups==group2))

    dds <- dds[rowSums(counts(dds)) > 10,]
    dds <- DESeq(dds)

    res_compare <- results(dds, contrast = c("groups", group1, group2))

    TPM_name_G1 = make.names(paste("TPM", group1, sep="_"))
    TPM_name_G2 = make.names(paste("TPM", group2, sep="_"))

    res_compare[TPM_name_G1] = rowMeans(df_tpm[rownames(res_compare), samples_group1])
    res_compare[TPM_name_G2] = rowMeans(df_tpm[rownames(res_compare), samples_group2])

    genes_down = subset(res_compare, log2FoldChange>=1 & padj<=0.05)[,c("log2FoldChange", "padj", TPM_name_G1, TPM_name_G2)]
    genes_up = subset(res_compare, log2FoldChange<=-1 & padj<=0.05)[,c("log2FoldChange", "padj", TPM_name_G1, TPM_name_G2)]

    return(list(genes_up, genes_down))

}


for (x in 1:dim(df_comp)[1]){
    diff_genes = compare_two_groups(df_comp[x, 1], df_comp[x, 2])
    genes_up_path = paste(out_path, "/genes_up.", df_comp[x, 1], "_", df_comp[x, 2], ".txt", sep="")
    genes_down_path = paste(out_path, "/genes_down.", df_comp[x, 1], "_", df_comp[x, 2], ".txt", sep="")
    print(genes_up_path)
    write.table(diff_genes[[1]], genes_up_path)
    write.table(diff_genes[[2]], genes_down_path)
}

diff_genes = compare_two_groups("C1T1", "C1T2")
diff_genes = compare_two_groups("C1T1", "C1T3")
diff_genes = compare_two_groups("C1T2", "C1T3")
