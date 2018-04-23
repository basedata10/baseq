library("ggplot2")
library("Cairo")
library("DESeq2")
library("jsonlite")
library("pheatmap")
library("clusterProfiler")

Args <- commandArgs()
countFile = Args[6]
tpmFile = Args[7]
configs = Args[8]
outPath = Args[9]

countFile="../../../datas/Bioinfo/RNA/Counts/Example/Example.Count.txt"
tpmFile="../../../datas/Bioinfo/RNA/TPM/Example/Example.TPM.txt"
configs="./config.json"
outPath="./output"

cfgs = read_json(configs, simplifyVector = T)
samples = cfgs$groups

#Sample Group Data.frames
df_sample_groups = do.call("rbind", apply(samples, 1,
                     function(x){data.frame(sample=x$samples, groups = x$name, stringsAsFactors = F)}))
rownames(df_sample_groups) = df_sample_groups$sample

df_counts <- as.matrix(read.csv(countFile, sep="\t", row.names="gene"))
df_tpm <- as.matrix(read.csv(tpmFile, sep="\t", row.names="gene"))
storage.mode(df_counts) = "integer"

#Build The DDS
dds <- DESeqDataSetFromMatrix(
        countData = df_counts[, rownames(df_sample_groups)],
        colData = df_sample_groups ,
        design = ~ groups)

keep <- rowSums(counts(dds)) >= 100
dds <- dds[keep,]
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)

##Transformation
vsd <- vst(dds, blind=FALSE)
df_PCA = plotPCA(vsd, intgroup=c("groups"), returnData = T)

#The PCA Plots...
file_MZ = paste(outPath, "/DEseq.PCA.svg", sep="")
Cairo(file_MZ, bg="white", type="svg", units="in",
        width=10, height=8, pointsize=13, dpi=120)
ggplot(df_PCA, aes(PC1, PC2, color=group))+
    geom_point(size=3)+
    theme_bw()
dev.off()

#Differntical Expressed Genes Testing...
diff_groups = cfgs$compare

volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot",
                         legendpos="bottomright", labelsig=FALSE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, cex = 0.5, main=main, ...))
  with(subset(res, padj<sigthresh ),
       points(log2FoldChange, -log10(pvalue), pch=20, cex = 0.5,  col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh),
       points(log2FoldChange, -log10(pvalue), pch=20, cex = 0.5, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh),
       points(log2FoldChange, -log10(pvalue), pch=20, cex = 0.5, col="green", ...))
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}

for(x in 1:dim(diff_groups)[1]){
  diff_name = diff_groups[x,]$name
  diff_group = diff_groups[x,]$groups[[1]]

  res <- results(dds, contrast = c("groups", diff_group[1], diff_group[2]))
  res <- lfcShrink(dds, coef=2, res = res)

  samples_1 = df_sample_groups$sample[df_sample_groups$groups == diff_group[1]]
  samples_2 = df_sample_groups$sample[df_sample_groups$groups == diff_group[2]]

  res$TPM_Group1 = rowMeans(df_tpm[rownames(res), samples_1])
  res$TPM_Group2 = rowMeans(df_tpm[rownames(res), samples_2])

  #MZ plot...
  file_MZ = paste(outPath, "/DEseq.MZ.", diff_name, ".png", sep="")
  CairoPNG(file_MZ, width=800, height=800, dpi = "auto", pointsize = 15)
  plotMA(res, ylim=c(-2,2))
  dev.off()

  #Valcano Plot...
  file_MZ = paste(outPath, "/DEseq.Volcano.", diff_name, ".png", sep="")
  CairoPNG(file_MZ, width=800, height=800, dpi = "auto", pointsize = 15)
  volcanoplot(res, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
  dev.off()

  genes_up = subset(res, log2FoldChange>=1 & padj<=0.05)[,c("log2FoldChange", "padj", "TPM_Group1", "TPM_Group2")]
  genes_down = subset(res, log2FoldChange<=-1 & padj<=0.05)[,c("log2FoldChange", "padj", "TPM_Group1", "TPM_Group2")]

  #HeatMap...
  file_MZ = paste(outPath, "/DEseq.HeatmapDiff.", diff_name, ".png", sep="")
  df_TPM_diffgenes = df_tpm[c(rownames(genes_up), rownames(genes_down)),
                            c(samples_1, samples_2)]
  CairoPNG(file_MZ, width=800, height=1000, dpi = "auto", pointsize = 15)
  pheatmap(df_TPM_diffgenes, cluster_cols = F, cluster_rows = T, scale = "row", labels_row = FALSE)
  dev.off()

  #GO Analysis....
  eg = bitr(rownames(genes_up), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  gene = eg$ENTREZID
  for(type in c("MF", "BP", "CC")){
    ego <- enrichGO(gene          = gene,
                    OrgDb         = org.Hs.eg.db,
                    ont           = type,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
    file_MZ = paste(outPath, "/DEseq.UP.", type, ".", diff_name, ".svg", sep="")
    CairoSVG(file_MZ, width=8, height=6, dpi = "auto", pointsize = 15)
    barplot(ego, showCategory=8)
    dev.off()

  }

  #Write The JSON FILE...
  file_JSON = paste(outPath, "/DEseq.genes.", diff_name, ".json", sep="")
  outs = toJSON(list(
    colnames = c("log2FoldChange", "padj", "TPM_Group1", "TPM_Group2"),
    down_regulate = as.data.frame(genes_down),
    up_regulate = as.data.frame(genes_up),
    groups = diff_name,
    samples_1 = samples_1,
    samples_2 = samples_2
    ), digits = 4, dataframe = "values")


  write(paste("window.", diff_name, "=", outs, sep=""), file_JSON)
}