templete_DESeq = """
library("DESeq2")
dds <- DESeqDataSetFromMatrix(
          countData = geneMatrix,
          colData = groupDatas,
          design = ~ groups)
dds <- dds[keep,]
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
res <- results(dds, contrast = c("groups", diff_group[1], diff_group[2]))
res <- lfcShrink(dds, coef=2, res = res)
genes_up = subset(res, log2FoldChange>=${FC} & ${pvalue}<=0.05)[,c("log2FoldChange", "padj", "TPM_Group1", "TPM_Group2")]
genes_down = subset(res, log2FoldChange<=-${FC} & ${pvalue}<=0.05)[,c("log2FoldChange", "padj", "TPM_Group1", "TPM_Group2")]
"""

templete_SampleFiltering = """
    keep <- rowSums(counts(dds)) >= 100
"""

tps_Valcano = """
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot",
                         legendpos="bottomright", labelsig=FALSE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, cex = 0.5, main=main, ...))
  with(subset(res, padj<sigthresh ),
       points(log2FoldChange, -log10(pvalue), pch=20, cex = 0.5,  col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh),
       points(log2FoldChange, -log10(pvalue), pch=20, cex = 0.5, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh),
       points(log2FoldChange, -log10(pvalue), pch=20, cex = 0.5, col="green", ...))
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""),
    paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
"""

tps_heatmap = """
df_TPM_diffgenes = df_tpm[c(rownames(genes_up), rownames(genes_down)),
                        c(samples_1, samples_2)]
CairoPNG(file_MZ, width=800, height=1000, dpi = "auto", pointsize = 15)
pheatmap(df_TPM_diffgenes, cluster_cols = F, cluster_rows = T, scale = "row", labels_row = FALSE)
dev.off()
"""

class DESeq2:

    def __init__(self):
        self.type = "deseq2"
        self.DE_FC = 1
        self.DE_pvalue = 0.05
        self.groups = []

    def add_groups(self):
        pass

    def addVisualize(self):
        pass

    def set_FC_PValue(self, FC, pvalue):
        self.DE_FC = FC
        self.DE_pvalue = pvalue