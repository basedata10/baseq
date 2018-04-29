library(Seurat)
library(dplyr)
library(Matrix)
library(matrixStats)

wang1.data <- Read10X(data.dir = "E:/exercise/cellranger/wang1/mm10/")
wang1 <- CreateSeuratObject(raw.data = wang1.data, min.cells = 3, min.genes = 250, 
                            project = "wang1_10X")
#add mt megadata
wang1.mito.genes <- grep(pattern = "^mm10_mt-", x = rownames(x = wang1@data), value = TRUE)
wang1.percent.mito <- Matrix::colSums(wang1@raw.data[wang1.mito.genes, ])/Matrix::colSums(wang1@raw.data)

wang1 <- AddMetaData(object = wang1, metadata = wang1.percent.mito, col.name = "percent.mito")
VlnPlot(object = wang1, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

#visualization of correlation
par(mfrow = c(1, 2))
GenePlot(object = wang1, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = wang1, gene1 = "nUMI", gene2 = "nGene")

wang1 <- FilterCells(object = wang1, subset.names = c("nGene", "percent.mito"),
                     low.thresholds = c(800, -Inf), high.thresholds = c(4000, 0.05))

wang1 <- NormalizeData(object = wang1, normalization.method = "LogNormalize", 
                       scale.factor = 10000)
par(mfrow = c(1, 1))
wang1 <- FindVariableGenes(object = wang1, mean.function = ExpMean, dispersion.function = LogVMR, 
                           x.low.cutoff = 0.005, x.high.cutoff = 3.5, y.cutoff = 0.5)

length(x = wang1@var.genes)

#PCA
wang1 <- ScaleData(object = wang1,vars.to.regress = c("nUMI", "percent.mito"))#, model.use = "negbinom"
wang1 <- RunPCA(object = wang1, pc.genes = wang1@var.genes, do.print = TRUE, pcs.print = 1:5, 
                genes.print = 5)
PCHeatmap(object = wang1, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE)
PCElbowPlot(object = wang1)

#TSNE
wang1 <- FindClusters(object = wang1, reduction.type = "pca", dims.use = 1:8, 
                      resolution = 0.6, print.output = 0, save.SNN = TRUE,force.recalc =TRUE)
wang1 <- RunTSNE(object = wang1, dims.use = 1:8, do.fast = TRUE)
TSNEPlot(object = wang1)

wang1.markers <- FindAllMarkers(object = wang1, only.pos = TRUE, min.pct = 0.25, 
                                thresh.use = 0.25)
top20_wang1 <- wang1.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
DoHeatmap(object = wang1, genes.use = top20_wang1$gene, slim.col.label = TRUE, remove.key = TRUE)