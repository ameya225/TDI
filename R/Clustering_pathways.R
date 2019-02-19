# Hierrarchical clustering on the filtered data
d <- dist(fpkm_gene)
hc <- hclust(d, method = "ward.D")
plot(hc, labels = diag_lab)

d.ts <- dist(fpkm.ts)
hc.ts <- hclust(d, method = "complete")
plot(hc.ts, labels = structure, cex.lab=0.5)


# Find differentially expressed genes by mean and variance
Gene_obj <- CreateSeuratObject(counts = as.data.frame(t(fpkm_gene)), min.cells = 180, min.features = 100)
Gene_obj$diagnosis <- diagnosis
Gene_obj$sex <- sex
Gene_obj$structure <- structure
Gene_obj$race <- race
Gene_obj$tbi <- tbi
VlnPlot(object = Gene_obj, features = c("nCount_RNA"), group.by = "diagnosis")
Gene_obj_norm <- NormalizeData(object = Gene_obj, normalization.method = "LogNormalize", scale.factor = 1e4)
Gene_obj_varfeat <- FindVariableFeatures(object = Gene_obj_norm, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, Inf), dispersion.cutoff = c(0.5, Inf))
length(x = VariableFeatures(object = Gene_obj_varfeat))
VariableFeaturePlot(Gene_obj_varfeat) + ggtitle("Mean-Dispersion plot: AD-genes")
Gene_obj_scale <- ScaleData(object = Gene_obj_varfeat, features = rownames(x = Gene_obj_varfeat), vars.to.regress = c("nCount_RNA", "structure", "sex", "race", "tbi"))
Gene_obj_pca <- RunPCA(object = Gene_obj_scale, features = VariableFeatures(object = Gene_obj_scale), verbose = FALSE)
DimPlot(Gene_obj_pca)
DimPlot(Gene_obj_pca, group.by = "diagnosis")
DimPlot(Gene_obj_pca, group.by = "sex")
DimPlot(Gene_obj_pca, group.by = "race")
DimPlot(Gene_obj_pca, group.by = "tbi")
DimPlot(Gene_obj_pca, group.by = "structure") + ggtitle("PCA Plot after regressing across tech variates")

x <- VariableFeatures(object = Gene_obj_varfeat)
eg = bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
kk <- enrichKEGG(gene         = eg,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
barplot(kk, drop=TRUE, showCategory=5)

# Differential expression analysis (log2 transformed FPKM- limma)

log_fpkm <- log2(t(fpkm_gene+1))
design <- model.matrix(~0+diagnosis+structure)
edat <- DGEList(log_fpkm)