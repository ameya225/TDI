require("here")
require("methods")
require("edgeR")
require("ggplot2")
require("ggpubr")
require("gridExtra")
require("grid")
require("Seurat")
require("clusterProfiler")


# Read data

protpath <- read.csv("AD/ProteinAndPathologyQuantifications.csv")
donor <- read.csv("AD/DonorInformation.csv")
stains <- read.csv("AD/DescriptionOfStains.csv")
tbi <- read.csv("AD/tbi_data_files.csv")
conversion_code <- read.csv("AD/ConversionCode_AllenAgingDementiaTBI.csv")

# Gene expression data

genes <- read.csv("AD/Gene_Expression_Matrix/rows-genes.csv")
samples <- read.csv("AD/Gene_Expression_Matrix/columns-samples.csv")
class(samples$rnaseq_profile_id) <- "character"
samples_relev <- samples[,c("rnaseq_profile_id", "specimen_name")]
samples_relev <- dplyr::left_join(samples_relev, conversion_code, by="specimen_name")
samples_relev <- samples_relev[,c("rnaseq_profile_id", "specimen_name", "folderName")]

fpkm <- read.csv("AD/Gene_Expression_Matrix/fpkm_table_normalized.csv")
colnames(fpkm) <- gsub(colnames(fpkm), pattern = "X", replacement = "")
colnames(fpkm)[1] <- "gene_id"
fpkm <- fpkm[complete.cases(fpkm),]

# Genes --> Gene symbols
fpkm <- dplyr::left_join(fpkm, genes, by="gene_id")
rownames(fpkm) <- fpkm$gene_symbol
fpkm <- fpkm[,-c(1,379:382)]

# Add Sample
fpkm <- as.data.frame(t(fpkm))
fpkm$rnaseq_profile_id <- rownames(fpkm)
fpkm <- dplyr::left_join(fpkm, samples_relev, by="rnaseq_profile_id")
# rownames(fpkm) <- fpkm$folderName

fpkm_numerical <- fpkm[,-c(50282,50283)]
fpkm_dup <- aggregate(fpkm_numerical,by=list(fpkm$folderName),data=fpkm_numerical,FUN=mean)
rownames(fpkm_dup) <- fpkm_dup$Group.1
foldernames <- rownames(fpkm_dup)
fpkm_dup <- fpkm_dup[,-c(1,50283)]
fpkm_dup <- fpkm_dup[complete.cases(fpkm_dup),] 
fpkm_t <- as.data.frame(t(fpkm_dup))
fpkm_t <- DGEList(fpkm_t)

# Remove lowly expressed genes

keep <- filterByExpr(fpkm_t)
fpkm_t <- fpkm_t[keep,] # 5419 genes; 377 samples
fpkm <- as.data.frame(t(fpkm_t$counts))

# Create factors for metadata and adjustment in modelling
name <- factor(unlist(lapply(strsplit(colnames(fpkm_t$counts), split = "[_]"), function(x){x[1]})))
name <- factor(unlist(lapply(strsplit(as.character(name), split = "[.]"), function(x){paste(x[1], x[2], x[3], sep = ".")})))
structure <- factor(unlist(lapply(strsplit(colnames(fpkm_t$counts), split = "[_]"), function(x){x[2]})))
fpkm <- data.frame(fpkm, name, structure)
fpkm <- dplyr::left_join(fpkm, donor, by="name") # This fpkm data frame contains all the metadata needed to build the model
rownames(fpkm) <- foldernames
sex <- fpkm$sex
diagnosis <- fpkm$act_demented
race <- fpkm$race
tbi <- fpkm$ever_tbi_w_loc


# PCA on the entire data
fpkm_gene <- fpkm[1:5419]
pca <- prcomp(fpkm_gene, center = T, scale. = T)
pca.score <- as.data.frame(pca$x)
pca.plot_1 <- ggplot(data = pca.score, aes(x = PC1, y = PC2, colour= structure)) + geom_point() + scale_shape_manual(1:24)
pca.plot_1 <- pca.plot_1 + ggtitle(ggtitle("PCA Plot Brain gene expression- All genes"))


genes.var <- apply(t(fpkm_gene), 1, var)
genes.var.select <- order(-genes.var)[1:200]
fpkm.t <- as.data.frame(t(fpkm))
fpkm.ts <- fpkm.t[genes.var.select,]
fpkm.ts <- as.matrix(fpkm.ts)
fpkm.ts <- apply(fpkm.ts, 1, as.numeric)

pca <- prcomp(fpkm.ts, center = T, scale. = T)
pca.score <- as.data.frame(pca$x)
pca.plot_2 <- ggplot(data = pca.score, aes(x = PC1, y = PC2, colour= structure)) + geom_point()
pca.plot_2 <- pca.plot_2 + ggtitle("PCA Plot Brain gene expression- Most variant genes")

pca <- prcomp(fpkm.ts, center = T, scale. = T)
pca.score <- as.data.frame(pca$x)
pca.plot_3 <- ggplot(data = pca.score, aes(x = PC1, y = PC2, colour= diagnosis)) + geom_point()
pca.plot_3 <- pca.plot_3 + ggtitle("PCA Plot Brain gene expression by diagnosis- Most variant genes")
pca.plot_3


# PCA plots

pca_plots <- grid.arrange(pca.plot_1, pca.plot_2, pca.plot_3, nrow=3)
ggsave("PCA_All_variant.tiff", pca_plots)
dev.off()


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




