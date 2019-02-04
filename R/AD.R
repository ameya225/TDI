require("here")
require("methods")
require("edgeR")
require("ggplot2")
require("ggpubr")

# Read data

protpath <- read.csv("Allen_Brain_Atlas/AD/ProteinAndPathologyQuantifications.csv")
donor <- read.csv("Allen_Brain_Atlas/AD/DonorInformation.csv")
stains <- read.csv("Allen_Brain_Atlas/AD/DescriptionOfStains.csv")
tbi <- read.csv("Allen_Brain_Atlas/AD/tbi_data_files.csv")
conversion_code <- read.csv("Allen_Brain_Atlas/AD/ConversionCode_AllenAgingDementiaTBI.csv")

# Gene expression data

genes <- read.csv("Allen_Brain_Atlas/AD/gene_expression_matrix_2016-03-03/rows-genes.csv")
samples <- read.csv("Allen_Brain_Atlas/AD/gene_expression_matrix_2016-03-03/columns-samples.csv")
class(samples$rnaseq_profile_id) <- "character"
samples_relev <- samples[,c("rnaseq_profile_id", "specimen_name")]
samples_relev <- dplyr::left_join(samples_relev, conversion_code, by="specimen_name")
samples_relev <- samples_relev[,c("rnaseq_profile_id", "specimen_name", "folderName")]

fpkm <- read.csv("AD/gene_expression_matrix_2016-03-03/fpkm_table_normalized.csv")
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
rownames(fpkm) <- fpkm$folderName

fpkm_numerical <- fpkm[,-c(50282,50283)]
fpkm_dup <- aggregate(fpkm_numerical,by=list(fpkm$folderName),data=fpkm_numerical,FUN=mean)
rownames(fpkm_dup) <- fpkm_dup$Group.1
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
fpkm <- dplyr::left_join(fpkm, donor, by="name")# This fpkm data frame contains all the metadata needed to build the model
sex <- fpkm$sex
diagnosis <- fpkm$act_demented
race <- fpkm$race
tbi <- fpkm$ever_tbi_w_loc


# PCA on the entire data
fpkm_gene <- fpkm[1:5419]
pca <- prcomp(fpkm_gene, center = T, scale. = T)
pca.score <- as.data.frame(pca$x)
pca.plot <- ggplot(data = pca.score, aes(x = PC1, y = PC2, colour= diagnosis)) + geom_point() + scale_shape_manual(1:24)
pca.plot


# Find differentially expressed genes by mean and variance
tmp <- CreateSeuratObject(counts = as.data.frame(t(fpkm_gene)))


meta <- read.csv("../akulkarn_1_30_2019_6_46_41.csv")




