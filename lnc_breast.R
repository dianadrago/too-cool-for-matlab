
setwd("/home/diana/Documents/yossi/tcga_breast/")

data <- read.table("unc.edu_BRCA_IlluminaHiSeq_RNASeqV2.geneExp.tsv", row.names = 1, header = T)

# filter protein coding genes and cases & normal adyacent tissue
library(NOISeq)
library(limma)
library(edgeR)
# phenoype
tumor <- grep('.01', colnames(data), value = TRUE)
hat <- grep('.11A|.11B', colnames(data), value = TRUE)
met <- grep('.06A', colnames(data), value = TRUE)
tumor <- cbind(tumor, "tumor")
hat <- cbind(hat, "HAT")
met <- cbind(met, "MET")
pheno <- as.data.frame(rbind(tumor,hat,met))

# data <- data[,pheno[,1]]
# group <- factor(pheno[,2])
# 
# y <- DGEList(counts=data,group=group)
# yTMM <- calcNormFactors(y, method="TMM")
# yUQ <- calcNormFactors(y, method="upperquartile")
# 
# TMM <- cpm(yTMM, normalized.lib.sizes=T)
# UQ <- cpm(yUQ, normalized.lib.sizes=T)
# 
# # density plot
# library(ggplot2)
# library(reshape2)
# 
# raw <- as.data.frame(data)
# TMM <- as.data.frame(TMM)
# 
# densityplot <- function(g){
#   M1 <- log10(g+1)
#   M2 <- cbind(rownames(M1), M1)
#   rownames(M2) <- NULL
#   colnames(M2)[1] <- "ID"
#   l2 = melt(M2)
#   l2
# }
# 
# plot_raw <- densityplot(raw)
# plot_TMM <- densityplot(TMM)
# 
# ggplot(plot_raw, aes(value)) + geom_density(aes(group = variable))
# ggplot(plot_TMM, aes(value)) + geom_density(aes(group = variable))
# 
# # RNAseq_TMM <- TMM[rowSums(x = TMM) == 0, , drop=FALSE]
# 
# ######################
# 
# RNAseq <- data[rowSums(x = data) >= length(colnames(data)), , drop=FALSE]
# 
# plot_RNAseq <- densityplot(RNAseq)
# ggplot(plot_RNAseq, aes(value)) + geom_density(aes(group = variable))
# 
# boxplot(log2(RNAseq))
# boxplot(log2(RNAseq_TMM))
################

non_pc <- read.table("non_cod_annot.txt", sep = '\t', header = T)
non_pc_id <- as.vector(non_pc[,2])
non_pc_id <- gsub('^ENTREZGENE_ACC:', '', non_pc_id)

breast_data <- cbind(rownames(data), data)
gene_entrez <- gsub('^.*[|]', '', rownames(data))
rownames(breast_data) <- gene_entrez
nc_breast <- breast_data[non_pc_id,]

library(DESeq2)

RNAseq <- round(RNAseq, 0)
