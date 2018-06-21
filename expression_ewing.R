# expression analysis cancer types brainarrays annotation

setwd("/home/diana/Documents/yossi/Swati/expression_data/")

# adapted from https://cran.r-project.org/web/packages/MM2S/vignettes/MM2S_FromRawData.pdf
# download.file(url = "http://mbni.org/customcdf/19.0.0/entrezg.download/hgu133plus2hsentrezgcdf_19.0.0.tar.gz",
#   method = "auto",destfile = "hgu133plus2hsentrezgcdf_19.0.0.tar.gz")
# install.packages("hgu133plus2hsentrezgcdf_19.0.0.tar.gz",type = "source",repos=NULL)

library(affy)
library(limma)
library(frma)
library(hgu133plus2hsentrezgcdf)
library(RColorBrewer)
library(plyr)
library(gplots)
library(RColorBrewer)

# input: raw affy cel files
# output entrez id expression matrix and normalization report in pdf format

normalize_custom <- function(affybatch, report_name){
  frmaData <- frma(affybatch, background="rma", normalize="quantile", summarize="robust_weighted_average")
  eset<-exprs(frmaData)
  rownames(eset)<-gsub(rownames(eset),pattern="_at",replacement="")
  pdf(file=report_name,width=7,height=5)
  N=length(affybatch@phenoData@data$sample)
  pm.mm=0
  for (i in 1:N) {pm.mm[i] = mean(mm(affybatch[,i])>pm(affybatch[,i]))}
  mycolors = colorRampPalette(brewer.pal(11,"Spectral"))(N)
  hist(affybatch, col=mycolors, main="Raw data distribution")
  boxplot(affybatch,col=mycolors, main="Raw data distribution")
  plot(100*pm.mm, type='h', main='Percent of MMs > PMs', ylab="%",xlab="Microarrays", ylim=c(0,50), col="red", lwd=5 )
  grid(nx = NULL, ny = 6, col = "blue", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)
  mycolors = colorRampPalette(brewer.pal(11,"Spectral"))(N)
  plotDensity(eset, col=mycolors, main="Data After normalization")
  boxplot(eset,col=mycolors, main="Data After normalization")
  dev.off()
  return(eset)
}

#sprouty 10252
# gluco 2908

osteosarcoma_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/osteosarcoma_geo/GSE14827_RAW", cdfname="hgu133plus2hsentrezgcdf")
osteosarcoma_report <- "osteosarcoma_norm_report.pdf"
osteosarcoma_exp <- normalize_custom(osteosarcoma_affybatch, osteosarcoma_report)
remove(osteosarcoma_affybatch)

ewing_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/ewing_sarcoma_geo/GSE34620_RAW", cdfname="hgu133plus2hsentrezgcdf")
ewing_report <- "ewing_norm_report.pdf"
ewing_exp <- normalize_custom(ewing_affybatch, ewing_report)
remove(ewing_affybatch)

cancers_names <- list(row.names(osteosarcoma_exp),row.names(ewing_exp))

all(sapply(cancers_names, FUN = identical, row.names(ewing_exp)))

cancers_todos <- cbind(osteosarcoma_exp,ewing_exp)

cancer_norm <- normalizeBetweenArrays(cancers_todos)


osteosarcoma_exp <- cancer_norm[,colnames(osteosarcoma_exp)]
ewing_exp <- cancer_norm[,colnames(ewing_exp)]

names <- read.table("mart_export.txt", sep='\t', header = T)
names_sort <- na.omit(names[order(names$NCBI.gene.ID),])
names_sort <- unique(names_sort[,2:3])

# genes
chip_genes <- read.table("genes_chip_gr.txt", sep = '\t')
common_genes <- unique(names_sort[names_sort$Gene.name %in% intersect(as.vector(names_sort$Gene.name), as.vector(chip_genes$V1)),])
rownames(common_genes) <- common_genes$NCBI.gene.ID

osteosarcoma <- t(osteosarcoma_exp[intersect(rownames(common_genes), rownames(osteosarcoma_exp)),])
ewing <- t(ewing_exp[intersect(rownames(common_genes), rownames(ewing_exp)),])

gene_ID <- as.vector(common_genes[intersect(rownames(common_genes), rownames(ewing_exp)),1])

cancers <- list(osteosarcoma,ewing)
cancer_matrix <- t(rbind.fill.matrix(cancers))
colnames(cancer_matrix) <- c(rownames(osteosarcoma), rownames(ewing))
rownames(cancer_matrix) <- gene_ID

osteosarcoma_subset <- cancer_matrix[,1:27]
ewing_subset <- cancer_matrix[,28:236]

scaled.dat.osteo <- scale(osteosarcoma_subset)
genes_up_osteo <- scaled.dat.osteo[rowMeans(x = scaled.dat.osteo) > 2, , drop=FALSE]
genes_down_osteo <- scaled.dat.osteo[rowMeans(x = scaled.dat.osteo)< -1.5, , drop=FALSE]

scaled.dat.ewing <- scale(ewing_subset)
genes_up_ewing <- scaled.dat.ewing[rowMeans(x = scaled.dat.ewing) > 2, , drop=FALSE]
genes_down_ewing <- scaled.dat.ewing[rowMeans(x = scaled.dat.ewing)< -1.5, , drop=FALSE]

x1 <- rownames(genes_up_osteo)
x2 <- rownames(genes_down_osteo)
x3 <- rownames(genes_up_ewing)
x4 <- rownames(genes_down_ewing)

#venn diagram
venn(list("Genes Up Osteosarcoma (z-score > 1.5)"=x1, "Genes Up Ewing Sarcoma (z-score > 1.5)"=x3))
venn(list("Genes Down Osteosarcoma (z-score < -1.5)"=x2, "Genes Down Ewing Sarcoma (z-score < -1.5)"=x4))

intersect(x1,x3)
intersect(x2,x4)

genes_es <- unique(c(intersect(x1,x3), intersect(x2,x4)))

heatmap_genes <- cancer_matrix[genes_es,]
codes <- substr(colnames(heatmap_genes), 1, 9)
osteosarcoma_codes<- paste(codes[1:27], "- Osteosarcoma")
ewing_codes <- paste(codes[28:236], "- Ewing Sarcoma")
colnames(heatmap_genes) <- c(osteosarcoma_codes, ewing_codes)


my_palette <- rev(colorRampPalette(brewer.pal(11, "Spectral"))(n = 500))
my_palette2 <- rev(colorRampPalette(brewer.pal(9, "YlOrRd"))(n = 1000))
#my_palette <- rev(colorRampPalette(brewer.pal(11, "PuOr"))(n = 1000))

colLabels <-  c(rep("#377EB8", 27), rep("#E41A1C", 209))

row.distance = dist(heatmap_genes, method = "euclidean")
row.cluster = hclust(row.distance, method = "ward.D2")

col.distance = dist(t(heatmap_genes), method = "euclidean")
col.cluster = hclust(col.distance, method = "ward.D2")
###############################################################################
## Assign Column labels (Optional)
###############################################################################

###############################################################################
## Plotting the Heatmap!! (where all colorful things happen...)
###############################################################################

heatmap.2(heatmap_genes,
          main = "" ,
          density.info = "none",
          key.xlab="z-score",
          margins = c(5,25),
          trace = "none",
          col = my_palette,
          Rowv = as.dendrogram(row.cluster),
          keysize = 1,
          key.par=list(mar=c(3,1,3,1)),
          key.title = "",
          lhei=c(2,8),
          cexRow=1,
          sepwidth=c(0.07,0.07),
          Colv=as.dendrogram(col.cluster),
          sepcolor="white",
          scale = "column",
          labCol = FALSE,
          ColSideColors = colLabels)

par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       legend = c("Osteosarcoma", "Ewing Sarcoma"), # category labels
       col = c("dodgerblue", "firebrick1"),  # color key
       lty= 1,          # line style
       lwd = 5)   # line width

ewing_genes_only <- c(x3,x4)


heatmap_genes <- cancer_matrix[ewing_genes_only,]
codes <- substr(colnames(heatmap_genes), 1, 9)
osteosarcoma_codes<- paste(codes[1:27], "- Osteosarcoma")
ewing_codes <- paste(codes[28:236], "- Ewing Sarcoma")
colnames(heatmap_genes) <- c(osteosarcoma_codes, ewing_codes)


my_palette <- rev(colorRampPalette(brewer.pal(11, "Spectral"))(n = 500))
my_palette2 <- rev(colorRampPalette(brewer.pal(9, "YlOrRd"))(n = 1000))
#my_palette <- rev(colorRampPalette(brewer.pal(11, "PuOr"))(n = 1000))

colLabels <-  c(rep("#377EB8", 27), rep("#E41A1C", 209))

row.distance = dist(heatmap_genes, method = "euclidean")
row.cluster = hclust(row.distance, method = "ward.D2")

col.distance = dist(t(heatmap_genes), method = "euclidean")
col.cluster = hclust(col.distance, method = "ward.D2")
###############################################################################
## Assign Column labels (Optional)
###############################################################################

###############################################################################
## Plotting the Heatmap!! (where all colorful things happen...)
###############################################################################

heatmap.2(heatmap_genes,
          main = "" ,
          density.info = "none",
          key.xlab="z-score",
          margins = c(5,25),
          trace = "none",
          col = my_palette,
          Rowv = as.dendrogram(row.cluster),
          keysize = 1,
          key.par=list(mar=c(3,1,3,1)),
          key.title = "",
          lhei=c(2,8),
          cexRow=1,
          sepwidth=c(0.07,0.07),
          Colv=as.dendrogram(col.cluster),
          sepcolor="white",
          scale = "column",
          labCol = FALSE,
          ColSideColors = colLabels)

par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       legend = c("Osteosarcoma", "Ewing Sarcoma"), # category labels
       col = c("dodgerblue", "firebrick1"),  # color key
       lty= 1,          # line style
       lwd = 5)   # line width

# ###############################################################################
# # boxplot
# ###############################################################################
# medians <- apply(cancer_matrix, 2, FUN = median, na.rm=TRUE)
# order(medians)
# cancer_matrix <- cancer_matrix[,order(medians)]
# N = length(colnames(cancer_matrix))
# mycolors = colorRampPalette(brewer.pal(11,"Spectral"))(N)
# op <- par(mar = c(5, 15, 4, 2) + 0.1)
# boxplot(cancer_matrix, horizontal = TRUE, las = 1, col=mycolors, main="NR0B1", log = "x", na.rm=TRUE)

###############################################################################
library(beanplot)
# bean plot
###############################################################################

cancer <- as.data.frame(cancer_matrix)
N = length(colnames(cancer))
mycolors = colorRampPalette(brewer.pal(3,"Spectral"))(N)
op <- par(mar = c(5, 15, 4, 2) + 0.1)
beanplot(cancer, horizontal = TRUE, las = 1, na.rm=TRUE, col = mycolors)

###############################################################################
# heatmap
###############################################################################

my_palette <- rev(colorRampPalette(brewer.pal(11, "Spectral"))(n = 1000))

row.distance = dist(ewing_exp, method = "euclidean")
row.cluster = hclust(row.distance, method = "ward.D2")

col.distance = dist(t(ewing_exp), method = "euclidean")
col.cluster = hclust(col.distance, method = "ward.D2")


###############################################################################
## Plotting the Heatmap!! (where all colorful things happen...)
###############################################################################

heatmap.2(ewing_exp,
          main = "" ,
          density.info = "none",
          key.xlab="log intensity",
          margins = c(4,6),
          trace = "none",
          col = my_palette,
          Rowv = as.dendrogram(row.cluster),
          Colv = as.dendrogram(col.cluster), # apply selected clustering method
          keysize = 1.5,
          key.par=list(mar=c(3,0,3,0)),
          key.title = "",
          lhei=c(2,7),
          labCol = F)

##################################################################################
# 
# library(hgu133plus2frmavecs)
# 
# rabdo_affybatch <- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/rhabdomyosarcoma_ebi/files")
# rabdo_report <- "rabdo_norm_report.pdf"
# rabdo_exp <- normalize_custom(rabdo_affybatch, rabdo_report)
# remove(rabdo_affybatch)
# colnames(rabdo_exp) <- sprintf("rabdo_%s",seq(1:length(colnames(rabdo_exp))))
# 
# stomach_affybatch <- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/stomach_cancer_geo/GSE22377_RAW")
# stomach_report <- "stomach_norm_report.pdf"
# stomach_exp <- normalize_custom(stomach_affybatch, stomach_report)
# remove(stomach_affybatch)
# colnames(stomach_exp) <- sprintf("stomach_%s",seq(1:length(colnames(stomach_exp))))
# 
# renal_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/renal_cancer_geo/GSE11151_RAW")
# renal_report <- "renal_norm_report.pdf"
# renal_exp <- normalize_custom(renal_affybatch, renal_report)
# remove(renal_affybatch)
# colnames(renal_exp) <- sprintf("renal_%s",seq(1:length(colnames(renal_exp))))
# 
# 
# pancreatic_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/pancreatic_cancer_geo/GSE32676_RAW_and_GSE17891_RAW")
# pancreatic_report <- "pancreatic_norm_report.pdf"
# pancreatic_exp <- normalize_custom(pancreatic_affybatch, pancreatic_report)
# remove(pancreatic_affybatch)
# colnames(pancreatic_exp) <- sprintf("pancreatic_%s",seq(1:length(colnames(pancreatic_exp))))
# 
# 
# ovarian_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/ovarian_cancer_geo/GSE18520_RAW_and_GSE14001_RAW")
# ovarian_report <- "ovarian_norm_report.pdf"
# ovarian_exp <- normalize_custom(ovarian_affybatch, ovarian_report)
# remove(ovarian_affybatch)
# 
# osteosarcoma_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/osteosarcoma_geo/GSE14827_RAW")
# osteosarcoma_report <- "osteosarcoma_norm_report.pdf"
# osteosarcoma_exp <- normalize_custom(osteosarcoma_affybatch, osteosarcoma_report)
# remove(osteosarcoma_affybatch)
# 
# NSCLC_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/NSCLC_geo/GSE37745_RAW")
# NSCLC_report <- "NSCLC_norm_report.pdf"
# NSCLC_exp <- normalize_custom(NSCLC_affybatch, NSCLC_report)
# remove(NSCLC_affybatch)
# 
# neuroblastoma_affybatch <- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/neuroblastoma_geo/GSE16476_RAW")
# neuroblastoma_report <- "neuroblastoma_norm_report.pdf"
# neuroblastoma_exp <- normalize_custom(neuroblastoma_affybatch, neuroblastoma_report)
# remove(neuroblastoma_affybatch)
# 
# melanoma_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/melanoma_geo/melanoma")
# melanoma_report <- "melanoma_norm_report.pdf"
# melanoma_exp <- normalize_custom(melanoma_affybatch, melanoma_report)
# remove(melanoma_affybatch)
# 
# liver_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/liver_cancer_geo/GSE9843_RAW")
# liver_report <- "liver_norm_report.pdf"
# liver_exp <- normalize_custom(liver_affybatch, liver_report)
# remove(liver_affybatch)
# 
# ewing_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/ewing_sarcoma_geo/GSE34620_RAW")
# ewing_report <- "ewing_norm_report.pdf"
# ewing_exp <- normalize_custom(ewing_affybatch, ewing_report)
# remove(ewing_affybatch)
# 
# aml_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/AML/GSE17855_RAW")
# aml_report <- "aml_norm_report.pdf"
# aml_exp <- normalize_custom(aml_affybatch, aml_report)
# remove(aml_affybatch)
# 
# all_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/ALL/GSE28460_RAW")
# all_report <- "all_norm_report.pdf"
# all_exp <- normalize_custom(all_affybatch, all_report)
# remove(all_affybatch)
# 
# breast_affybatch <- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/breast_cancer_geo/GSE36774_RAW")
# breast_report <- "breast_norm_report.pdf"
# breast_exp <- normalize_custom(breast_affybatch, breast_report)
# remove(breast_affybatch)
# 
# colon_affybatch <- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/colon_cancer_geo/GSE17536_RAW")
# colon_report <- "colon_norm_report.pdf"
# colon_exp <- normalize_custom(colon_affybatch, colon_report)
# remove(colon_affybatch)
# 
# glioblastoma_affybatch <- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/glioblastoma_geo/GSE7696_RAW")
# glioblastoma_report <- "glioblastoma_norm_report.pdf"
# glioblastoma_exp <- normalize_custom(glioblastoma_affybatch, glioblastoma_report)
# remove(glioblastoma_affybatch)

###############################################################################
###############################################################################









