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

rabdo_affybatch <- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/rhabdomyosarcoma_ebi/files", cdfname="hgu133plus2hsentrezgcdf")
rabdo_report <- "rabdo_norm_report.pdf"
rabdo_exp <- normalize_custom(rabdo_affybatch, rabdo_report)
remove(rabdo_affybatch)
colnames(rabdo_exp) <- sprintf("rabdo_%s",seq(1:length(colnames(rabdo_exp))))

stomach_affybatch <- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/stomach_cancer_geo/GSE22377_RAW", cdfname="hgu133plus2hsentrezgcdf")
stomach_report <- "stomach_norm_report.pdf"
stomach_exp <- normalize_custom(stomach_affybatch, stomach_report)
remove(stomach_affybatch)
colnames(stomach_exp) <- sprintf("stomach_%s",seq(1:length(colnames(stomach_exp))))

renal_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/renal_cancer_geo/GSE11151_RAW", cdfname="hgu133plus2hsentrezgcdf")
renal_report <- "renal_norm_report.pdf"
renal_exp <- normalize_custom(renal_affybatch, renal_report)
remove(renal_affybatch)
colnames(renal_exp) <- sprintf("renal_%s",seq(1:length(colnames(renal_exp))))


pancreatic_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/pancreatic_cancer_geo/GSE32676_RAW_and_GSE17891_RAW", cdfname="hgu133plus2hsentrezgcdf")
pancreatic_report <- "pancreatic_norm_report.pdf"
pancreatic_exp <- normalize_custom(pancreatic_affybatch, pancreatic_report)
remove(pancreatic_affybatch)
colnames(pancreatic_exp) <- sprintf("pancreatic_%s",seq(1:length(colnames(pancreatic_exp))))


ovarian_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/ovarian_cancer_geo/GSE18520_RAW_and_GSE14001_RAW", cdfname="hgu133plus2hsentrezgcdf")
ovarian_report <- "ovarian_norm_report.pdf"
ovarian_exp <- normalize_custom(ovarian_affybatch, ovarian_report)
remove(ovarian_affybatch)

osteosarcoma_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/osteosarcoma_geo/GSE14827_RAW", cdfname="hgu133plus2hsentrezgcdf")
osteosarcoma_report <- "osteosarcoma_norm_report.pdf"
osteosarcoma_exp <- normalize_custom(osteosarcoma_affybatch, osteosarcoma_report)
remove(osteosarcoma_affybatch)

NSCLC_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/NSCLC_geo/GSE37745_RAW", cdfname="hgu133plus2hsentrezgcdf")
NSCLC_report <- "NSCLC_norm_report.pdf"
NSCLC_exp <- normalize_custom(NSCLC_affybatch, NSCLC_report)
remove(NSCLC_affybatch)

neuroblastoma_affybatch <- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/neuroblastoma_geo/GSE16476_RAW", cdfname="hgu133plus2hsentrezgcdf")
neuroblastoma_report <- "neuroblastoma_norm_report.pdf"
neuroblastoma_exp <- normalize_custom(neuroblastoma_affybatch, neuroblastoma_report)
remove(neuroblastoma_affybatch)

melanoma_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/melanoma_geo/melanoma", cdfname="hgu133plus2hsentrezgcdf")
melanoma_report <- "melanoma_norm_report.pdf"
melanoma_exp <- normalize_custom(melanoma_affybatch, melanoma_report)
remove(melanoma_affybatch)

liver_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/liver_cancer_geo/GSE9843_RAW", cdfname="hgu133plus2hsentrezgcdf")
liver_report <- "liver_norm_report.pdf"
liver_exp <- normalize_custom(liver_affybatch, liver_report)
remove(liver_affybatch)

ewing_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/ewing_sarcoma_geo/GSE34620_RAW", cdfname="hgu133plus2hsentrezgcdf")
ewing_report <- "ewing_norm_report.pdf"
ewing_exp <- normalize_custom(ewing_affybatch, ewing_report)
remove(ewing_affybatch)

aml_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/AML/GSE17855_RAW", cdfname="hgu133plus2hsentrezgcdf")
aml_report <- "aml_norm_report.pdf"
aml_exp <- normalize_custom(aml_affybatch, aml_report)
remove(aml_affybatch)

all_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/ALL/GSE28460_RAW", cdfname="hgu133plus2hsentrezgcdf")
all_report <- "all_norm_report.pdf"
all_exp <- normalize_custom(all_affybatch, all_report)
remove(all_affybatch)

breast_affybatch <- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/breast_cancer_geo/GSE36774_RAW", cdfname="hgu133plus2hsentrezgcdf")
breast_report <- "breast_norm_report.pdf"
breast_exp <- normalize_custom(breast_affybatch, breast_report)
remove(breast_affybatch)

colon_affybatch <- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/colon_cancer_geo/GSE17536_RAW", cdfname="hgu133plus2hsentrezgcdf")
colon_report <- "colon_norm_report.pdf"
colon_exp <- normalize_custom(colon_affybatch, colon_report)
remove(colon_affybatch)

glioblastoma_affybatch <- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/glioblastoma_geo/GSE7696_RAW", cdfname="hgu133plus2hsentrezgcdf")
glioblastoma_report <- "glioblastoma_norm_report.pdf"
glioblastoma_exp <- normalize_custom(glioblastoma_affybatch, glioblastoma_report)
remove(glioblastoma_affybatch)

cancers_names <- list(row.names(rabdo_exp), row.names(stomach_exp),row.names(renal_exp),row.names(pancreatic_exp),row.names(ovarian_exp),row.names(osteosarcoma_exp),row.names(NSCLC_exp),row.names(neuroblastoma_exp),row.names(melanoma_exp),row.names(liver_exp),row.names(ewing_exp),row.names(aml_exp),row.names(all_exp),row.names(breast_exp),row.names(colon_exp),row.names(glioblastoma_exp))

all(sapply(cancers_names, FUN = identical, row.names(rabdo_exp)))

cancers_todos <- cbind(rabdo_exp, stomach_exp,renal_exp,pancreatic_exp,ovarian_exp,osteosarcoma_exp,NSCLC_exp,neuroblastoma_exp,melanoma_exp,liver_exp,ewing_exp,aml_exp,all_exp,breast_exp,colon_exp,glioblastoma_exp)

cancer_norm <- normalizeBetweenArrays(cancers_todos)

rabdo_exp <- cancer_norm[,colnames(rabdo_exp)]
stomach_exp <- cancer_norm[,colnames(stomach_exp)]
pancreatic_exp <- cancer_norm[,colnames(pancreatic_exp)]
renal_exp <- cancer_norm[,colnames(renal_exp)]
ovarian_exp <- cancer_norm[,colnames(ovarian_exp)]
osteosarcoma_exp <- cancer_norm[,colnames(osteosarcoma_exp)]
NSCLC_exp <- cancer_norm[,colnames(NSCLC_exp)]
neuroblastoma_exp <- cancer_norm[,colnames(neuroblastoma_exp)]
melanoma_exp <- cancer_norm[,colnames(melanoma_exp)]
liver_exp <- cancer_norm[,colnames(liver_exp)]
glioblastoma_exp <- cancer_norm[,colnames(glioblastoma_exp)]
ewing_exp <- cancer_norm[,colnames(ewing_exp)]
breast_exp <- cancer_norm[,colnames(breast_exp)]
all_exp <- cancer_norm[,colnames(all_exp)]
aml_exp <- cancer_norm[,colnames(aml_exp)]
colon_exp <- cancer_norm[,colnames(colon_exp)]

rabdo <- t(rabdo_exp["190",])
stomach <- t(stomach_exp["190",])
pancreatic <- t(pancreatic_exp["190",])
renal <- t(renal_exp["190",])
ovarian <- t(ovarian_exp["190",])
osteosarcoma <- t(osteosarcoma_exp["190",])
NSCLC <- t(NSCLC_exp["190",])
neuroblastoma <- t(neuroblastoma_exp["190",])
melanoma <- t(melanoma_exp["190",])
liver <- t(liver_exp["190",])
glioblastoma <- t(glioblastoma_exp["190",])
ewing <- t(ewing_exp["190",])
breast <- t(breast_exp["190",])
all <- t(all_exp["190",])
aml <- t(aml_exp["190",])
colon <- t(colon_exp["190",])


cancers <- list(rabdo, stomach,renal,pancreatic,ovarian,osteosarcoma,NSCLC,neuroblastoma,melanoma,liver,ewing,aml,all,breast,colon,glioblastoma)
cancer_matrix <- t(rbind.fill.matrix(cancers))
colnames(cancer_matrix) <- c("Rhabdomyosarcoma (101)", "Stomach (43)","Renal (67)","Pancreatic (59)","Ovarian (83)","Osteosarcoma (27)","NSCLC (196)","Neuroblastoma (88)","Melanoma (89)","Liver (91)","Ewing Sarcoma (209)","AML (237)","ALL (98)","Breast (107)","Colon (177)","Glioblastoma (84)")

###############################################################################
# boxplot
###############################################################################
medians <- apply(cancer_matrix, 2, FUN = median, na.rm=TRUE)
order(medians)
cancer_matrix <- cancer_matrix[,order(medians)]
N = length(colnames(cancer_matrix))
mycolors = colorRampPalette(brewer.pal(11,"Spectral"))(N)
op <- par(mar = c(5, 15, 4, 2) + 0.1)
boxplot(cancer_matrix, horizontal = TRUE, las = 1, col=mycolors, main="NR0B1", log = "x", na.rm=TRUE)

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









