setwd("/home/diana/Documents/yossi/Swati/expression_data/")

# adapted from https://cran.r-project.org/web/packages/MM2S/vignettes/MM2S_FromRawData.pdf
download.file(url = "http://mbni.org/customcdf/19.0.0/entrezg.download/hgu133plus2hsentrezgcdf_19.0.0.tar.gz",
  method = "auto",destfile = "hgu133plus2hsentrezgcdf_19.0.0.tar.gz")
install.packages("hgu133plus2hsentrezgcdf_19.0.0.tar.gz",type = "source",repos=NULL)

library(affy)
library(limma)
library(frma)
library(hgu133plus2hsentrezgcdf)

breast_affybatch <- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/breast_cancer_geo/GSE36774_RAW", cdfname="hgu133plus2hsentrezgcdf")
colon_affybatch <- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/colon_cancer_geo/GSE17536_RAW", cdfname="hgu133plus2hsentrezgcdf")
glioblastoma_affybatch <- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/glioblastoma_geo/GSE7696_RAW", cdfname="hgu133plus2hsentrezgcdf")
rabdo_affybatch <- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/rhabdomyosarcoma_ebi/files", cdfname="hgu133plus2hsentrezgcdf")
stomach_affybatch <- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/stomach_cancer_geo/GSE22377_RAW", cdfname="hgu133plus2hsentrezgcdf")
renal_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/renal_cancer_geo/GSE11151_RAW", cdfname="hgu133plus2hsentrezgcdf")
pancreatic_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/pancreatic_cancer_geo/GSE32676_RAW_and_GSE17891_RAW", cdfname="hgu133plus2hsentrezgcdf")
ovarian_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/ovarian_cancer_geo/GSE18520_RAW_and_GSE14001_RAW", cdfname="hgu133plus2hsentrezgcdf")
osteosarcoma_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/osteosarcoma_geo/GSE14827_RAW", cdfname="hgu133plus2hsentrezgcdf")
NSCLC_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/NSCLC_geo/GSE37745_RAW", cdfname="hgu133plus2hsentrezgcdf")
neuroblastoma_affybatch <- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/neuroblastoma_geo/GSE16476_RAW", cdfname="hgu133plus2hsentrezgcdf")
melanoma_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/melanoma_geo/melanoma", cdfname="hgu133plus2hsentrezgcdf")
liver_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/liver_cancer_geo/GSE9843_RAW", cdfname="hgu133plus2hsentrezgcdf")
ewing_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/ewing_sarcoma_geo/GSE34620_RAW", cdfname="hgu133plus2hsentrezgcdf")
aml_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/AML/GSE17855_RAW", cdfname="hgu133plus2hsentrezgcdf")
all_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/ALL/GSE28460_RAW", cdfname="hgu133plus2hsentrezgcdf")

normalize_custom <- function(affybatch, report_name){
  frmaData <- frma(affybatch, background="rma", normalize="quantile", summarize="robust_weighted_average")
  eset<-exprs(frmaData)
  pdf(file=report_name,width=7,height=5)
  N=length(affybatch@phenoData@data$sample)
  pm.mm=0
  for (i in 1:N) {pm.mm[i] = mean(mm(affybatch[,i])>pm(affybatch[,i]))}
  mycolors = rep(c("blue","red","green", "magenta"), each = 2)
  hist(affybatch, col=mycolors, main="Raw data distribution")
  boxplot(affybatch,col=mycolors, main="Raw data distribution")
  plot(100*pm.mm, type='h', main='Percent of MMs > PMs', ylab="%",xlab="Microarrays", ylim=c(0,50), col="red", lwd=5 )
  grid(nx = NULL, ny = 6, col = "blue", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)
  mycolors = rep(c("blue","red","green", "magenta"), each = 2)
  plotDensity(eset, col=mycolors, main="Data After normalization")
  boxplot(eset,col=mycolors, main="Data After normalization")
  dev.off()
  return(eset)
}

breast_report <- "breast_norm_report.pdf"
breast_test <- normalize_custom(breast_affybatch, breast_report)


