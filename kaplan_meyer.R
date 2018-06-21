# survival plots with R

# https://cran.r-project.org/web/packages/survminer/vignettes/Informative_Survival_Plots.html
setwd("/home/diana/Documents/yossi/Swati/expression_data/")

library(affy)
library(limma)
library(frma)
library(hgu133plus2hsentrezgcdf)
library(RColorBrewer)
library(plyr)
library(survminer)
library(survival)

ewing_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/ewing_sarcoma_geo/GSE34620_RAW", cdfname="hgu133plus2hsentrezgcdf")
ewing_report <- "ewing_norm_report.pdf"
ewing_exp <- normalize_custom(ewing_affybatch, ewing_report)
remove(ewing_affybatch)
ewing_exp <- normalizeBetweenArrays(ewing_exp)

names <- read.table("mart_export.txt", sep='\t', header = T)
names_sort <- na.omit(names[order(names$NCBI.gene.ID),])
names_sort <- unique(names_sort[,2:3])

km_patients <- colnames(ewing_exp)[grep("._R[0-9]", colnames(ewing_exp))]
ewing_exp_km <- ewing_exp[,km_patients]
colnames(ewing_exp_km) <- gsub('^.{10}|.{4}$', '', km_patients)

# define gene of interest
ewin_survival <- read.table("ewing_sarcoma_geo/GSE17618_Survival_data_Savola.txt", sep = '\t', row.names = 1, header = T,stringsAsFactors=FALSE)
ewin_survival <- ewin_survival[,1:7]
ewin_survival[ewin_survival=="Dead"] <-1
ewin_survival[ewin_survival=="NED" | ewin_survival=="AWD"]<-0


candidates <- c("A2M","ACTB","ANXA1","ANXA2","ANXA5","CALM2","COL1A2",
               "COL6A3","DAD1","DPYSL2","IARS","IGFBP7","MT2A","PPA1","PSMA6",
               "PSMB4","RPL3","RPL37A","S100A10","SELENOP","SRP9","UBE2E1","YWHAZ",
               "ZFAND5","EIF3D","EIF3H","COX7A2L","TMSB10","MRPL33","PRDX6","MAFB",
               "PDIA6","ATP6AP2","BASP1","PNRC1","TMED10","USP22","ISCU","SEC61G",
               "MMADHC","HSD17B12","NDUFA12","GOLPH3","TUBB6","MRFAP1","ZNF664",
               "IGLL1","PDE6C","PF4V1","SLC13A1","USH2A","MCF2L2","MORC1","PCDH15","ABCA13",
               "KIAA0825","EYS","IYD")

ewing_surv_plot <- function(gene){
gene_code <- as.character(names_sort[names_sort$Gene.name == gene,2])
gene_ewing <- t(ewing_exp_km[gene_code,])
gene_ewing <- t(colMeans(data.matrix(gene_ewing)))
rownames(gene_ewing) <- gene
dataset_clin <- ewin_survival[(as.vector(colnames(gene_ewing))),]
survival <- cbind(dataset_clin, exp = as.numeric(t(gene_ewing)))
colnames(survival)[8] <- gene
cut_ewing<- surv_cutpoint(survival, time = "OVS..months.", event="Status", variables = gene)
cut_ewing_plot <- surv_categorize(cut_ewing)
fit <- do.call(survival::survfit, list(Surv(OVS..months., Status == "1")~ cut_ewing_plot[,3], data = cut_ewing_plot))
ggsurvplot(fit = fit, risk.table = TRUE,  pval = TRUE,
           conf.int = T,
           risk.table.y.text.col = T, palette =c("darkorange", "mediumorchid4"),
           risk.table.y.text = FALSE,
           legend.labs=c(paste(gene,"= high"), paste(gene,"= low")))
}


lapply(candidates, ewing_surv_plot)







