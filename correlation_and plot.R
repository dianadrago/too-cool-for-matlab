###############################################################################
# Correlation co expressed genes
###############################################################################
setwd("/home/diana/Documents/yossi/Swati/expression_data/")

library(affy)
library(limma)
library(frma)
library(hgu133plus2hsentrezgcdf)
library(RColorBrewer)
library(plyr)
library(stats)

###############################################################################
# Ewing sarcoma
###############################################################################


ewing_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/ewing_sarcoma_geo/GSE34620_RAW", cdfname="hgu133plus2hsentrezgcdf")
ewing_report <- "ewing_norm_report.pdf"
ewing_exp <- normalize_custom(ewing_affybatch, ewing_report)
remove(ewing_affybatch)

ewing_data <- normalizeBetweenArrays(t(ewing_exp))

gluco <- as.vector(ewing_data[,"2908"])
#0.4024772 fli1 ews
fli1 <- as.vector(ewing_data[,"2313"])
ews <- as.vector(ewing_data[,"2130"])
# foxo1 fli1 0.3258583
foxo1 <- as.vector(ewing_data[,"2308"])

pearson_val <- round(cor(ews, fli1, method = "pearson"),digits = 2)
kendall_val <- round(cor(ews, fli1, method = "kendall"),digits = 2)
test_p <- cor.test(ews, fli1, method = "pearson")
test_s <- cor.test(fli1, ews, method = "spearman")
mycolors = colorRampPalette(brewer.pal(11,"Spectral"))(length(ews))
par(bg="black", fg="white", col.lab="white", col.main="white", font.lab=2)
plot(x=fli1, y=ews, main="EWS - FLI1 expression", log = "xy", xlab="NR3C1", ylab="FLI1", col=mycolors, pch=16, bty="n")
legend("bottomright", inset=.05, title="Ewing Sarcoma (209)",bty = "n",
       c(paste("Pearson:", pearson_val, " p-val:", signif(test_p$p.value, digits = 3), sep = " "),paste("Spearman:", sparman_val, " p-val:", signif(test_s$p.value, digits = 3),sep = " ")))
abline(lm(ews ~ fli1), untf=TRUE)

#########################
NSCLC
#########################


NSCLC_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/NSCLC_geo/GSE37745_RAW", cdfname="hgu133plus2hsentrezgcdf")
NSCLC_report <- "NSCLC_norm_report.pdf"
NSCLC_exp <- normalize_custom(NSCLC_affybatch, NSCLC_report)
remove(NSCLC_affybatch)

lung_data <- normalizeBetweenArrays(t(NSCLC_exp))

pdl1 <- as.vector(lung_data[,"29126"])
rela <- as.vector(lung_data[,"5970"])
kb1 <- as.vector(lung_data[,"4790"])

plot(pdl1, kb1, main="PDL1 vs KB1",  log = "xy")
abline(lm(kb1 ~ pdl1), untf=TRUE)
cor(kb1, pdl1, method = "spearman")

###############################################################################
# non lineal correlation measures- empirical mutual information
###############################################################################

library(minet)


lung_mim <- build.mim(dataset, estimator = "mi.mm", disc = "none", nbins = sqrt(NROW(dataset)))