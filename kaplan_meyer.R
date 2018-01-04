# survival plots with R

# https://cran.r-project.org/web/packages/survminer/vignettes/Informative_Survival_Plots.html

library(affy)
library(limma)
library(frma)
library(hgu133plus2hsentrezgcdf)
library(RColorBrewer)
library(plyr)

ewing_affybatch<- ReadAffy(celfile.path="/home/diana/Documents/yossi/Swati/expression_data/ewing_sarcoma_geo/GSE34620_RAW", cdfname="hgu133plus2hsentrezgcdf")
ewing_report <- "ewing_norm_report.pdf"
ewing_exp <- normalize_custom(ewing_affybatch, ewing_report)
remove(ewing_affybatch)
ewing_exp <- normalizeBetweenArrays(ewing_exp)

km_patients <- colnames(ewing_exp)[grep("._R[0-9]", colnames(ewing_exp))]
ewing_exp_km <- ewing_exp[,km_patients]
colnames(ewing_exp_km) <- gsub('^.{10}|.{4}$', '', km_patients)

# define gene of interest
gluco_ewing <- t(ewing_exp_km["2313",])

#hist(gluco_ewing)

library(survminer)
library(survival)
# library(RTCGA.clinical)
# survivalTCGA(BRCA.clinical, OV.clinical,
#             extract.cols = "admin.disease_code") -> BRCAOV.survInfo

# strate patients by median expression of gene of interest
ewin <- read.table("ewing_sarcoma_geo/GSE17618_Survival_data_Savola.txt", sep = '\t', row.names = 1, header = T,stringsAsFactors=FALSE)
ewin <- ewin[,1:7]
high <- gluco_ewing[,gluco_ewing[1,]>=median(gluco_ewing)]
high <- replace(high, 1:(length(high)), "High")
low <- gluco_ewing[,gluco_ewing[1,]<=median(gluco_ewing)]
low <- replace(low, 1:(length(low)), "Low")

Expression <- c(high, low)
ewin <- ewin[names(Expression),]
ewin <- cbind(ewin,Expression)
ewin[ewin=="Dead"]<-2
ewin[ewin=="NED" | ewin=="AWD"]<-1

#ewin_primary <- ewin[ewin$State == "Primary",]
#ewin <- ewin_primary

ewin$SurvObj <- with(ewin, Surv(OVS..months., Status == "2"))

km.as.one <- survfit(ewin$SurvOb ~ 1, data = ewin, conf.type = "log-log")
km.as.one
plot(km.as.one)

# library(rms)
km.by.exp <- survfit(SurvObj ~ Expression, data = ewin, conf.type = "log-log")
km.by.exp
ggsurvplot((km.by.exp), pval = T, conf.int = TRUE, title= "FLI1 expression OVS survival plot")

##################################

ewin$SurvObj_2 <- with(ewin, Surv(EFS..moths., Status == "2"))

km.as.one_EFS <- survfit(ewin$SurvObj_2 ~ 1, data = ewin, conf.type = "log-log")
km.as.one_EFS
plot(km.as.one_EFS)

km.by.exp_EFS <- survfit(SurvObj_2 ~ Expression, data = ewin, conf.type = "log-log")
km.by.exp_EFS
ggsurvplot((km.by.exp_EFS),pval = T, conf.int = TRUE, title= "FLI1 expression EFS survival plot")

###################################

# fit <- survfit(Surv(times, patient.vital_status) ~ admin.disease_code,
#                data = BRCAOV.survInfo)
# # Visualize with survminer
# ggsurvplot(fit, data = BRCAOV.survInfo, risk.table = TRUE)
# 
# plot(ewing_exp_km)


###################################

# heatmap
library(gplots)

gluco <- rbind(gluco_ewing,gluco_ewing)

my_palette <- rev(colorRampPalette(brewer.pal(11, "Spectral"))(n = 1000))

row.distance = dist(gluco, method = "euclidean")
row.cluster = hclust(row.distance, method = "ward.D2")

col.distance = dist(t(gluco), method = "euclidean")
col.cluster = hclust(col.distance, method = "ward.D2")


###############################################################################
## Plotting the Heatmap!! (where all colorful things happen...)
###############################################################################

heatmap.2(gluco,
          main = "" ,
          density.info = "none",
          key.xlab="log intensity",
          margins = c(4,6),
          trace = "none",
          col = my_palette,
          Rowv = as.dendrogram(row.cluster),
          keysize = 1.5,
          key.par=list(mar=c(3,0,3,0)),
          key.title = "",
          lhei=c(2,7))
