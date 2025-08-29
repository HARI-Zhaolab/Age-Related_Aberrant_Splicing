setwd('/home/cmx/senes_splicing/')
rm(list=ls())
load("./outdata/step1_cluster_os.Rdata")
load("./outdata/step2_count_log2.Rdata")
load("./outdata/step2_tpm_log2.Rdata")
load("./indata/GSVA_out.Rdata")
.libPaths("/home/cmx/R/x86_64-pc-linux-gnu-library/4.2",include.site = TRUE)
library(limma)
library(clusterProfiler)
library(estimate)
library(GDCRNATools)
library(TCGAbiolinks)
library(ComplexHeatmap)
library(GSVA)
library(Rtsne)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(DESeq2)
library(edgeR)
# res1 <- nmf(PSI_data_top,3, 'Brunet', maxIter = 500,seed=211)
consensusmap(res1)

group1 <- predict(res1)
group1 <- as.data.frame(group1)
group1$group1 <- paste0('Cluster',group1$group1)
group1$sample <- rownames(group1)
group1<- group1[order(group1$group1),]
table(group1$group1)
Surv_4group<-merge(group1,clinical_age,by.x ="sample",by.y ="bcr_patient_barcode" ) 
clinical_OS<-Surv_4group
clinical_OS$OS.time<-ifelse(clinical_OS$vital_status=='Alive',clinical_OS$days_to_last_followup,clinical_OS$days_to_death)
clinical_OS$OS<- ifelse(clinical_OS$vital_status=="Alive",0,1)
clinical_OS$OS.time<-as.numeric(clinical_OS$OS.time)
#clinical_OS[167,8]<-0
clinical_OS<-clinical_OS[-grep("TCGA-B6-A0RQ", clinical_OS$sample), ]
# load("./outdata/step1_cluster_os.Rdata")
#BRCA_clin<-BRCA_clin%>%mutate(SA=ifelse(activ_sum>median(BRCA_clin$activ_sum),"High","Low"))
fit <- survfit(Surv(OS.time,OS) ~ group1, data = clinical_OS)
# save(clinical_OS,file = "./clinical_OS.Rdata")
cox <- coxph(Surv(OS.time,OS) ~ group1, data = clinical_OS)
coxSummary <- summary(cox)
sdiff <- pairwise_survdiff(Surv(OS.time,OS) ~ group1, data = clinical_OS)
sdiff

ggsurvplot(
  fit = fit,data = clinical_OS,#fit为通过survfit确定生存分析的对象，data 为生存分析所用的数据
  title="Elderly patients",
  legend=c(0.85,0.85),#更改图例的位置
  legend.title="Cluster",#删掉图例标题，或者定义图例标题为legend.title="SEX"
  # legend.labs=c("Cluster1 (n = 85)","Cluster2 (n = 156)","Cluster3 (n = 99)"),
  palette =
    c("#e07854", "#51b1b7","#CD3333"), ###  自定义颜色
  ggtheme = theme_survminer(),
  # surv.median.line = "hv",
  risk.table = F,
  axes.offset=TRUE,#设置使生存曲线从原点开始
  xlab ="Time(days)",#横轴描述
  break.y.by = 0.2,
  # xlim=c(0,4050),
  conf.int=FALSE,
  ylab ="Overall Survival",#纵轴描述
  size=1.5,##线条粗细
  censor = TRUE,##显示删失值
  pval = TRUE,
  # pval.size = 5,
  # pval.coord = c(11, 0.12)
)
