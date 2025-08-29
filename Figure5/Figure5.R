setwd('/home/cmx/senes_splicing/')
rm(list=ls())
library(estimate)
library(GDCRNATools)
library(TCGAbiolinks)
library(ComplexHeatmap)
library(GSVA)
library(Rtsne)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(Hmisc)
library(survival)
library(survminer)
load("./outdata/step1_cluster_os.Rdata")
load("./outdata/step2_count_log2.Rdata")
load("./outdata/step2_tpm_log2.Rdata")
load("./outdata/old_patient_stage.Rdata")

PSI_data_group<-PSI_data_top[,clinical_OS$sample]
PSI_data_group <- PSI_data_group[, order(match(colnames(PSI_data_group), old_patient_stage$Sample_ID))]
####扰动程度评级####
# # 下载的时候：a percentage of samples with a PSI value >=75 was set as the cutoff value
# PSI_data <- read.table(file="/home/cmx/My_Project/indata/PSI_download_BRCA.txt",header=TRUE,check.names=F,quote="")
# # PSI_range0.75以下的,认为不是那么重要的可变剪接事件
# PSI_data <- PSI_data[PSI_data$psi_range >= 0.75,]
# colnames(PSI_data)[11:1217] <- gsub("_", "-", colnames(PSI_data)[11:1217])
# # 不去除正常样本
# #PSI_data <- subset(PSI_data, select = !grepl("-Norm", colnames(PSI_data)))#得到1104个样本
# # mean value of PSI less than 0.01 or standard deviations value of PSI less than 0.01 were ruled out.另一篇文献有说
# PSI_data <- PSI_data[PSI_data$std_psi >= 0.01 ,]
# # 其实起始点记录不是那么重要，还是先不删除了
# # PSI_data <- PSI_data[!(PSI_data$from_exon == "null" & PSI_data$to_exon == "null"),]
# 
# PSI_data[,11:1217] <- PSI_data[,11:1217] %>%
#   mutate(across(where(is.character), as.numeric))
# # 使用K近邻（KNN）算法填充PDUI或PSI值的缺失值
# #PSI_data_filled <- knnImputation(PSI_data[11:1217], k=10)#跑太久了，不想再跑了
# #write.csv(PSI_data_filled,"./outdata/PSI_data_filled.csv")
# 
# PSI_data_filled<-read.csv("./outdata/PSI_data_filled.csv",header = TRUE,sep = ",")
# 
# PSI_data_filled<-PSI_data_filled[,-1]
# colnames(PSI_data_filled) <- gsub("\\.","-", colnames(PSI_data_filled))
# 
# PSI_data_data<-PSI_data[,1:10]
# PSI_data<-cbind(PSI_data_data,PSI_data_filled)
# 
# # 创建新的行名，并删除空格
# Splice_event <- apply(PSI_data[, 1:3], 1, function(x) gsub("\\s", "", paste(x, collapse = "-")))
# 
# # 将新的行名赋值给数据框的行名
# rownames(PSI_data) <- Splice_event
# 
# PSI_data<-PSI_data[,c(11:1217)]
# PSI_data<-as.data.frame(t(PSI_data))
# PSI_data$Sample_ID<-row.names(PSI_data)
# 
# PSI_data <- PSI_data %>%
#   dplyr::select(Sample_ID, everything())
# 
# colnames(PSI_data) <- gsub(" ", "", colnames(PSI_data))
# 
# col_names <- rownames(PSI_data)
# 
# # 按照包含"-Norm"的列排序
# new_order <- c(
#   col_names[grepl("-Norm", col_names)],
#   col_names[!grepl("-Norm", col_names)]
# )
# 
# # 使用order函数重新排列列
# PSI_data <- PSI_data[new_order, ]
# PSI_data_normal<-PSI_data[1:113,]
# 
# PSI_data_cancer<-PSI_data[114:1207,]
# PSI_data_normal<-as.data.frame(t(PSI_data_normal))
# PSI_data_normal<-PSI_data_normal[-1,]
# PSI_data_cancer<-as.data.frame(t(PSI_data_cancer))
# PSI_data_cancer<-PSI_data_cancer[-1,]
# 
# PSI_data_group1<-PSI_data_group
# PSI_data_group1$as_id<-row.names(PSI_data_group1)
# 
# PSI_data_normal <- PSI_data_normal[rownames(PSI_data_normal) %in% PSI_data_group1$as_id, ]
# PSI_data_normal <- PSI_data_normal %>%
#   mutate(across(where(is.character), as.numeric))
# PSI_data_normal$PSI_NOR<-rowSums(PSI_data_normal)/113
# PSI_data_cancer<-PSI_data_cancer[rownames(PSI_data_cancer) %in% PSI_data_group1$as_id, ]
# PSI_data_cancer <- PSI_data_cancer %>%
#   mutate(across(where(is.character), as.numeric))

# save(PSI_data_cancer,PSI_data_normal,file="./outdata/step3_perturbation_level.Rdata")
load("./outdata/step3_perturbation_level.Rdata")
BRCA_perturbation<-PSI_data_cancer-PSI_data_normal$PSI_NOR
BRCA_perturbation<-as.data.frame(t(BRCA_perturbation))
BRCA_perturbation <- BRCA_perturbation[rownames(BRCA_perturbation) %in% clinical_OS$sample, ]
BRCA_perturbation<-as.data.frame(t(BRCA_perturbation))

BRCA_perturbation$acti<-rowSums(BRCA_perturbation)/340
BRCA_perturbation$acti_level<-ifelse(BRCA_perturbation$acti > 0,"act","res")

BRCA_perturbation_acti<-BRCA_perturbation[BRCA_perturbation$acti_level=='act',]
BRCA_perturbation_acti$event<-rownames(BRCA_perturbation_acti)
BRCA_perturbation_res<-BRCA_perturbation[BRCA_perturbation$acti_level=='res',]
BRCA_perturbation_res$event<-rownames(BRCA_perturbation_res)

AS_acti<-BRCA_perturbation_acti[,341:343]
AS_res<-BRCA_perturbation_res[,341:343]

AS_acti_matrix<-PSI_data_cancer[AS_acti$event,]
AS_acti_matrix <- AS_acti_matrix[,colnames(AS_acti_matrix) %in% clinical_OS$sample ]
AS_res_matrix<-PSI_data_cancer[AS_res$event,]
AS_res_matrix <- AS_res_matrix[, colnames(AS_res_matrix) %in% clinical_OS$sample]

Ontology <- read_csv("/home/cmx/senes_splicing/indata//papers_Ontology.csv")#RBP list
genes<-Ontology$HGNC_symbol
TCGA_KIRC1<-old_patient_count_log[genes,]
TCGA_KIRC1<-na.omit(TCGA_KIRC1)
TCGA_KIRC1<-as.data.frame(t(TCGA_KIRC1))
AS_acti_matrix<-as.data.frame(t(AS_acti_matrix))
AS_res_matrix<-as.data.frame(t(AS_res_matrix))

####找出正向调控者####
TCGA_KIRC1 <- TCGA_KIRC1[order(match(rownames(TCGA_KIRC1), rownames(AS_acti_matrix))), ]

All_corr<-rcorr(as.matrix(TCGA_KIRC1),as.matrix(AS_acti_matrix),type = "spearman")##总的相关性

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

flattenCorrMatrix(All_corr$r,All_corr$P)
CHOL_All_cor<-flattenCorrMatrix(All_corr$r,All_corr$P)
fdr=p.adjust(as.numeric(as.vector(CHOL_All_cor$p)),method = "BH")
acti_All_cor<-cbind(CHOL_All_cor,fdr)
acti_All_cor <- acti_All_cor[grepl("-AA|-AD|-AP|-AT|-ES|-ME|-RI", acti_All_cor$column), ]
acti_All_cor <- acti_All_cor[!grepl("-AA|-AD|-AP|-AT|-ES|-ME|-RI", acti_All_cor$row), ]
acti_All_cor<-acti_All_cor[acti_All_cor$fdr < 0.1,]
acti_All_cor<-acti_All_cor[acti_All_cor$cor > 0.3,]
h_acti<-as.data.frame(table(acti_All_cor$row))
h_acti_30<-h_acti[h_acti$Freq >= 30,]

####找出负向调控者####
TCGA_KIRC1 <- TCGA_KIRC1[order(match(rownames(TCGA_KIRC1), rownames(AS_res_matrix))), ]

All_corr<-rcorr(as.matrix(TCGA_KIRC1),as.matrix(AS_res_matrix),type = "spearman")##总的相关性

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

flattenCorrMatrix(All_corr$r,All_corr$P)
CHOL_All_cor<-flattenCorrMatrix(All_corr$r,All_corr$P)
fdr=p.adjust(as.numeric(as.vector(CHOL_All_cor$p)),method = "BH")
res_All_cor<-cbind(CHOL_All_cor,fdr)
res_All_cor <- res_All_cor[grepl("-AA|-AD|-AP|-AT|-ES|-ME|-RI", res_All_cor$column), ]
res_All_cor <- res_All_cor[!grepl("-AA|-AD|-AP|-AT|-ES|-ME|-RI", res_All_cor$row), ]
res_All_cor<-res_All_cor[res_All_cor$fdr < 0.1,]
res_All_cor<-res_All_cor[res_All_cor$cor > 0.3,]
h_res<-as.data.frame(table(res_All_cor$row))
h_res_30<-h_res[h_res$Freq >= 30,]
#有没有交集
res_only<-setdiff(h_res_30$Var1,h_acti_30$Var1)
acti_only<-setdiff(h_acti_30$Var1,h_res_30$Var1)
res_act_RBP<-union(h_res_30$Var1,h_acti_30$Var1)
resWact_RBP<-intersect(h_res_30$Var1,h_acti_30$Var1)
save(h_res_30,h_acti_30,file = "./outdata/h_30.Rdata")