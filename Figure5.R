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
