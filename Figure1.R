setwd('/home/cmx/senes_splicing/')
rm(list=ls())
library(dbplyr)
library(clusterProfiler)
library(dplyr)
library(TCGAbiolinks)
library(DMwR2)
library(XML)
library(limma)
library(Rtsne)
library(ggplot2)
library(NMF)
library(tidyverse)
library(ConsensusClusterPlus)
library(survival)
library(survminer)
.libPaths("/home/cmx/R/x86_64-pc-linux-gnu-library/4.2",include.site = TRUE)
load("./outdata/step1_cluster_os.Rdata")
###整理一下这些剪接事件的患者年龄信息####
###整理可变剪接数据####
# 下载的时候：a percentage of samples with a PSI value >=75 was set as the cutoff value
PSI_data <- read.table(file="/home/cmx/My_Project/indata/TCGA_cancer/BRCA/1_raw_download_data4step1/PSI_download_BRCA.txt",header=TRUE,check.names=F,quote="")
# PSI_range0.75以下的,认为不是那么重要的可变剪接事件
PSI_data <- PSI_data[PSI_data$psi_range >= 0.75,]
colnames(PSI_data)[11:1217] <- gsub("_", "-", colnames(PSI_data)[11:1217])
# 不去除正常样本
#PSI_data <- subset(PSI_data, select = !grepl("-Norm", colnames(PSI_data)))#得到1104个样本
# mean value of PSI less than 0.01 or standard deviations value of PSI less than 0.01 were ruled out.另一篇文献有说
PSI_data <- PSI_data[PSI_data$std_psi >= 0.01 ,]


####整理一下临床数据####
setwd('/home/cmx/Others_Analysis/KeH/InputData/')
length(dir("/home/cmx/Others_Analysis/KeH/InputData/clinical/"))
#[1] 537
length(dir("/home/cmx/Others_Analysis/KeH/InputData/expdata/"))
#[1] 613
library(XML)
xmls = dir("clinical/",pattern = "*.xml$",recursive = T)
cl = list()
for(i in 1:length(xmls)){
  result = xmlParse(paste0("clinical/",xmls[[i]]))
  rootnode = xmlRoot(result)
  cl[[i]] = xmlToDataFrame(rootnode[2])
}
setwd('/home/cmx/senes_splicing/')
clinical = do.call(rbind,cl)
clinical[1:3,1:3]
clinical_age<-clinical[,c(12,22,9,6,10)]
# save(clinical_age,file="./outdata/clinical_age.Rdata")
clinical_young<-clinical_age[clinical_age$age_at_initial_pathologic_diagnosis<=45,]
clinical_old<-clinical_age[clinical_age$age_at_initial_pathologic_diagnosis>=65,]
####年轻正常样本有多少####
Sample_young_normol<-merge(clinical_young,PSI_data_normal,by.x = "bcr_patient_barcode",by.y="Sample_ID")
# 如果Sample_ID不包含"-Norm"，则加上
Sample_young_normol$bcr_patient_barcode <- ifelse(!grepl("-Norm", Sample_young_normol$bcr_patient_barcode), 
                                                  paste0(Sample_young_normol$bcr_patient_barcode, "-Norm"), Sample_young_normol$bcr_patient_barcode)
Sample_young_normol<-Sample_young_normol[,-c(2:5)]
####年轻癌症样本有多少####
Sample_young_cancer<-merge(clinical_young,PSI_data_cancer,by.x = "bcr_patient_barcode",by.y="Sample_ID")
Sample_young_cancer<-Sample_young_cancer[,-c(2:5)]
####老年正常样本有多少####
Sample_old_normol<-merge(clinical_old,PSI_data_normal,by.x = "bcr_patient_barcode",by.y="Sample_ID")
Sample_old_normol$bcr_patient_barcode <- ifelse(!grepl("-Norm", Sample_old_normol$bcr_patient_barcode), 
                                                paste0(Sample_old_normol$bcr_patient_barcode, "-Norm"), Sample_old_normol$bcr_patient_barcode)
Sample_old_normol<-Sample_old_normol[,-c(2:5)]
####老年癌症样本有多少####
Sample_old_cancer<-merge(clinical_old,PSI_data_cancer,by.x = "bcr_patient_barcode",by.y="Sample_ID")
Sample_old_cancer<-Sample_old_cancer[,-c(2:5)]