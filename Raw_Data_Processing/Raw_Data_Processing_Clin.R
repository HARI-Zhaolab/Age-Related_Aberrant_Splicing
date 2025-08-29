setwd('/home/cmx/senes_splicing/')
.libPaths("/home/cmx/R/x86_64-pc-linux-gnu-library/4.2",include.site = TRUE)
rm(list=ls())
library(dbplyr)
library(clusterProfiler)
library(dplyr)
library(TCGAbiolinks)
library(DMwR2)
library(limma)
library(Rtsne)
library(ggplot2)
library(NMF)
library(tidyverse)
library(survival)
library(survminer)

setwd('/home/cmx/senes_splicing/')
load("./clinical_age.Rdata")
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

####整理最终的四个grop的PSI矩阵
Sample_young_normol<-as.data.frame(t(Sample_young_normol))
colnames(Sample_young_normol)<-Sample_young_normol[1,]
Sample_young_normol<-Sample_young_normol[-1,]

Sample_young_cancer<-as.data.frame(t(Sample_young_cancer))
colnames(Sample_young_cancer)<-Sample_young_cancer[1,]
Sample_young_cancer<-Sample_young_cancer[-1,]


Sample_old_normol<-as.data.frame(t(Sample_old_normol))
colnames(Sample_old_normol)<-Sample_old_normol[1,]
Sample_old_normol<-Sample_old_normol[-1,]

Sample_old_cancer<-as.data.frame(t(Sample_old_cancer))
colnames(Sample_old_cancer)<-Sample_old_cancer[1,]
Sample_old_cancer<-Sample_old_cancer[-1,]
