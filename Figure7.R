setwd('/home/cmx/senes_splicing/')
rm(list=ls())
load("/home/cmx/senes_splicing/outdata/step1_cluster_os.Rdata")
load("/home/cmx/senes_splicing/outdata/step2_count_log2.Rdata")
library(oncoPredict)
library(rGREAT)
library(niboR)
library(tidyverse)
#####CTRP####
dir='/home/cmx/Reproduce_and_become_mine/BLCA_Oncogene/InputData/DataFiles/Training Data'
CTRP1_Res = readRDS(file = file.path(dir,"CTRP2_Res.rds"))
CTRP1_Res <- exp(CTRP1_Res)
CTRP1_Expr = readRDS(file=file.path(dir,'CTRP2_Expr (TPM, not log transformed).rds'))
CTRP1_Expr<-log2(CTRP1_Expr+1)
th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
tumor_logcount<-old_patient_count_log[,clinical_OS$sample]
tumor_logcount<-as.data.frame(tumor_logcount)
expr_batch<-tumor_logcount
testExpr <- expr_batch####运算很慢，这还只是十个数据
calcPhenotype(trainingExprData = CTRP1_Expr,
              trainingPtype = CTRP1_Res,
              testExprData = as.matrix(testExpr),
              batchCorrect = 'eb',  #   "eb" for ComBat
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10,
              printOutput = TRUE,
              removeLowVaringGenesFrom = 'rawData' )