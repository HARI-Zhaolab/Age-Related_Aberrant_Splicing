setwd('/home/cmx/senes_splicing/')
rm(list=ls())
load("/home/cmx/senes_splicing/outdata/step1_cluster_os.Rdata")
load("/home/cmx/senes_splicing/outdata/step2_count_log2.Rdata")
library(oncoPredict)
library(rGREAT)
library(niboR)
library(tidyverse)
# #####整合好GDSC####
# dir='/home/cmx/Reproduce_and_become_mine/BLCA_Oncogene/InputData/DataFiles/Training Data/'
# GDSC1_Res = readRDS(file = file.path(dir,"GDSC1_Res.rds"))
# GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
# GDSC_Res_list<-intersect(rownames(GDSC1_Res),rownames(GDSC2_Res))
# GDSC1_Res <-GDSC1_Res[GDSC_Res_list,]
# GDSC2_Res <-GDSC1_Res[GDSC_Res_list,]
# GDSC_Res<-cbind(GDSC1_Res,GDSC2_Res)
# GDSC_Res <- exp(GDSC_Res)
# 
# GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
# GDSC2_Expr<-GDSC2_Expr[,GDSC_Res_list]
# GDSC1_Expr = readRDS(file=file.path(dir,'GDSC1_Expr (RMA Normalized and Log Transformed).rds'))
# GDSC1_Expr<-GDSC1_Expr[,GDSC_Res_list]
# GDSC_Expr<-GDSC1_Expr
# 
# 
# th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
# tumor_logcount<-old_patient_count_log[,clinical_OS$sample]
# tumor_logcount<-as.data.frame(tumor_logcount)
# expr_batch<-tumor_logcount
# testExpr <- expr_batch####运算很慢，这还只是十个数据
# calcPhenotype(trainingExprData = GDSC_Expr,
#               trainingPtype = GDSC_Res,
#               testExprData = as.matrix(testExpr),
#               batchCorrect = 'eb',  #   "eb" for ComBat
#               powerTransformPhenotype = TRUE,
#               removeLowVaryingGenes = 0.2,
#               minNumSamples = 10,
#               printOutput = TRUE,
#               removeLowVaringGenesFrom = 'rawData' )
# 
#####CTRP####
# dir='/home/cmx/Reproduce_and_become_mine/BLCA_Oncogene/InputData/DataFiles/Training Data'
# CTRP1_Res = readRDS(file = file.path(dir,"CTRP2_Res.rds"))
# CTRP1_Res <- exp(CTRP1_Res)
# CTRP1_Expr = readRDS(file=file.path(dir,'CTRP2_Expr (TPM, not log transformed).rds'))
# CTRP1_Expr<-log2(CTRP1_Expr+1)
# th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
# tumor_logcount<-old_patient_count_log[,clinical_OS$sample]
# tumor_logcount<-as.data.frame(tumor_logcount)
# expr_batch<-tumor_logcount
# testExpr <- expr_batch####运算很慢，这还只是十个数据
# calcPhenotype(trainingExprData = CTRP1_Expr,
#               trainingPtype = CTRP1_Res,
#               testExprData = as.matrix(testExpr),
#               batchCorrect = 'eb',  #   "eb" for ComBat
#               powerTransformPhenotype = TRUE,
#               removeLowVaryingGenes = 0.2,
#               minNumSamples = 10,
#               printOutput = TRUE,
#               removeLowVaringGenesFrom = 'rawData' )

####整理CCLE数据库####
# # 读取 GCT 文件
# CCLE_Expr <- read.gct("/home/cmx/Reproduce_and_become_mine/BLCA_Oncogene/InputData/DataFiles/Training Data/CCLE_RNAseq_genes_counts_20180929.gct")
# CCLE_Expr<-as.data.frame(CCLE_Expr)
# rownames(CCLE_Expr) <- substr(rownames(CCLE_Expr),1,15)
# ##处理count
# CCLE_Expr$ENSEMBL<-rownames(CCLE_Expr)
# 
# # for test
# ENSEMBL <-CCLE_Expr$ENSEMBL
# df <- data.frame(ENSEMBL)
# 
# # ENSG转换Symbol
# df$symbol <- mapIds(org.Hs.eg.db,
#                     keys=ENSEMBL,
#                     column="SYMBOL",
#                     keytype="ENSEMBL",
#                     multiVals="first")
# CCLE_Expr<-merge(CCLE_Expr,df,by="ENSEMBL")
# CCLE_Expr <- CCLE_Expr[!is.na(CCLE_Expr$symbol), ]
# #将gene_name列去除重复的基因，保留每个基因最大表达量结果
# CCLE_Expr <- aggregate( . ~ symbol,data=CCLE_Expr, max)
# rownames(CCLE_Expr)<-CCLE_Expr$symbol
# CCLE_Expr<-CCLE_Expr[,-1]
# CCLE_Expr<-CCLE_Expr[,-1]
# CCLE_Expr <- CCLE_Expr %>%
#   dplyr::mutate(across(where(is.character), as.numeric))
# CCLE_Expr<-CCLE_Expr[!rowSums(CCLE_Expr)==0,]
# CCLE_Expr<-log2(CCLE_Expr+1)
# CCLE_Res<-read.csv2("/home/cmx/Reproduce_and_become_mine/BLCA_Oncogene/InputData/DataFiles/Training Data/CCLE_NP24.2009_Drug_data_2015.02.24.csv",
#                     sep = ",", check.names=F)
# CCLE_Res<-CCLE_Res[,c(1,3,11)]
# # 使用 spread 转换为矩阵
# CCLE_Matrix <- CCLE_Res %>%
#   spread(Compound, `IC50 (uM)`)
# 
# # 将行名设置为 'Primary Cell Line Name'
# rownames(CCLE_Matrix) <- CCLE_Matrix$`CCLE Cell Line Name`
# CCLE_Matrix<-CCLE_Matrix[,-1]
# CCLE_Res<-CCLE_Matrix
# CCLE_Res <- CCLE_Res %>%
#   dplyr::mutate(across(where(is.character), as.numeric))
# CCLE_Expr<-CCLE_Expr[,colnames(CCLE_Expr) %in% rownames(CCLE_Res)]
# CCLE_Res<-CCLE_Res[rownames(CCLE_Res) %in% colnames(CCLE_Expr),]
# CCLE_Res <- exp(CCLE_Res)
# 
# CCLE_Expr<-as.matrix(CCLE_Expr)
# CCLE_Res<-as.matrix(CCLE_Res)
# save(CCLE_Res,CCLE_Expr,file = "/home/cmx/Reproduce_and_become_mine/BLCA_Oncogene/InputData/DataFiles/Training Data/CCLE_EXP_AND_DRUG.Rdata")
# load("/home/cmx/Reproduce_and_become_mine/BLCA_Oncogene/InputData/DataFiles/Training Data/CCLE_EXP_AND_DRUG.Rdata")
# load("/home/cmx/senes_splicing/outdata/step1_cluster_os.Rdata")
# th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
# tumor_logcount<-old_patient_count_log[,clinical_OS$sample]
# tumor_logcount<-as.data.frame(tumor_logcount)
# expr_batch<-tumor_logcount
# testExpr <- expr_batch####运算很慢，这还只是十个数据
# calcPhenotype(trainingExprData = CCLE_Expr,
#               trainingPtype = CCLE_Res,
#               testExprData = as.matrix(testExpr),
#               batchCorrect = 'eb',  #   "eb" for ComBat
#               powerTransformPhenotype = TRUE,
#               removeLowVaryingGenes = 0.2,
#               minNumSamples = 10,
#               printOutput = TRUE,
#               removeLowVaringGenesFrom = 'rawData' )
# CCLE_predict<-read.csv2("./calcPhenotype_Output/DrugPredictions.csv",sep = ",")
# save(CCLE_predict, file = "./outdata/CCLE_predict.Rdata")
# load("./outdata/CCLE_predict.Rdata")
# # 
# # 
# # GDSC_predict<-read.csv2("./calcPhenotype_Output/DrugPredictions.csv",sep = ",")
# # save(GDSC_predict, file = "./outdata/GDSC_predict.Rdata")
# load("./outdata/GDSC_predict.Rdata")
# 
# # CTRP_predict<-read.csv2("./calcPhenotype_Output/DrugPredictions.csv",sep = ",")
# # save(CTRP_predict, file = "./outdata/GTRP_predict.Rdata")
# load("./outdata/GTRP_predict.Rdata")
# 
# ####Metabric整合好GDSC####
# load("/home/cmx/Others_Analysis/KeH/InputData/METABRIC-master/data/metabric_clinical.Rdata")
# load("/home/cmx/Others_Analysis/KeH/InputData/METABRIC-master/data/metabric_expression.Rdata")
# # 将点号替换为破折号 
# expr<-as.data.frame(expr)
# colnames(expr) <- sub("\\.", "-", colnames(expr))
# clin_TAR_NEED<-clin[,c(1,2,3,7,16)]
# clin_TAR_NEED<-clin_TAR_NEED[clin_TAR_NEED$AGE_AT_DIAGNOSIS>=65,]#只要65的
# 
# expr<-as.data.frame(t(expr))
# expr$sample<-rownames(expr)
# expr = dplyr::select(expr,24361,everything())
# expr <- expr[expr$sample %in% clin_TAR_NEED$PATIENT_ID,]
# 
# dir='/home/cmx/Reproduce_and_become_mine/BLCA_Oncogene/InputData/DataFiles/Training Data/'
# GDSC1_Res = readRDS(file = file.path(dir,"GDSC1_Res.rds"))
# GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
# GDSC_Res_list<-intersect(rownames(GDSC1_Res),rownames(GDSC2_Res))
# GDSC1_Res <-GDSC1_Res[GDSC_Res_list,]
# GDSC2_Res <-GDSC1_Res[GDSC_Res_list,]
# GDSC_Res<-cbind(GDSC1_Res,GDSC2_Res)
# GDSC_Res <- exp(GDSC_Res)
# 
# GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
# GDSC2_Expr<-GDSC2_Expr[,GDSC_Res_list]
# GDSC1_Expr = readRDS(file=file.path(dir,'GDSC1_Expr (RMA Normalized and Log Transformed).rds'))
# GDSC1_Expr<-GDSC1_Expr[,GDSC_Res_list]
# GDSC_Expr<-GDSC1_Expr
# 
# th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
# testExpr <- expr####运算很慢，这还只是十个数据
# calcPhenotype(trainingExprData = GDSC_Expr,
#               trainingPtype = GDSC_Res,
#               testExprData = as.matrix(expr),
#               batchCorrect = 'eb',  #   "eb" for ComBat
#               powerTransformPhenotype = TRUE,
#               removeLowVaryingGenes = 0.2,
#               minNumSamples = 10,
#               printOutput = TRUE,
#               removeLowVaringGenesFrom = 'rawData' )
# ####Metabric整合好CTRP####
# load("/home/cmx/Others_Analysis/KeH/InputData/METABRIC-master/data/metabric_clinical.Rdata")
# load("/home/cmx/Others_Analysis/KeH/InputData/METABRIC-master/data/metabric_expression.Rdata")
# library(oncoPredict)
# library(rGREAT)
# library(niboR)
# library(tidyverse)
# # 将点号替换为破折号 
# expr<-as.data.frame(expr)
# colnames(expr) <- sub("\\.", "-", colnames(expr))
# clin_TAR_NEED<-clin[,c(1,2,3,7,16)]
# clin_TAR_NEED<-clin_TAR_NEED[clin_TAR_NEED$AGE_AT_DIAGNOSIS>=65,]#只要65的
# 
# expr<-as.data.frame(t(expr))
# expr$sample<-rownames(expr)
# expr = dplyr::select(expr,24361,everything())
# expr <- expr[expr$sample %in% clin_TAR_NEED$PATIENT_ID,]
# expr<-expr[,-1]
# expr<-as.data.frame(t(expr))
# dir='/home/cmx/Reproduce_and_become_mine/BLCA_Oncogene/InputData/DataFiles/Training Data'
# CTRP1_Res = readRDS(file = file.path(dir,"CTRP2_Res.rds"))
# CTRP1_Res <- exp(CTRP1_Res)
# CTRP1_Expr = readRDS(file=file.path(dir,'CTRP2_Expr (TPM, not log transformed).rds'))
# CTRP1_Expr<-log2(CTRP1_Expr+1)
# th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
# testExpr <- expr####运算很慢，这还只是十个数据
# calcPhenotype(trainingExprData = CTRP1_Expr,
#               trainingPtype = CTRP1_Res,
#               testExprData = as.matrix(testExpr),
#               batchCorrect = 'eb',  #   "eb" for ComBat
#               powerTransformPhenotype = TRUE,
#               removeLowVaryingGenes = 0.2,
#               minNumSamples = 10,
#               printOutput = TRUE,
#               removeLowVaringGenesFrom = 'rawData' )
# 
# 
# ####Metabric整理CCLE数据库####
# # # 读取 GCT 文件
# # CCLE_Expr <- read.gct("/home/cmx/Reproduce_and_become_mine/BLCA_Oncogene/InputData/DataFiles/Training Data/CCLE_RNAseq_genes_counts_20180929.gct")
# # CCLE_Expr<-as.data.frame(CCLE_Expr)
# # rownames(CCLE_Expr) <- substr(rownames(CCLE_Expr),1,15)
# # ##处理count
# # CCLE_Expr$ENSEMBL<-rownames(CCLE_Expr)
# # 
# # # for test
# # ENSEMBL <-CCLE_Expr$ENSEMBL
# # df <- data.frame(ENSEMBL)
# # 
# # # ENSG转换Symbol
# # df$symbol <- mapIds(org.Hs.eg.db,
# #                     keys=ENSEMBL,
# #                     column="SYMBOL",
# #                     keytype="ENSEMBL",
# #                     multiVals="first")
# # CCLE_Expr<-merge(CCLE_Expr,df,by="ENSEMBL")
# # CCLE_Expr <- CCLE_Expr[!is.na(CCLE_Expr$symbol), ]
# # #将gene_name列去除重复的基因，保留每个基因最大表达量结果
# # CCLE_Expr <- aggregate( . ~ symbol,data=CCLE_Expr, max)
# # rownames(CCLE_Expr)<-CCLE_Expr$symbol
# # CCLE_Expr<-CCLE_Expr[,-1]
# # CCLE_Expr<-CCLE_Expr[,-1]
# # CCLE_Expr <- CCLE_Expr %>%
# #   dplyr::mutate(across(where(is.character), as.numeric))
# # CCLE_Expr<-CCLE_Expr[!rowSums(CCLE_Expr)==0,]
# # CCLE_Expr<-log2(CCLE_Expr+1)
# # CCLE_Res<-read.csv2("/home/cmx/Reproduce_and_become_mine/BLCA_Oncogene/InputData/DataFiles/Training Data/CCLE_NP24.2009_Drug_data_2015.02.24.csv",
# #                     sep = ",", check.names=F)
# # CCLE_Res<-CCLE_Res[,c(1,3,11)]
# # # 使用 spread 转换为矩阵
# # CCLE_Matrix <- CCLE_Res %>%
# #   spread(Compound, `IC50 (uM)`)
# # 
# # # 将行名设置为 'Primary Cell Line Name'
# # rownames(CCLE_Matrix) <- CCLE_Matrix$`CCLE Cell Line Name`
# # CCLE_Matrix<-CCLE_Matrix[,-1]
# # CCLE_Res<-CCLE_Matrix
# # CCLE_Res <- CCLE_Res %>%
# #   dplyr::mutate(across(where(is.character), as.numeric))
# # CCLE_Expr<-CCLE_Expr[,colnames(CCLE_Expr) %in% rownames(CCLE_Res)]
# # CCLE_Res<-CCLE_Res[rownames(CCLE_Res) %in% colnames(CCLE_Expr),]
# # CCLE_Res <- exp(CCLE_Res)
# # 
# # CCLE_Expr<-as.matrix(CCLE_Expr)
# # CCLE_Res<-as.matrix(CCLE_Res)
# # save(CCLE_Res,CCLE_Expr,file = "/home/cmx/Reproduce_and_become_mine/BLCA_Oncogene/InputData/DataFiles/Training Data/CCLE_EXP_AND_DRUG.Rdata")
# # load("/home/cmx/Reproduce_and_become_mine/BLCA_Oncogene/InputData/DataFiles/Training Data/CCLE_EXP_AND_DRUG.Rdata")
# # load("/home/cmx/senes_splicing/outdata/step1_cluster_os.Rdata")
# # th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
# # testExpr <- expr####运算很慢，这还只是十个数据
# # calcPhenotype(trainingExprData = CCLE_Expr,
# #               trainingPtype = CCLE_Res,
# #               testExprData = as.matrix(testExpr),
# #               batchCorrect = 'eb',  #   "eb" for ComBat
# #               powerTransformPhenotype = TRUE,
# #               removeLowVaryingGenes = 0.2,
# #               minNumSamples = 10,
# #               printOutput = TRUE,
# #               removeLowVaringGenesFrom = 'rawData' )
# # CCLE_predict_metabric<-read.csv2("./calcPhenotype_Output/DrugPredictions.csv",sep = ",")
# # save(CCLE_predict_metabric, file = "./outdata/CCLE_predict_metabric.Rdata")
# load("./outdata/CCLE_predict_metabric.Rdata")
# 
# # GDSC_predict_metabric<-read.csv2("./calcPhenotype_Output/DrugPredictions.csv",sep = ",")
# # save(GDSC_predict_metabric, file = "./outdata/GDSC_predict_metabric.Rdata")
# load("./outdata/GDSC_predict_metabric.Rdata")
# 
# # CTRP_predict_metabric<-read.csv2("./calcPhenotype_Output/DrugPredictions.csv",sep = ",")
# # save(CTRP_predict_metabric, file = "./outdata/GTRP_predict_metabric.Rdata")
# load("./outdata/GTRP_predict_metabric.Rdata")
# 
# 

####加载后台运行得到的所有结果####
load("/home/cmx/senes_splicing/outdata/step1_cluster_os.Rdata")
load("/home/cmx/senes_splicing/outdata/step2_count_log2.Rdata")
load("./outdata/CCLE_predict.Rdata")
load("./outdata/GDSC_predict.Rdata")
load("./outdata/GTRP_predict.Rdata")
load("./outdata/CCLE_predict_metabric.Rdata")
load("./outdata/GDSC_predict_metabric.Rdata")
load("./outdata/GTRP_predict_metabric.Rdata")
####TCGA-GDSC####
rownames(GDSC_predict)<-GDSC_predict$X
GDSC_predict<-GDSC_predict[,-1]
GDSC_predict<-as.data.frame(t(GDSC_predict))
GDSC_predict <- GDSC_predict %>%
  dplyr::mutate(across(where(is.character), as.numeric))
GDSC_predict<-as.matrix(t(scale(t(GDSC_predict))))

TCGA_KIRC1<-GDSC_predict

p_value <- data.frame(gene = rownames(data.frame(TCGA_KIRC1)))

TCGA_KIRC1 <- as.data.frame(t(TCGA_KIRC1))
TCGA_KIRC1$Sample_ID <- rownames(TCGA_KIRC1)
TCGA_KIRC1<-merge(TCGA_KIRC1,old_patient_stage,by="Sample_ID")
colnames(TCGA_KIRC1)[747]<-"Cluster"
#HNRNPA3,HNRNPH2,HNRNPH3,RBM15,RBM23,RBM43,RBM7,RBMS2,SRSF1
TCGA_KIRC1$Cluster <- factor(TCGA_KIRC1$Cluster,levels = c("Cluster1","Cluster2","Cluster3"))

anno_pvalue <- data.frame()
for(i in p_value$gene){
  data <- dplyr::select(TCGA_KIRC1,gene = i,Cluster = Cluster)
  out_results <- compare_means(gene~Cluster ,data = data, method ="kruskal.test")
  out_results$.y. <- i
  anno_pvalue <- rbind(anno_pvalue,out_results)
}

##重新画图
anno_pvalue<-anno_pvalue[!anno_pvalue$p.signif=="ns",]
TCGA_KIRC1<-GDSC_predict[anno_pvalue$.y.,]
TCGA_KIRC1<-as.matrix(t(scale(t(TCGA_KIRC1))))
TCGA_KIRC1<-as.matrix(t(scale(t(TCGA_KIRC1))))

# TCGA_KIRC1[TCGA_KIRC1 > 2] <- 2
# TCGA_KIRC1[TCGA_KIRC1 < -2] <- -2

TCGA_KIRC1 <- TCGA_KIRC1[, order(match(colnames(TCGA_KIRC1), old_patient_stage$Sample_ID))]
sfd = HeatmapAnnotation(Subtype = old_patient_stage$group1)

row_ha2 = rowAnnotation(foo = anno_text(anno_pvalue$p.signif  , location = 0.5, just = "center"))

# 绘制热图
library(circlize)
ComplexHeatmap::Heatmap(TCGA_KIRC1,
                        col =colorRamp2(c(-4,0,4), c("#4DBBD5FF", "#FFFFFF", "#ED0000FF")),
                        column_split  = factor(old_patient_stage$group1, levels = c("Cluster1","Cluster2","Cluster3")),
                        name = " ",
                        row_names_side = "left",
                        right_annotation = row_ha2,
                        top_annotation = sfd,
                        show_column_names = F,
                        cluster_rows = F,
                        cluster_columns = F)

rownames(GDSC_predict)<-GDSC_predict$X
GDSC_predict<-GDSC_predict[,-1]
GDSC_predict<-as.data.frame(t(GDSC_predict))
GDSC_predict <- GDSC_predict %>%
  dplyr::mutate(across(where(is.character), as.numeric))
GDSC_predict<-as.matrix(t(scale(t(GDSC_predict))))

TCGA_KIRC1<-GDSC_predict

p_value <- data.frame(gene = rownames(data.frame(TCGA_KIRC1)))

TCGA_KIRC1 <- as.data.frame(t(TCGA_KIRC1))
TCGA_KIRC1$Sample_ID <- rownames(TCGA_KIRC1)
TCGA_KIRC1<-merge(TCGA_KIRC1,old_patient_stage,by="Sample_ID")
colnames(TCGA_KIRC1)[747]<-"Cluster"
#HNRNPA3,HNRNPH2,HNRNPH3,RBM15,RBM23,RBM43,RBM7,RBMS2,SRSF1
TCGA_KIRC1$Cluster <- factor(TCGA_KIRC1$Cluster,levels = c("Cluster1","Cluster2","Cluster3"))

anno_pvalue <- data.frame()
for(i in p_value$gene){
  data <- dplyr::select(TCGA_KIRC1,gene = i,Cluster = Cluster)
  out_results <- compare_means(gene~Cluster ,data = data, method ="kruskal.test")
  out_results$.y. <- i
  anno_pvalue <- rbind(anno_pvalue,out_results)
}

##重新画图
anno_pvalue<-anno_pvalue[!anno_pvalue$p.signif=="ns",]
TCGA_KIRC1<-GDSC_predict[anno_pvalue$.y.,]
TCGA_KIRC1<-as.matrix(t(scale(t(TCGA_KIRC1))))
TCGA_KIRC1<-as.matrix(t(scale(t(TCGA_KIRC1))))

# TCGA_KIRC1[TCGA_KIRC1 > 2] <- 2
# TCGA_KIRC1[TCGA_KIRC1 < -2] <- -2

TCGA_KIRC1 <- TCGA_KIRC1[, order(match(colnames(TCGA_KIRC1), old_patient_stage$Sample_ID))]
sfd = HeatmapAnnotation(Subtype = old_patient_stage$group1)

row_ha2 = rowAnnotation(foo = anno_text(anno_pvalue$p.signif  , location = 0.5, just = "center"))

# 绘制热图
library(circlize)
ComplexHeatmap::Heatmap(TCGA_KIRC1,
                        col =colorRamp2(c(-4,0,4), c("#4DBBD5FF", "#FFFFFF", "#ED0000FF")),
                        column_split  = factor(old_patient_stage$group1, levels = c("Cluster1","Cluster2","Cluster3")),
                        name = " ",
                        row_names_side = "left",
                        right_annotation = row_ha2,
                        top_annotation = sfd,
                        show_column_names = F,
                        cluster_rows = F,
                        cluster_columns = F)

####TCGA-CCLE####
rownames(CCLE_predict)<-CCLE_predict$X
CCLE_predict<-CCLE_predict[,-1]
CCLE_predict<-as.data.frame(t(CCLE_predict))
CCLE_predict <- CCLE_predict %>%
  dplyr::mutate(across(where(is.character), as.numeric))
CCLE_predict<-as.matrix(t(scale(t(CCLE_predict))))

TCGA_KIRC1<-CCLE_predict

p_value <- data.frame(gene = rownames(data.frame(TCGA_KIRC1)))

TCGA_KIRC1 <- as.data.frame(t(TCGA_KIRC1))
TCGA_KIRC1$Sample_ID <- rownames(TCGA_KIRC1)
TCGA_KIRC1<-merge(TCGA_KIRC1,old_patient_stage,by="Sample_ID")
colnames(TCGA_KIRC1)[37]<-"Cluster"
#HNRNPA3,HNRNPH2,HNRNPH3,RBM15,RBM23,RBM43,RBM7,RBMS2,SRSF1
TCGA_KIRC1$Cluster <- factor(TCGA_KIRC1$Cluster,levels = c("Cluster1","Cluster2","Cluster3"))

anno_pvalue <- data.frame()
for(i in p_value$gene){
  data <- dplyr::select(TCGA_KIRC1,gene = i,Cluster = Cluster)
  out_results <- compare_means(gene~Cluster ,data = data, method ="kruskal.test")
  out_results$.y. <- i
  anno_pvalue <- rbind(anno_pvalue,out_results)
}

##重新画图
# anno_pvalue<-anno_pvalue[!anno_pvalue$p.signif=="ns",]
TCGA_KIRC1<-CCLE_predict[anno_pvalue$.y.,]
TCGA_KIRC1<-as.matrix(t(scale(t(TCGA_KIRC1))))
TCGA_KIRC1<-as.matrix(t(scale(t(TCGA_KIRC1))))

# TCGA_KIRC1[TCGA_KIRC1 > 2] <- 2
# TCGA_KIRC1[TCGA_KIRC1 < -2] <- -2

TCGA_KIRC1 <- TCGA_KIRC1[, order(match(colnames(TCGA_KIRC1), old_patient_stage$Sample_ID))]
sfd = HeatmapAnnotation(Subtype = old_patient_stage$group1)

row_ha2 = rowAnnotation(foo = anno_text(anno_pvalue$p.signif  , location = 0.5, just = "center"))

# 绘制热图
library(circlize)
ComplexHeatmap::Heatmap(TCGA_KIRC1,
                        col =colorRamp2(c(-4,0,4), c("#4DBBD5FF", "#FFFFFF", "#ED0000FF")),
                        column_split  = factor(old_patient_stage$group1, levels = c("Cluster1","Cluster2","Cluster3")),
                        name = " ",
                        row_names_side = "left",
                        right_annotation = row_ha2,
                        top_annotation = sfd,
                        show_column_names = F,
                        cluster_rows = F,
                        cluster_columns = F)

rownames(CCLE_predict)<-CCLE_predict$X
CCLE_predict<-CCLE_predict[,-1]
CCLE_predict<-as.data.frame(t(CCLE_predict))
CCLE_predict <- CCLE_predict %>%
  dplyr::mutate(across(where(is.character), as.numeric))
CCLE_predict<-as.matrix(t(scale(t(CCLE_predict))))

TCGA_KIRC1<-CCLE_predict

p_value <- data.frame(gene = rownames(data.frame(TCGA_KIRC1)))

TCGA_KIRC1 <- as.data.frame(t(TCGA_KIRC1))
TCGA_KIRC1$Sample_ID <- rownames(TCGA_KIRC1)
TCGA_KIRC1<-merge(TCGA_KIRC1,old_patient_stage,by="Sample_ID")
colnames(TCGA_KIRC1)[37]<-"Cluster"
#HNRNPA3,HNRNPH2,HNRNPH3,RBM15,RBM23,RBM43,RBM7,RBMS2,SRSF1
TCGA_KIRC1$Cluster <- factor(TCGA_KIRC1$Cluster,levels = c("Cluster1","Cluster2","Cluster3"))

anno_pvalue <- data.frame()
for(i in p_value$gene){
  data <- dplyr::select(TCGA_KIRC1,gene = i,Cluster = Cluster)
  out_results <- compare_means(gene~Cluster ,data = data, method ="kruskal.test")
  out_results$.y. <- i
  anno_pvalue <- rbind(anno_pvalue,out_results)
}

##重新画图
anno_pvalue<-anno_pvalue[!anno_pvalue$p.signif=="ns",]
TCGA_KIRC1<-CCLE_predict[anno_pvalue$.y.,]
TCGA_KIRC1<-as.matrix(t(scale(t(TCGA_KIRC1))))
TCGA_KIRC1<-as.matrix(t(scale(t(TCGA_KIRC1))))

# TCGA_KIRC1[TCGA_KIRC1 > 2] <- 2
# TCGA_KIRC1[TCGA_KIRC1 < -2] <- -2

TCGA_KIRC1 <- TCGA_KIRC1[, order(match(colnames(TCGA_KIRC1), old_patient_stage$Sample_ID))]
sfd = HeatmapAnnotation(Subtype = old_patient_stage$group1)

row_ha2 = rowAnnotation(foo = anno_text(anno_pvalue$p.signif  , location = 0.5, just = "center"))

# 绘制热图
library(circlize)
ComplexHeatmap::Heatmap(TCGA_KIRC1,
                        col =colorRamp2(c(-4,0,4), c("#4DBBD5FF", "#FFFFFF", "#ED0000FF")),
                        column_split  = factor(old_patient_stage$group1, levels = c("Cluster1","Cluster2","Cluster3")),
                        name = " ",
                        row_names_side = "left",
                        right_annotation = row_ha2,
                        top_annotation = sfd,
                        show_column_names = F,
                        cluster_rows = F,
                        cluster_columns = F)


####TCGA-CTRP####
rownames(CTRP_predict)<-CTRP_predict$X
CTRP_predict<-CTRP_predict[,-1]
CTRP_predict<-as.data.frame(t(CTRP_predict))
CTRP_predict <- CTRP_predict %>%
  dplyr::mutate(across(where(is.character), as.numeric))
CTRP_predict<-as.matrix(t(scale(t(CTRP_predict))))

TCGA_KIRC1<-CTRP_predict

p_value <- data.frame(gene = rownames(data.frame(TCGA_KIRC1)))

TCGA_KIRC1 <- as.data.frame(t(TCGA_KIRC1))
TCGA_KIRC1$Sample_ID <- rownames(TCGA_KIRC1)
TCGA_KIRC1<-merge(TCGA_KIRC1,old_patient_stage,by="Sample_ID")
colnames(TCGA_KIRC1)[558]<-"Cluster"
#HNRNPA3,HNRNPH2,HNRNPH3,RBM15,RBM23,RBM43,RBM7,RBMS2,SRSF1
TCGA_KIRC1$Cluster <- factor(TCGA_KIRC1$Cluster,levels = c("Cluster1","Cluster2","Cluster3"))

anno_pvalue <- data.frame()
for(i in p_value$gene){
  data <- dplyr::select(TCGA_KIRC1,gene = i,Cluster = Cluster)
  out_results <- compare_means(gene~Cluster ,data = data, method ="kruskal.test")
  out_results$.y. <- i
  anno_pvalue <- rbind(anno_pvalue,out_results)
}

##重新画图
# anno_pvalue<-anno_pvalue[!anno_pvalue$p.signif=="ns",]
TCGA_KIRC1<-CTRP_predict[anno_pvalue$.y.,]
TCGA_KIRC1<-as.matrix(t(scale(t(TCGA_KIRC1))))
TCGA_KIRC1<-as.matrix(t(scale(t(TCGA_KIRC1))))

# TCGA_KIRC1[TCGA_KIRC1 > 2] <- 2
# TCGA_KIRC1[TCGA_KIRC1 < -2] <- -2

TCGA_KIRC1 <- TCGA_KIRC1[, order(match(colnames(TCGA_KIRC1), old_patient_stage$Sample_ID))]
sfd = HeatmapAnnotation(Subtype = old_patient_stage$group1)

row_ha2 = rowAnnotation(foo = anno_text(anno_pvalue$p.signif  , location = 0.5, just = "center"))

# 绘制热图
library(circlize)
ComplexHeatmap::Heatmap(TCGA_KIRC1,
                        col =colorRamp2(c(-4,0,4), c("#4DBBD5FF", "#FFFFFF", "#ED0000FF")),
                        column_split  = factor(old_patient_stage$group1, levels = c("Cluster1","Cluster2","Cluster3")),
                        name = " ",
                        row_names_side = "left",
                        right_annotation = row_ha2,
                        top_annotation = sfd,
                        show_column_names = F,
                        cluster_rows = F,
                        cluster_columns = F)

rownames(CTRP_predict)<-CTRP_predict$X
CTRP_predict<-CTRP_predict[,-1]
CTRP_predict<-as.data.frame(t(CTRP_predict))
CTRP_predict <- CTRP_predict %>%
  dplyr::mutate(across(where(is.character), as.numeric))
CTRP_predict<-as.matrix(t(scale(t(CTRP_predict))))

TCGA_KIRC1<-CTRP_predict

p_value <- data.frame(gene = rownames(data.frame(TCGA_KIRC1)))

TCGA_KIRC1 <- as.data.frame(t(TCGA_KIRC1))
TCGA_KIRC1$Sample_ID <- rownames(TCGA_KIRC1)
TCGA_KIRC1<-merge(TCGA_KIRC1,old_patient_stage,by="Sample_ID")
colnames(TCGA_KIRC1)[37]<-"Cluster"
#HNRNPA3,HNRNPH2,HNRNPH3,RBM15,RBM23,RBM43,RBM7,RBMS2,SRSF1
TCGA_KIRC1$Cluster <- factor(TCGA_KIRC1$Cluster,levels = c("Cluster1","Cluster2","Cluster3"))

anno_pvalue <- data.frame()
for(i in p_value$gene){
  data <- dplyr::select(TCGA_KIRC1,gene = i,Cluster = Cluster)
  out_results <- compare_means(gene~Cluster ,data = data, method ="kruskal.test")
  out_results$.y. <- i
  anno_pvalue <- rbind(anno_pvalue,out_results)
}

##重新画图
anno_pvalue<-anno_pvalue[!anno_pvalue$p.signif=="ns",]
TCGA_KIRC1<-CTRP_predict[anno_pvalue$.y.,]
TCGA_KIRC1<-as.matrix(t(scale(t(TCGA_KIRC1))))
TCGA_KIRC1<-as.matrix(t(scale(t(TCGA_KIRC1))))

# TCGA_KIRC1[TCGA_KIRC1 > 2] <- 2
# TCGA_KIRC1[TCGA_KIRC1 < -2] <- -2

TCGA_KIRC1 <- TCGA_KIRC1[, order(match(colnames(TCGA_KIRC1), old_patient_stage$Sample_ID))]
sfd = HeatmapAnnotation(Subtype = old_patient_stage$group1)

row_ha2 = rowAnnotation(foo = anno_text(anno_pvalue$p.signif  , location = 0.5, just = "center"))

# 绘制热图
library(circlize)
ComplexHeatmap::Heatmap(TCGA_KIRC1,
                        col =colorRamp2(c(-4,0,4), c("#4DBBD5FF", "#FFFFFF", "#ED0000FF")),
                        column_split  = factor(old_patient_stage$group1, levels = c("Cluster1","Cluster2","Cluster3")),
                        name = " ",
                        row_names_side = "left",
                        right_annotation = row_ha2,
                        top_annotation = sfd,
                        show_column_names = F,
                        cluster_rows = F,
                        cluster_columns = F)

