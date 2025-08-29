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

# From Raw_Data_Processing
load("./step1_cluster_os.Rdata")
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

####1.1、t-sne看一下四种样本的一个聚散程度（所有）####
set.seed(422)
expr <- cbind(Sample_old_cancer,Sample_old_normol ,Sample_young_cancer,Sample_young_normol)
All<-t(expr)
# save(All,file = "./outdata/All.Rdata")
# load("./outdata/All.Rdata")
All<-as.data.frame(All)
All <- All %>%
  dplyr::mutate(across(where(is.character), as.numeric))

tsne_out <- Rtsne(All,pca=FALSE,check_duplicates = FALSE,
                  perplexity = 20,theta=0.5)
#save(tsne_out,file = "./outdata/tsne_out.Rdata")
#load("/home/cmx/Pan-cancer/outdata/tsne_out.Rdata")

head(tsne_out)
plot(tsne_out$Y) 
# 使用ggplot2包可视化tSNE降维的结果

tsne_res <- as.data.frame(tsne_out$Y)
colnames(tsne_res) <- c("tSNE1","tSNE2")
head(tsne_res)
##      tSNE1    tSNE2
## 1 7.715389 17.70517
## 2 6.714282 13.76982
## 3 8.157742 14.12467
## 4 7.809143 13.69793
## 5 7.987368 17.52846
## 6 6.791518 20.15391

# 使用ggplot2可视化tSNE降维的结果
PSI_data_group<-rep(c("Sample_old_cancer","Sample_old_normol","Sample_young_cancer","Sample_young_normol"),
                    c(342,31,186,31))
PSI_data_group<-as.data.frame(PSI_data_group)
colnames(PSI_data_group)<-"group"
library(ggplot2)

df_label <- summarise(group_by(tsne_res, PSI_data_group$group), 
                      x_pos = mean(tSNE2), 
                      y_pos = mean(tSNE1))
p<-ggplot(tsne_res, aes(x = tSNE1 , y =tSNE2, color = PSI_data_group$group)) +
  geom_point(size = 1, shape = 19) +
  theme_classic()+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', 
                                                                      fill = 'transparent'))+
  theme(legend.title =element_text(size=14,face = "bold"))+
  theme(legend.position = "top")+
  labs(x = "TSNE-1", y = "TSNE-2",color = " ") +
  scale_color_manual(values = c("#191970","#9ACD32","#1E90FF", "#EE6AA7"), 
                     labels = c("old_cancer","old_normol","young_cancer","young_normol")) +
  theme(legend.text = element_text(size = 8, color = "black",face = "bold"))
p
####t-sne看一下两种样本的一个聚散程度（老年）####
expr <- cbind(Sample_old_cancer,Sample_old_normol)
All<-t(expr)
# save(All,file = "./outdata/All.Rdata")
# load("./outdata/All.Rdata")
All<-as.data.frame(All)
All <- All %>%
  dplyr::mutate(across(where(is.character), as.numeric))

tsne_out <- Rtsne(All,pca=FALSE,check_duplicates = FALSE,
                  perplexity = 20,theta=0.5)
#save(tsne_out,file = "./outdata/tsne_out.Rdata")
#load("/home/cmx/Pan-cancer/outdata/tsne_out.Rdata")

head(tsne_out)
plot(tsne_out$Y) 
# 使用ggplot2包可视化tSNE降维的结果

tsne_res <- as.data.frame(tsne_out$Y)
colnames(tsne_res) <- c("tSNE1","tSNE2")
head(tsne_res)
##      tSNE1    tSNE2
## 1 7.715389 17.70517
## 2 6.714282 13.76982
## 3 8.157742 14.12467
## 4 7.809143 13.69793
## 5 7.987368 17.52846
## 6 6.791518 20.15391

# 使用ggplot2可视化tSNE降维的结果
PSI_data_group<-rep(c("Sample_old_cancer","Sample_old_normol"),
                    c(342,31))

# 假设expr是你的表达矩阵，包含所有样本和基因
# 提取需要的两列（ARHGAP19的两个剪接事件）

PSI_data_group<-as.data.frame(PSI_data_group)
colnames(PSI_data_group)<-"group"
library(ggplot2)

df_label <- summarise(group_by(tsne_res, PSI_data_group$group), 
                      x_pos = mean(tSNE2), 
                      y_pos = mean(tSNE1))
p<-ggplot(tsne_res, aes(x = tSNE1 , y =tSNE2, color = PSI_data_group$group)) +
  geom_point(size = 1, shape = 19) +
  theme_classic()+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', 
                                                                      fill = 'transparent'))+
  theme(legend.title =element_text(size=14,face = "bold"))+
  theme(legend.position = "top")+
  labs(x = "TSNE-1", y = "TSNE-2",color = " ") +
  scale_color_manual(values = c("#191970","#9ACD32"), 
                     labels = c("old_cancer","old_normol")) +
  theme(legend.text = element_text(size = 8, color = "black",face = "bold"))
p
####删掉一个老的正常样本，一个癌症样本再来####
#Sample_old_normol$`TCGA-AC-A2FM-Norm`
#Sample_old_cancer$`TCGA-A1-A0SB`
Sample_old_normol<-Sample_old_normol[,-5]
Sample_old_cancer<-Sample_old_cancer[,-1]

####1.2、再用t-sne看一下两种样本的一个聚散程度（老年）####
expr <- cbind(Sample_old_cancer,Sample_old_normol)
All<-t(expr)
# save(All,file = "./outdata/All.Rdata")
# load("./outdata/All.Rdata")
All<-as.data.frame(All)
All <- All %>%
  dplyr::mutate(across(where(is.character), as.numeric))

tsne_out <- Rtsne(All,pca=FALSE,check_duplicates = FALSE,
                  perplexity = 20,theta=0.5)
#save(tsne_out,file = "./outdata/tsne_out.Rdata")
#load("/home/cmx/Pan-cancer/outdata/tsne_out.Rdata")

head(tsne_out)
plot(tsne_out$Y) 
# 使用ggplot2包可视化tSNE降维的结果

tsne_res <- as.data.frame(tsne_out$Y)
colnames(tsne_res) <- c("tSNE1","tSNE2")
head(tsne_res)
##      tSNE1    tSNE2
## 1 7.715389 17.70517
## 2 6.714282 13.76982
## 3 8.157742 14.12467
## 4 7.809143 13.69793
## 5 7.987368 17.52846
## 6 6.791518 20.15391

# 使用ggplot2可视化tSNE降维的结果
PSI_data_group<-rep(c("Sample_old_cancer","Sample_old_normol"),
                    c(341,30))
PSI_data_group<-as.data.frame(PSI_data_group)
colnames(PSI_data_group)<-"group"
library(ggplot2)

df_label <- summarise(group_by(tsne_res, PSI_data_group$group), 
                      x_pos = mean(tSNE2), 
                      y_pos = mean(tSNE1))
p<-ggplot(tsne_res, aes(x = tSNE1 , y =tSNE2, color = PSI_data_group$group)) +
  geom_point(size = 1, shape = 19) +
  theme_classic()+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', 
                                                                      fill = 'transparent'))+
  theme(legend.title =element_text(size=14,face = "bold"))+
  theme(legend.position = "top")+
  labs(x = "TSNE-1", y = "TSNE-2",color = " ") +
  scale_color_manual(values = c("#191970","#9ACD32"), 
                     labels = c("old_cancer","old_normol")) +
  theme(legend.text = element_text(size = 8, color = "black",face = "bold"))
p
####t-sne看一下两种样本的一个聚散程度（年轻）####
expr <- cbind(Sample_young_cancer,Sample_young_normol)
All<-t(expr)
# save(All,file = "./outdata/All.Rdata")
# load("./outdata/All.Rdata")
All<-as.data.frame(All)
All <- All %>%
  dplyr::mutate(across(where(is.character), as.numeric))

tsne_out <- Rtsne(All,pca=FALSE,check_duplicates = FALSE,
                  perplexity = 20,theta=0.5)
#save(tsne_out,file = "./outdata/tsne_out.Rdata")
#load("/home/cmx/Pan-cancer/outdata/tsne_out.Rdata")

head(tsne_out)
plot(tsne_out$Y) 
# 使用ggplot2包可视化tSNE降维的结果

tsne_res <- as.data.frame(tsne_out$Y)
colnames(tsne_res) <- c("tSNE1","tSNE2")
head(tsne_res)
##      tSNE1    tSNE2
## 1 7.715389 17.70517
## 2 6.714282 13.76982
## 3 8.157742 14.12467
## 4 7.809143 13.69793
## 5 7.987368 17.52846
## 6 6.791518 20.15391

# 使用ggplot2可视化tSNE降维的结果
PSI_data_group<-rep(c("Sample_young_cancer","Sample_young_normol"),
                    c(186,31))
PSI_data_group<-as.data.frame(PSI_data_group)
colnames(PSI_data_group)<-"group"
library(ggplot2)

df_label <- summarise(group_by(tsne_res, PSI_data_group$group), 
                      x_pos = mean(tSNE2), 
                      y_pos = mean(tSNE1))
p<-ggplot(tsne_res, aes(x = tSNE1 , y =tSNE2, color = PSI_data_group$group)) +
  geom_point(size = 1, shape = 19) +
  theme_classic()+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', 
                                                                      fill = 'transparent'))+
  theme(legend.title =element_text(size=14,face = "bold"))+
  theme(legend.position = "top")+
  labs(x = "TSNE-1", y = "TSNE-2",color = " ") +
  scale_color_manual(values = c("#1E90FF", "#EE6AA7"), 
                     labels = c("young_cancer","young_normol")) +
  theme(legend.text = element_text(size = 8, color = "black",face = "bold"))
p
####删掉年轻的一个癌症样本，一个正常再来#####
# Sample_young_normol$`TCGA-AC-A2FF-Norm`
Sample_young_cancer<-Sample_young_cancer[,-38]
Sample_young_normol<-Sample_young_normol[,-3]
Sample_young_normol<-Sample_young_normol[,-12]

####1.3、再用t-sne看一下两种样本的一个聚散程度（年轻）####
expr <- cbind(Sample_young_cancer,Sample_young_normol)
All<-t(expr)
# save(All,file = "./outdata/All.Rdata")
# load("./outdata/All.Rdata")
All<-as.data.frame(All)
All <- All %>%
  dplyr::mutate(across(where(is.character), as.numeric))

tsne_out <- Rtsne(All,pca=FALSE,check_duplicates = FALSE,
                  perplexity = 20,theta=0.5)
#save(tsne_out,file = "./outdata/tsne_out.Rdata")
#load("/home/cmx/Pan-cancer/outdata/tsne_out.Rdata")

head(tsne_out)
plot(tsne_out$Y) 
# 使用ggplot2包可视化tSNE降维的结果

tsne_res <- as.data.frame(tsne_out$Y)
colnames(tsne_res) <- c("tSNE1","tSNE2")
head(tsne_res)
##      tSNE1    tSNE2
## 1 7.715389 17.70517
## 2 6.714282 13.76982
## 3 8.157742 14.12467
## 4 7.809143 13.69793
## 5 7.987368 17.52846
## 6 6.791518 20.15391

# 使用ggplot2可视化tSNE降维的结果
PSI_data_group<-rep(c("Sample_young_cancer","Sample_young_normol"),
                    c(185,29))
PSI_data_group<-as.data.frame(PSI_data_group)
colnames(PSI_data_group)<-"group"
library(ggplot2)

df_label <- summarise(group_by(tsne_res, PSI_data_group$group), 
                      x_pos = mean(tSNE2), 
                      y_pos = mean(tSNE1))
p<-ggplot(tsne_res, aes(x = tSNE1 , y =tSNE2, color = PSI_data_group$group)) +
  geom_point(size = 1, shape = 19) +
  theme_classic()+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', 
                                                                      fill = 'transparent'))+
  theme(legend.title =element_text(size=14,face = "bold"))+
  theme(legend.position = "top")+
  labs(x = "TSNE-1", y = "TSNE-2",color = " ") +
  scale_color_manual(values = c("#1E90FF", "#EE6AA7"), 
                     labels = c("young_cancer","young_normol")) +
  theme(legend.text = element_text(size = 8, color = "black",face = "bold"))
p

####年轻正常vs年轻癌症####
young_group<-cbind(Sample_young_normol,Sample_young_cancer) 


inputFile=young_group
conNum=185#cancer组织样品数目

treaNum=29#normal组织样品数目

#读取输入文件

outTab=data.frame()

grade=c(rep(1,treaNum),rep(2,conNum))

#将肿瘤组织变为1，正常组织为2

rt=as.matrix(young_group)

exp=rt[,1:ncol(rt)]#去除第一列基因名

head(exp)[,1:10]

head(rt)[,1:10]

#这里是要制作表达基因矩阵

dimnames=list(rownames(exp),colnames(exp))

data=matrix(as.numeric(as.matrix(exp)),nrow = nrow(exp),dimnames = dimnames)


#差异分析,Wilcoxon检验方法

for (i in row.names(data)) {
  
  geneName=unlist(strsplit(i,"\\|",))[1]
  
  geneName=gsub("\\/","_",geneName)
  
  rt=rbind(expression=data[i,],grade=grade)
  
  rt=as.matrix(t(rt))
  
  wilcoxTest<-t.test(expression~grade,data=rt)
  
  conGeneMean=mean(data[i,1:conNum])
  
  treatGeneMean=mean(data[i,(conNum+1):ncol(data)])
  
  logFC=log2(treatGeneMean)-log2(conGeneMean)
  FC=2^(log2(treatGeneMean)-log2(conGeneMean))
  
  pValue=wilcoxTest$p.value
  
  conMed=median(data[i,1:conNum])
  
  treatMed=median(data[i,(conNum+1):ncol(data)])
  
  diffMed=treatMed-conMed
  
  if(((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){
    
    outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMean,treaMean=treatGeneMean,logFC=logFC,pValue=pValue
                              
    ))  }
  
}

pValue=outTab[,"pValue"]

fdr=p.adjust(as.numeric(as.vector(pValue)),method = "BH")

outTab=cbind(outTab,fdr=fdr)

outTab$logFC<-as.numeric(outTab$logFC)

#outTab1.2 <- outTab[abs(outTab[,4]) >= 1.2, ]
outTab_young <- outTab[(outTab$fdr < 0.05 & (outTab$logFC>log2(1.2) | outTab$logFC < (-log2(1.2)))),]
outTab_young$logFC<- -outTab_young$logFC
# 提取最后两个字母并创建新列"type" 
outTab_young$type <- substr(outTab_young$gene, nchar(outTab_young$gene) - 1, nchar(outTab_young$gene))
table(outTab_young$type)
outTab_young$up_down<-ifelse(outTab_young$logFC >0,"Up",
                             ifelse(outTab_young$logFC < 0,"Down","Unchanged"))
####老年正常vs老年癌症####
old_group<-cbind(Sample_old_normol,Sample_old_cancer) 

inputFile=old_group
conNum=341#cancer组织样品数目

treaNum=30#normal组织样品数目

#读取输入文件

outTab=data.frame()

grade=c(rep(1,treaNum),rep(2,conNum))

#将肿瘤组织变为1，正常组织为2

rt=as.matrix(old_group)

exp=rt[,1:ncol(rt)]#去除第一列基因名

head(exp)[,1:10]

head(rt)[,1:10]

#这里是要制作表达基因矩阵

dimnames=list(rownames(exp),colnames(exp))

data=matrix(as.numeric(as.matrix(exp)),nrow = nrow(exp),dimnames = dimnames)


#差异分析,Wilcoxon检验方法

for (i in row.names(data)) {
  
  geneName=unlist(strsplit(i,"\\|",))[1]
  
  geneName=gsub("\\/","_",geneName)
  
  rt=rbind(expression=data[i,],grade=grade)
  
  rt=as.matrix(t(rt))
  
  wilcoxTest<-t.test(expression~grade,data=rt)
  
  conGeneMean=mean(data[i,1:conNum])
  
  treatGeneMean=mean(data[i,(conNum+1):ncol(data)])
  
  logFC=log2(treatGeneMean)-log2(conGeneMean)
  FC=2^(log2(treatGeneMean)-log2(conGeneMean))
  
  pValue=wilcoxTest$p.value
  
  conMed=median(data[i,1:conNum])
  
  treatMed=median(data[i,(conNum+1):ncol(data)])
  
  diffMed=treatMed-conMed
  
  if(((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){
    
    outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMean,treaMean=treatGeneMean,logFC=logFC,pValue=pValue
                              
    ))  }
  
}

pValue=outTab[,"pValue"]

fdr=p.adjust(as.numeric(as.vector(pValue)),method = "BH")

outTab=cbind(outTab,fdr=fdr)

outTab$logFC<-as.numeric(outTab$logFC)

#outTab1.2 <- outTab[abs(outTab[,4]) >= 1.2, ]
outTab_old <- outTab[(outTab$fdr < 0.05 & (outTab$logFC>log2(1.2) | outTab$logFC < (-log2(1.2)))),]
outTab_old$logFC<- -outTab_old$logFC
# 提取最后两个字母并创建新列"type" 
outTab_old$type <- substr(outTab_old$gene, nchar(outTab_old$gene) - 1, nchar(outTab_old$gene))
table(outTab_old$type)
outTab_old$up_down<-ifelse(outTab_old$logFC > 0,"Up",
                           ifelse(outTab_old$logFC < 0,"Down","Unchanged"))
Sample_old_cancer1<-Sample_old_cancer[outTab_old$gene,]
Sample_old_normol1<-Sample_old_normol[outTab_old$gene,]
ALL_EXP<-cbind(Sample_old_cancer1,Sample_old_normol1)
library(pheatmap)
library(RColorBrewer)
# 列名含-Norm的样本标记为Normal，其余为Tumor
sample_types <- ifelse(grepl("-Norm$", colnames(ALL_EXP)), "Normal", "Tumor")

# 验证分组数量
table(sample_types)
# 输出示例：
# Tumor  Normal 
#   330     20 
ALL_EXP <- All %>%
  dplyr::mutate(across(where(is.character), as.numeric))
ALL_EXP_scaled <- as.data.frame(scale(ALL_EXP))  # 按行（基因）标准化
# 假设ALL_EXP是你的数据框，行名是基因，列名是样本
# 此处直接使用用户提供的ALL_EXP数据（已修正列名中的-Norm）
# 创建注释数据框
annotation <- data.frame(
  Sample_Type = factor(sample_types, levels = c("Normal", "Tumor"))
)
rownames(annotation) <- colnames(ALL_EXP_scaled)

# 颜色方案
color_annotation <- list(
  Sample_Type = c(Normal = "#4477AA", Tumor = "#CC6677")  # 蓝=正常，红=肿瘤
)

# 热图颜色渐变（蓝低表达，红高表达）
color_heatmap <- colorRampPalette(brewer.pal(9, "RdBu"))(255)
pheatmap::pheatmap(
  ALL_EXP_scaled,
  scale="row",
  # 聚类参数
  clustering_row = TRUE,          # 行聚类（基因）
  clustering_col = TRUE,          # 列聚类（样本）
  # method = "ward.D2",             # 聚类算法
  distance = "euclidean",         # 距离度量
  # 注释参数
  annotation_col = annotation,
  annotation_colors = color_annotation,
  annotation_names_col = FALSE,   # 隐藏注释标题
  # 可视化参数
  show_rownames = F,           # 显示基因名（若基因数少）
  show_colnames = FALSE,          # 隐藏样本ID
  color = color_heatmap,
  na_col = "gray",                # 缺失值颜色
  # main = "Clustered Heatmap of AS Events\n(Normal Samples Marked with -Norm)",
  # filename = "Figure_S2C_Clustered_Heatmap.png",
  # width = 10, height = 8
)

save(outTab_young,outTab_old,
     file = "./outdata/outTab_YwO.Rdata")
####1.4画图表示一下他们的差异剪接事件类型####
AllType <-
  data.frame(
    'Type' = c(
      "old",
      "young"),
    'AA' = c(57,30),
    'AD' = c(89,51),
    'AP' = c(699,475),
    'AT' = c(338,232),
    'ES' = c(499,392),
    'ME' = c(10,9),
    'RI' = c(133,61)
  )

measure_name=setdiff(colnames(AllType),
                     c('Type'))
library(data.table)
#id.vars:合并的过程中通过那一列进行合并
#measure.vars：将那些列合并在一起
#variable.name：将那些列合并在一起后，新列的名称
#value.name：矩阵内部的数值，赋予一个新的名称
data1=reshape2::melt(AllType,
                     id.vars='Type', 
                     measure.vars=measure_name,
                     variable.name = "AS_Type", 
                     value.name = "numbers")
head(data1)

#Group转换为因子，并排序
group = factor(data1$Type,levels=unique(data1$Type),order=TRUE)
#设置颜色
mycolors<-c("#EEB422","#6495ED","#3A5FCD","#FFF68F","#8FBC8F","#228B22","#A52A2A")

#作图
p <- ggplot(data=data1,aes(x=Type,y=numbers,fill=AS_Type)) + 
  geom_bar(stat="identity",position="fill") + 
  scale_fill_manual(values=mycolors)+
  scale_y_continuous(expand = expansion(mult=c(0.01,0.02)),
                     labels = scales::percent_format())+
  labs(x=" ",y="Proportion of AS events",
       fill=" ",title="")+
  theme_classic()+
  theme(axis.title.y=element_text(size=14))+
  theme(legend.text=element_text(size=10))+
  theme(axis.text.x = element_text(size = 12, color = "black"))+
  theme(axis.text.y = element_text(size = 12, color = "black"))+
  theme(axis.ticks.length=unit(0.3,"cm"))+
  theme(axis.text.x=element_text(angle=0,vjust=1,hjust=1,size=11))
p
####1.5 画图表示一下他们的差异事件数量####
# 示例数据
dat <- data.frame(Type = c("old", "young"), Total = c(1825, 1250))

# 绘制柱状图
custom_colors <- c("#3A5FCD","#8FBC8F")
# 绘制柱状图
ggplot(dat, aes(x = Type, y = Total, fill = Type)) + 
  geom_col(color = "white") +
  geom_text(aes(label = Total), position = position_stack(vjust = 1.05), size = 4, color = "black") +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(x = " ", y = "Total Number of AS events", fill = "", title = "") +
  theme_classic() +
  theme(axis.title.y = element_text(size = 14)) +
  theme(legend.text = element_text(size = 10)) +
  theme(axis.text.x = element_text(size = 12, color = "black")) +
  theme(axis.text.y = element_text(size = 12, color = "black")) +
  theme(axis.ticks.length = unit(0.3, "cm")) +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1, size = 11))
dat <-
  data.frame(
    'Type' = c(
      "old",
      "young"),
    'Up' = c(643,607),
    'Down' = c(1123,702)
  )

table <- gather(dat, key = Change, value = Count,-Type) #数据转换，宽长数据
table$Type <- factor(table$Type)
table[which(table$Change == 'Down'), c('Count')] <-
  table[which(table$Change == 'Down'), c('Count')] * -1   #将另一个样本的数据转化为负数，这是必须的一步
p <-
  ggplot(table, aes(Type, Count, fill = Change)) + geom_col(width = 0.85) +
  geom_bar(
    stat = 'identity',
    fill = ifelse(table$Change == "Up", '#A52A2A', '#458B74'),
    # 根据y值的正负设置颜色
    width = 0.85
  ) +
  labs(x = " ", y = "Number of differential AS") +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(color = 'black',
                                    fill = 'transparent'),
    legend.title = element_blank(),
    axis.text.x  = element_text(
      angle = 0,
      hjust = 1,
      color = "black",
      size = 12
    ),
    axis.text.y  = element_text(
      angle =0,
      hjust = 1,
      color = "black",
      size = 12
    ),
    axis.title = element_text(color = "black", size = 15)
  )  + #添加中间的线，或者不添加也行
  annotate('text', label = 'Up in cancer',2, 1000) + annotate('text', label = 'Down in cancer',2, -1000)  #添加注释信息
p
#### t-sne看一下四种样本的一个聚散程度（交集差异事件）####
#只要他们扰动事件的
pase_Sample_old_cancer<- Sample_old_cancer[rownames(Sample_old_cancer) %in% outTab_old$gene, ]
pase_Sample_old_normol<- Sample_old_normol[rownames(Sample_old_normol) %in% outTab_old$gene, ]

pase_Sample_young_cancer<- Sample_young_cancer[rownames(Sample_young_cancer) %in% outTab_young$gene, ]
pase_Sample_young_normol<- Sample_young_normol[rownames(Sample_young_normol) %in% outTab_young$gene, ]

###intersect函数取交集
library(ComplexUpset)
library(tidyverse)
as_names <- intersect(rownames(pase_Sample_old_cancer),rownames(pase_Sample_young_cancer))
# 提取基因名和事件类型
as_df <- tibble(
  full_name = as_names,
  gene = str_extract(as_names, "^[^-]+"),
  type = str_extract(as_names, "[A-Z]+$")
)

# 创建 wide format：每行一个基因，列为事件类型，值为是否有该事件
as_wide <- as_df %>%
  distinct(gene, type) %>%
  pivot_wider(names_from = type, values_from = type,
              values_fn = length, values_fill = 0) %>%
  mutate(across(-gene, ~ . > 0))

# UpSet 图
ComplexUpset::upset(as_wide, intersect = colnames(as_wide)[-1],
                    name = 'AS Types',
                    width_ratio = 0.2)

pase_Sample_old_cancer<- pase_Sample_old_cancer[rownames(pase_Sample_old_cancer) %in% as_names, ]
pase_Sample_old_normol<- pase_Sample_old_normol[rownames(pase_Sample_old_normol) %in% as_names, ]

pase_Sample_young_cancer<- pase_Sample_young_cancer[rownames(pase_Sample_young_cancer) %in% as_names, ]
pase_Sample_young_normol<- pase_Sample_young_normol[rownames(pase_Sample_young_normol) %in% as_names, ]

##cbind进行合并
#BLCA,BRCA,CHOL,COAD,ESCA,GBM,HNSC,KICH,KIRC,KIRP,LIHC,LUAD,LUSC,PRAD,READ,STAD,THCA,UCEC 

expr <- cbind(pase_Sample_old_cancer,pase_Sample_old_normol ,pase_Sample_young_cancer,pase_Sample_young_normol)
All<-t(expr)
# save(All,file = "./outdata/All.Rdata")
# load("./outdata/All.Rdata")
All<-as.data.frame(All)
All <- All %>%
  dplyr::mutate(across(where(is.character), as.numeric))

tsne_out <- Rtsne(All,pca=FALSE,check_duplicates = FALSE,
                  perplexity = 20,theta=0.5)
#save(tsne_out,file = "./outdata/tsne_out.Rdata")
#load("/home/cmx/Pan-cancer/outdata/tsne_out.Rdata")

head(tsne_out)
plot(tsne_out$Y) 
# 使用ggplot2包可视化tSNE降维的结果

tsne_res <- as.data.frame(tsne_out$Y)
colnames(tsne_res) <- c("tSNE1","tSNE2")
head(tsne_res)
##      tSNE1    tSNE2
## 1 7.715389 17.70517
## 2 6.714282 13.76982
## 3 8.157742 14.12467
## 4 7.809143 13.69793
## 5 7.987368 17.52846
## 6 6.791518 20.15391

# 使用ggplot2可视化tSNE降维的结果
PSI_data_group<-rep(c("Sample_old_cancer","Sample_old_normol","Sample_young_cancer","Sample_young_normol"),
                    c(341,30,185,29))
PSI_data_group<-as.data.frame(PSI_data_group)
colnames(PSI_data_group)<-"group"
library(ggplot2)

df_label <- summarise(group_by(tsne_res, PSI_data_group$group), 
                      x_pos = mean(tSNE2), 
                      y_pos = mean(tSNE1))
p<-ggplot(tsne_res, aes(x = tSNE1 , y =tSNE2, color = PSI_data_group$group)) +
  geom_point(size = 1, shape = 19) +
  theme_classic()+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', 
                                                                      fill = 'transparent'))+
  theme(legend.title =element_text(size=14,face = "bold"))+
  theme(legend.position = "top")+
  labs(x = "TSNE-1", y = "TSNE-2",color = " ") +
  scale_color_manual(values = c("#191970","#9ACD32","#1E90FF", "#EE6AA7"), 
                     labels = c("old_cancer","old_normol","young_cancer","young_normol")) +
  theme(legend.text = element_text(size = 8, color = "black",face = "bold"))
p
####1.6t-sne看一下四种样本的一个聚散程度（并集差异事件）####
###intersect函数取交集
all_names <- union(outTab_old$gene,outTab_young$gene)

pase_Sample_old_cancer<- Sample_old_cancer[rownames(Sample_old_cancer) %in% all_names, ]
pase_Sample_old_normol<- Sample_old_normol[rownames(Sample_old_normol) %in% all_names, ]

pase_Sample_young_cancer<- Sample_young_cancer[rownames(Sample_young_cancer) %in% all_names, ]
pase_Sample_young_normol<- Sample_young_normol[rownames(Sample_young_normol) %in% all_names, ]

##cbind进行合并
#BLCA,BRCA,CHOL,COAD,ESCA,GBM,HNSC,KICH,KIRC,KIRP,LIHC,LUAD,LUSC,PRAD,READ,STAD,THCA,UCEC 

expr <- cbind(pase_Sample_old_cancer,pase_Sample_old_normol ,pase_Sample_young_cancer,pase_Sample_young_normol)
All<-t(expr)
# save(All,file = "./outdata/All.Rdata")
# load("./outdata/All.Rdata")
All<-as.data.frame(All)
All <- All %>%
  dplyr::mutate(across(where(is.character), as.numeric))

tsne_out <- Rtsne(All,pca=FALSE,check_duplicates = FALSE,
                  perplexity = 20,theta=0.5)
#save(tsne_out,file = "./outdata/tsne_out.Rdata")
#load("/home/cmx/Pan-cancer/outdata/tsne_out.Rdata")

head(tsne_out)
plot(tsne_out$Y) 
# 使用ggplot2包可视化tSNE降维的结果

tsne_res <- as.data.frame(tsne_out$Y)
colnames(tsne_res) <- c("tSNE1","tSNE2")
head(tsne_res)
##      tSNE1    tSNE2
## 1 7.715389 17.70517
## 2 6.714282 13.76982
## 3 8.157742 14.12467
## 4 7.809143 13.69793
## 5 7.987368 17.52846
## 6 6.791518 20.15391

# 使用ggplot2可视化tSNE降维的结果
PSI_data_group<-rep(c("Sample_old_cancer","Sample_old_normol","Sample_young_cancer","Sample_young_normol"),
                    c(341,30,185,29))
PSI_data_group<-as.data.frame(PSI_data_group)
colnames(PSI_data_group)<-"group"
library(ggplot2)

df_label <- summarise(group_by(tsne_res, PSI_data_group$group), 
                      x_pos = mean(tSNE2), 
                      y_pos = mean(tSNE1))
p<-ggplot(tsne_res, aes(x = tSNE1 , y =tSNE2, color = PSI_data_group$group)) +
  geom_point(size = 1, shape = 19) +
  theme_classic()+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', 
                                                                      fill = 'transparent'))+
  theme(legend.title =element_text(size=14,face = "bold"))+
  theme(legend.position = "top")+
  labs(x = "TSNE-1", y = "TSNE-2",color = " ") +
  scale_color_manual(values = c("#191970","#9ACD32","#1E90FF", "#EE6AA7"), 
                     labels = c("old_cancer","old_normol","young_cancer","young_normol")) +
  theme(legend.text = element_text(size = 8, color = "black",face = "bold"))
p
#### t-sne看一下老年这些差异差异剪接事件确实能区分癌症和正常####
###intersect函数取交集
all_names <- outTab_old$gene

pase_Sample_old_cancer<- Sample_old_cancer[rownames(Sample_old_cancer) %in% all_names, ]
pase_Sample_old_normol<- Sample_old_normol[rownames(Sample_old_normol) %in% all_names, ]

##cbind进行合并
#BLCA,BRCA,CHOL,COAD,ESCA,GBM,HNSC,KICH,KIRC,KIRP,LIHC,LUAD,LUSC,PRAD,READ,STAD,THCA,UCEC 

expr <- cbind(pase_Sample_old_cancer,pase_Sample_old_normol)
All<-t(expr)
# save(All,file = "./outdata/All.Rdata")
# load("./outdata/All.Rdata")
All<-as.data.frame(All)
All <- All %>%
  dplyr::mutate(across(where(is.character), as.numeric))

tsne_out <- Rtsne(All,pca=FALSE,check_duplicates = FALSE,
                  perplexity = 20,theta=0.5)
#save(tsne_out,file = "./outdata/tsne_out.Rdata")
#load("/home/cmx/Pan-cancer/outdata/tsne_out.Rdata")

head(tsne_out)
plot(tsne_out$Y) 
# 使用ggplot2包可视化tSNE降维的结果

tsne_res <- as.data.frame(tsne_out$Y)
colnames(tsne_res) <- c("tSNE1","tSNE2")
head(tsne_res)
##      tSNE1    tSNE2
## 1 7.715389 17.70517
## 2 6.714282 13.76982
## 3 8.157742 14.12467
## 4 7.809143 13.69793
## 5 7.987368 17.52846
## 6 6.791518 20.15391

# 使用ggplot2可视化tSNE降维的结果
PSI_data_group<-rep(c("Sample_old_cancer","Sample_old_normol"),
                    c(341,30))
PSI_data_group<-as.data.frame(PSI_data_group)
colnames(PSI_data_group)<-"group"
library(ggplot2)

df_label <- summarise(group_by(tsne_res, PSI_data_group$group), 
                      x_pos = mean(tSNE2), 
                      y_pos = mean(tSNE1))
p<-ggplot(tsne_res, aes(x = tSNE1 , y =tSNE2, color = PSI_data_group$group)) +
  geom_point(size = 1, shape = 19) +
  theme_classic()+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', 
                                                                      fill = 'transparent'))+
  theme(legend.title =element_text(size=14,face = "bold"))+
  theme(legend.position = "top")+
  labs(x = "TSNE-1", y = "TSNE-2",color = " ") +
  scale_color_manual(values = c("#191970","#9ACD32"), 
                     labels = c("old_cancer","old_normol")) +
  theme(legend.text = element_text(size = 8, color = "black",face = "bold"))
p

####1.7、t-sne看一下老人特异剪接的一个聚散程度####
load("/home/cmx/senes_splicing/indata/390as.Rdata")
###intersect函数取交集
all_names <-rownames(PSI_data_group)
pase_Sample_old_cancer<- Sample_old_cancer[rownames(Sample_old_cancer) %in% all_names, ]
pase_Sample_old_normol<- Sample_old_normol[rownames(Sample_old_normol) %in% all_names, ]

pase_Sample_young_cancer<- Sample_young_cancer[rownames(Sample_young_cancer) %in% all_names, ]
pase_Sample_young_normol<- Sample_young_normol[rownames(Sample_young_normol) %in% all_names, ]
##cbind进行合并
#BLCA,BRCA,CHOL,COAD,ESCA,GBM,HNSC,KICH,KIRC,KIRP,LIHC,LUAD,LUSC,PRAD,READ,STAD,THCA,UCEC 

expr <- cbind(pase_Sample_old_cancer,pase_Sample_old_normol)
All<-t(expr)
# load("./outdata/All.Rdata")
# save(All,file = "./outdata/All.Rdata")
####下面得到的390个高变特异剪接事件####
All<-as.data.frame(All)
All <- All %>%
  dplyr::mutate(across(where(is.character), as.numeric))
# PSI_data_top$event<-rownames(PSI_data_top)
# All<-All[,PSI_data_top$event]
tsne_out <- Rtsne(All,pca=FALSE,check_duplicates = FALSE,
                  perplexity = 50,theta=0.5)
#save(tsne_out,file = "./outdata/tsne_out.Rdata")
#load("/home/cmx/Pan-cancer/outdata/tsne_out.Rdata")

head(tsne_out)
plot(tsne_out$Y) 
# 使用ggplot2包可视化tSNE降维的结果

tsne_res <- as.data.frame(tsne_out$Y)
colnames(tsne_res) <- c("tSNE1","tSNE2")
head(tsne_res)
##      tSNE1    tSNE2
## 1 7.715389 17.70517
## 2 6.714282 13.76982
## 3 8.157742 14.12467
## 4 7.809143 13.69793
## 5 7.987368 17.52846
## 6 6.791518 20.15391

# 使用ggplot2可视化tSNE降维的结果
PSI_data_group<-rep(c("Sample_old_cancer","Sample_old_normol"),
                    c(341,30))
PSI_data_group<-as.data.frame(PSI_data_group)
colnames(PSI_data_group)<-"group"
library(ggplot2)

df_label <- summarise(group_by(tsne_res, PSI_data_group$group), 
                      x_pos = mean(tSNE2), 
                      y_pos = mean(tSNE1))
p<-ggplot(tsne_res, aes(x = tSNE1 , y =tSNE2, color = PSI_data_group$group)) +
  geom_point(size = 1, shape = 19) +
  theme_classic()+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', 
                                                                      fill = 'transparent'))+
  theme(legend.title =element_text(size=14,face = "bold"))+
  theme(legend.position = "top")+
  labs(x = "TSNE-1", y = "TSNE-2",color = " ") +
  scale_color_manual(values = c("#191970","#9ACD32"), 
                     labels = c("old_cancer","old_normol")) +
  theme(legend.text = element_text(size = 8, color = "black",face = "bold"))
p
####看看老年特异的剪接事件是不是在年轻人里也看看是不是不发挥作用####
expr <- cbind(pase_Sample_young_cancer,pase_Sample_young_normol)
All<-t(expr)
# save(All,file = "./outdata/All.Rdata")
# load("./outdata/All.Rdata")
All<-as.data.frame(All)
All <- All %>%
  dplyr::mutate(across(where(is.character), as.numeric))
All<-All[,1:20]
# All<-All[,PSI_data_top$event]
# All<-All[,1:25]
tsne_out <- Rtsne(All,pca=FALSE,check_duplicates = FALSE,
                  perplexity = 50,theta=0.5)
# save(tsne_out,file = "./outdata/tsne_out.Rdata")
# load("./outdata/tsne_out.Rdata")

head(tsne_out)
plot(tsne_out$Y) 
# 使用ggplot2包可视化tSNE降维的结果

tsne_res <- as.data.frame(tsne_out$Y)
colnames(tsne_res) <- c("tSNE1","tSNE2")
head(tsne_res)
##      tSNE1    tSNE2
## 1 7.715389 17.70517
## 2 6.714282 13.76982
## 3 8.157742 14.12467
## 4 7.809143 13.69793
## 5 7.987368 17.52846
## 6 6.791518 20.15391

# 使用ggplot2可视化tSNE降维的结果
PSI_data_group<-rep(c("Sample_young_cancer","Sample_young_normol"),
                    c(185,29))
PSI_data_group<-as.data.frame(PSI_data_group)
colnames(PSI_data_group)<-"group"
library(ggplot2)

df_label <- summarise(group_by(tsne_res, PSI_data_group$group), 
                      x_pos = mean(tSNE2), 
                      y_pos = mean(tSNE1))
p<-ggplot(tsne_res, aes(x = tSNE1 , y =tSNE2, color = PSI_data_group$group)) +
  geom_point(size = 1, shape = 19) +
  theme_classic()+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', 
                                                                      fill = 'transparent'))+
  theme(legend.title =element_text(size=14,face = "bold"))+
  theme(legend.position = "top")+
  labs(x = "TSNE-1", y = "TSNE-2",color = " ") +
  scale_color_manual(values = c("#1E90FF", "#EE6AA7"), 
                     labels = c("young_cancer","young_normol")) +
  theme(legend.text = element_text(size = 8, color = "black",face = "bold"))
p

####按照方差最大的前1000个AS事件对这些患者进行NMF聚类,循环####
genes_only_in_old <- setdiff(outTab_old$gene, outTab_young$gene)
Sample_old_cancer2<-Sample_old_cancer[genes_only_in_old,]
Sample_old_normol2<-Sample_old_normol1[genes_only_in_old,]

ALL_EXP <- cbind(Sample_old_cancer2,Sample_old_normol2)
library(RColorBrewer)
# 列名含-Norm的样本标记为Normal，其余为Tumor
sample_types <- ifelse(grepl("-Norm$", colnames(ALL_EXP)), "Normal", "Tumor")

# 验证分组数量
table(sample_types)
# 输出示例：
# Tumor  Normal 
#   330     20 
ALL_EXP <- ALL_EXP %>%
  dplyr::mutate(across(where(is.character), as.numeric))
ALL_EXP_scaled <- as.data.frame(scale(ALL_EXP))  # 按行（基因）标准化
# 假设ALL_EXP是你的数据框，行名是基因，列名是样本
# 此处直接使用用户提供的ALL_EXP数据（已修正列名中的-Norm）
# 创建注释数据框
annotation <- data.frame(
  Sample_Type = factor(sample_types, levels = c("Normal", "Tumor"))
)
rownames(annotation) <- colnames(ALL_EXP_scaled)

ComplexHeatmap::Heatmap(as.matrix(ALL_EXP),
                        col = colorRampPalette(c("navy","white","firebrick3"))(100),
                        column_split  = factor(annotation$Sample_Type, levels = c("Tumor","Normal")),
                        name = " ",
                        row_names_side = "left",
                        # right_annotation = row_ha2,
                        top_annotation = annotation,
                        show_column_names = F,
                        cluster_rows = F,
                        cluster_columns = F)


# 如果还没处理数值类型，可直接强转（确保ALL_EXP是纯数值矩阵）
ALL_EXP <- apply(ALL_EXP, 2, as.numeric)
rownames(ALL_EXP) <- genes_only_in_old  # 确保行名是基因

# 确保列名就是样本名
colnames(ALL_EXP) <- colnames(cbind(Sample_old_cancer2, Sample_old_normol2))

# 分组：标记Normal与Tumor样本
sample_types <- ifelse(grepl("-Norm$", colnames(ALL_EXP)), "Normal", "Tumor")

# 创建注释对象（使用 ComplexHeatmap::HeatmapAnnotation）
annotation_ha <- HeatmapAnnotation(
  Sample_Type = factor(sample_types, levels = c("Tumor", "Normal")),
  col = list(Sample_Type = c("Normal" = "#4477AA", "Tumor" = "#CC6677")),
  annotation_name_side = "left"
)

# 画热图（无聚类，按分组显示）
ALL_EXP_scaled <- as.data.frame(t(scale(t(ALL_EXP)))) # 按行（基因）标准化

Heatmap(as.matrix(ALL_EXP_scaled),
        
        # name = "Expression",
        # col = colorRamp2(c(min(ALL_EXP), 0, max(ALL_EXP)), 
        #                  c("navy", "white", "firebrick3")),
        top_annotation = annotation_ha,
        column_split = factor(sample_types, levels = c("Tumor", "Normal")),
        show_column_names = FALSE,
        show_row_names = FALSE,
        cluster_rows = T,
        cluster_columns = FALSE,
        row_names_side = "left")
