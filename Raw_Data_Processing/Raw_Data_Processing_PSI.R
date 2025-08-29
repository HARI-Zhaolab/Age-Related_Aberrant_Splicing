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
# PSI_data <- PSI_data[!(PSI_data$from_exon == "null" & PSI_data$to_exon == "null"),]
# # 检查每一列是否可以成功转换为数值
# check_conversion <- function(x) {
#   suppressWarnings(as.numeric(x)) %>%
#     is.na() %>%
#     sum() == 0
# }
# 
# # 仅转换可以成功转换为数值的列
# PSI_data[, 11:1217] <- PSI_data[, 11:1217] %>%
#   mutate(across(where(is.character) & where(check_conversion), as.numeric))
# 
# # 对于无法成功转换的列，可以单独处理，比如用NA替换
# PSI_data[, 11:1217] <- PSI_data[, 11:1217] %>%
#   mutate(across(where(is.character) & !where(check_conversion), ~as.numeric(gsub("[^0-9.-]", "", .))))
# 
# # 检查转换结果
# head(PSI_data[, 11:1217])

# 使用K近邻（KNN）算法填充PDUI或PSI值的缺失值
#PSI_data_filled <- knnImputation(PSI_data[11:1217], k=10)#跑太久了，不想再跑了
#write.csv(PSI_data_filled,"./outdata/PSI_data_filled.csv")

PSI_data_filled<-read.csv("./outdata/PSI_data_filled.csv",header = TRUE,sep = ",")

PSI_data_filled<-PSI_data_filled[,-1]
colnames(PSI_data_filled) <- gsub("\\.","-", colnames(PSI_data_filled))

PSI_data_data<-PSI_data[,1:10]
PSI_data<-cbind(PSI_data_data,PSI_data_filled)

# 创建新的行名，并删除空格
Splice_event <- apply(PSI_data[, 1:3], 1, function(x) gsub("\\s", "", paste(x, collapse = "-")))

# 将新的行名赋值给数据框的行名
rownames(PSI_data) <- Splice_event

PSI_data$as<-rownames(PSI_data)
# PSI_data_top$as<-rownames(PSI_data_top)
# dssd<-merge(PSI_data,PSI_data_top,by="as")

PSI_data<-PSI_data[,c(11:1217)]
PSI_data<-as.data.frame(t(PSI_data))
PSI_data$Sample_ID<-row.names(PSI_data)

PSI_data <- PSI_data %>%
  dplyr::select(Sample_ID, everything())

colnames(PSI_data) <- gsub(" ", "", colnames(PSI_data))

col_names <- rownames(PSI_data)

# 按照包含"-Norm"的列排序
new_order <- c(
  col_names[grepl("-Norm", col_names)],
  col_names[!grepl("-Norm", col_names)]
)

# 使用order函数重新排列列
PSI_data <- PSI_data[new_order, ]
PSI_data_normal<-PSI_data[1:113,]

# 删除"-Norm"部分
PSI_data_normal$Sample_ID <- sub("-Norm", "", PSI_data_normal$Sample_ID)

PSI_data_cancer<-PSI_data[114:1207,]
