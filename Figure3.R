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
#看一下这三个表达模式的T-sne#
PSI_data_group<-PSI_data_top[,clinical_OS$sample]
PSI_data_group <- PSI_data_group[,order(match(colnames(PSI_data_group), clinical_OS$sample))]
# 更改perplexity和theta值
tsne_out <- Rtsne(as.data.frame(t(PSI_data_group)),pca=FALSE)#,perplexity = 20,theta=0.5

head(tsne_out)

# 使用ggplot2包可视化tSNE降维的结果
tsne_res <- as.data.frame(tsne_out$Y)
colnames(tsne_res) <- c("tSNE1","tSNE2")
rownames(tsne_res) <- clinical_OS$sample
head(tsne_res)
##      tSNE1    tSNE2
## 1 7.715389 17.70517
## 2 6.714282 13.76982
## 3 8.157742 14.12467
## 4 7.809143 13.69793
## 5 7.987368 17.52846
## 6 6.791518 20.15391

#### 1.使用ggplot2可视化tSNE降维的结果####
ggplot(tsne_res, aes(x = tSNE2, y = tSNE1, colour = clinical_OS$group1)) +
  geom_point(size = 2, shape = 16) +
  stat_ellipse(aes(fill=clinical_OS$group1,colour=clinical_OS$group1),geom ="polygon", type = "norm",
               alpha=0, level=0.90, show.legend = FALSE, linetype = 'dashed', linewidth=0.5)+
  theme_classic() +
  labs(x = "t-SNE 1", y = "t-SNE 2",colour = "Cluster") +
  scale_color_manual(values = c("#FF0000", "#1874CD", "#cc33CD"), 
                     labels = c("Cluster1", "Cluster2", "Cluster3")) +
  theme(legend.position = c( ))
ggplot(tsne_res, aes(x = tSNE2, y = tSNE1, colour = clinical_OS$group1)) +
  geom_point(size = 2, shape = 16) +
  stat_ellipse(aes(fill=clinical_OS$group1,colour=clinical_OS$group1),geom ="polygon", type = "norm",
               alpha=0, level=0.90, show.legend = FALSE, linetype = 'dashed', linewidth=0.5)+
  stat_ellipse(aes(fill=clinical_OS$group1,colour=clinical_OS$group1),geom ="polygon", type = "norm",
               alpha=0, level=0.85, show.legend = FALSE, linetype = 'dashed', linewidth=0.5)+
  theme_classic() +
  labs(x = "t-SNE 1", y = "t-SNE 2",colour = "Cluster") +
  scale_color_manual(values = c("#b93a38", "#5c7197", "#4f7f58"), 
                     labels = c("Cluster1", "Cluster2", "Cluster3")) +
  theme(legend.position = c( ))
# rownames(tsne_res)<- clinical_OS$group1
#PCA作图,更不好#
pca_result <- prcomp(t(PSI_data_group))
pca_result <- data.frame(Sample_ID = colnames(PSI_data_group), pca_result$x)

pca_result <- merge(pca_result,clinical_OS[,c(1,2)],by.x="Sample_ID",by.y = "sample")
ggplot(pca_result, aes(x = PC1, y = PC2, colour = group1)) +
  geom_point(size = 2, shape = 16) +
  #stat_ellipse(level = 0.99,lwd = 0.8,linetype = 2 ) +
  stat_ellipse(level = 0.95,lwd = 0.8,linetype = 2) +
  theme_classic() +
  labs(x = "PC1", y = "PC2",colour = "Cluster") +
  scale_color_manual(values = c("#FF7744", "#00BBFF", "#00AA55"), 
                     labels = c("Cluster1", "Cluster2", "Cluster3")) +
  theme(legend.position = c(0.85, 0.9),
        legend.title = element_blank())
# setwd('/home/cmx/Others_Analysis/KeH/InputData/')
# count_files = dir("expdata/",pattern = "*.tsv$",recursive = T)
# 
# exp = list()
# for(i in 1:length(count_files)){
#   exp[[i]] = read.table(paste0("expdata/",count_files[[i]]),header=T,sep="\t")
#   exp[[i]] = exp[[i]][-(1:4),]  # The first 4 rows contain unneeded information
#   exp[[i]] = exp[[i]]$unstranded  # the fourth column (unstranded) was the count value
# }
# exp = as.data.frame(do.call(cbind,exp))

# dim(exp)
# #[1] 60660   613
# exp[1:4,1:4]
# #     V1   V2    V3   V4
# # 1 2157 4455 10949 2375
# # 2   26   13    33   13
# # 3  978 1610  1621 1264
# # 4  604  831   462  764
# 
# #TCGA ID and file name match
# meta = jsonlite::fromJSON("metadata.cart.json")
# ID = sapply(meta$associated_entities,
#             function(x){x$entity_submitter_id})
# file2id = data.frame(file_name = meta$file_name,
#                      ID = ID)
# 
# count_files2 = stringr::str_split(count_files,"/",simplify = T)[,2]
# table(count_files2 %in% file2id$file_name)
# # TRUE 
# #  613 
# 
# file2id = file2id[match(count_files2,file2id$file_name),]
# identical(file2id$file_name,count_files2)
# #[1] TRUE
# 
# colnames(exp) = file2id$ID
# 
# gene_name = data.table::fread(paste0("expdata/",count_files[1]))$gene_name
# gene_name = gene_name[-seq(1,4)] # The first 4 rows contain unneeded information
# exp = cbind(gene_name=gene_name,exp)
# dim(exp)
# #[1] 60660   614
# exp = exp[!duplicated(exp$gene_name),]
# rownames(exp) = exp$gene_name
# exp = exp[,-1]
# dim(exp)
# #[1] 59427   613
# 
# #gene filter
# exp = exp[apply(exp, 1, function(x) sum(x > 0) > 0.5*ncol(exp)), ] # only keep genes express in half of samples
# dim(exp)
# #[1] 32809   613
# 
# #group information
# table(str_sub(colnames(exp),14,15))
# # Assuming 'exp' is your expression matrix
# tumor_columns <- colnames(exp)[as.numeric(str_sub(colnames(exp), 14, 15)) < 10]
# normal_columns <- colnames(exp)[as.numeric(str_sub(colnames(exp), 14, 15)) > 10]
# exp_normal<- exp[, !colnames(exp) %in% tumor_columns ]
# # Remove columns with 'normal' group
# exp <- exp[, !colnames(exp) %in% normal_columns]
# data<-exp
# # TCGA中有前面12个字母重复的样本，所以需要将这一部分去除
# duplicated <- substr(colnames(data),1,12)[which(duplicated(substr(colnames(data),1,12)))]
# unique(duplicated)
# repetition_last <- c()
# for (i in 1:length(unique(duplicated))) {
#   # i <- length(unique(duplicated))[1]
#   repetition <- colnames(data)[which(str_detect(colnames(data), unique(duplicated)[i]))]
#   repetition_last <- cbind(repetition,repetition_last)
# }
# 
# # 保留重复样本中的一个
# data <- data[, !duplicated(substr(colnames(data), 1, 12))]
# exp_normal <- exp_normal[, !duplicated(substr(colnames(exp_normal), 1, 12))]
# setwd('/home/cmx/senes_splicing/')
# 
# # 修改列名，方便与表达数据整合
# colnames(data) <- substr(colnames(data),1,12)
# colnames(exp_normal) <- substr(colnames(exp_normal),1,12)
# 
# data_log2<-log2(data+1)
# old_patient_count<-data[,colnames(data) %in% clinical_OS$sample ]
# old_patient_count_log<-data_log2[,colnames(data_log2) %in% clinical_OS$sample ]
#save(old_patient_count,old_patient_count_log,file = "./outdata/step2_count_log2.Rdata")
# load("./outdata/step2_count_log2.Rdata")
# load("/home/cmx/Pan-cancer/outdata/exp/BRCA.Rdata")####之前复现课题里下载的泛癌数据
# tumor_TPM_log2<-log2(tumor_TPM+1)
# old_patient_tpm<-tumor_TPM[,colnames(tumor_TPM) %in% clinical_OS$sample ]
# old_patient_tpm_log<-tumor_TPM_log2[,colnames(tumor_TPM_log2) %in% clinical_OS$sample ]
# save(old_patient_tpm,old_patient_tpm_log,file = "./outdata/step2_tpm_log2.Rdata")
# load("./outdata/step2_tpm_log2.Rdata")

# write.table(old_patient_tpm_log,"./outdata/FPKM_tumor_log2.txt",sep = "\t",
#             row.names = T,col.names = NA,quote = F)
# 
# ##转换GCT表达谱格式
# filterCommonGenes(input.f = "./outdata/FPKM_tumor_log2.txt",   #输入文件名
#                   output.f = "./outdata/FPKM_tumor_log2.gct",   #输出文件名
#                   id = "GeneSymbol")   #行名为gene symbol
# ##进行estimate分析
# estimateScore("./outdata/FPKM_tumor_log2.gct",   
#               "./outdata/FPKM_tumor_estimate_score.txt",  
#               platform="affymetrix")   #默认平台
est <- read.table("./outdata/FPKM_tumor_estimate_score.txt", 
                  sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
est <- est[,-1]   #移除第一列
colnames(est) <- est[1,]   #设置列名
est <- as.data.frame(t(est[-1,]))
rownames(est) <- gsub("\\.","-", rownames(est))
est <- est %>%
  mutate(across(where(is.character), as.numeric))
est$sample <- rownames(est)## 行名转换为第一列
est<-merge(est,group1,by="sample")###用来tsne的数据
colnames(est)
# est<-est[!est$group1=="Cluster3",]
# clinical_OS$OS.time<-ifelse(clinical_OS$vital_status=='Alive',clinical_OS$days_to_last_followup,clinical_OS$days_to_death)
# est$group1<-ifelse(est$group1=="Cluster2","High","Low")
####2.Immunescore####
# list<-c("TCGA-BH-A18G","TCGA-AC-A23H","TCGA-AN-A046","TCGA-AC-A5XS","TCGA-A2-A0EV","TCGA-BH-A0HP","TCGA-E9-A1R3",
#         "TCGA-C8-A1HM","TCGA-AC-A2B8","TCGA-LD-A74U","TCGA-EW-A1P5","TCGA-D8-A1J8",
#         "TCGA-B6-A0RQ")
# est <- est[ !est$sample %in% list,]
ggplot(est,aes(x=factor(group1),y=ImmuneScore,fill=factor(group1), color=factor(group1)))+
  geom_violin(trim=FALSE)+  
  scale_fill_manual(values = c("#388E8E", "#EE4000","#aa1000"))+
  scale_color_manual(values=c("black","black","black"))+
  geom_boxplot(fill="white", width=0.2, size=0.5)+
  stat_compare_means(comparisons = list(c("Cluster1", "Cluster2"),c("Cluster2", "Cluster3"),c("Cluster1", "Cluster3")))+
  # stat_compare_means(comparisons = list(c("High","Low")))+
  # stat_compare_means(method="wilcox.test",label.y =4000,label.x =1.25) +
  ggtitle("")+
  theme_classic()+
  theme(legend.position = "None")+
  theme(axis.text.x=element_text(colour="black",family="Times",size=14),
        axis.text.y=element_text(family="Times",size=14,face="plain"),
        axis.title.y=element_text(family="Times",size = 14,face="plain"),
        axis.title.x=element_text(family="Times",size = 14,face="plain"),
        plot.title = element_text(family="Times",size=15,face="bold",hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ylab("ImmuneScore")+xlab("Subgroup")
colnames(est)
####3.StromalScore####
ggplot(est,aes(x=factor(group1),y=StromalScore,fill=factor(group1), color=factor(group1)))+
  geom_violin(trim=FALSE)+  
  scale_fill_manual(values = c("#388E8E", "#EE4000","#aa1000"))+
  scale_color_manual(values=c("black","black","black"))+
  geom_boxplot(fill="white", width=0.2, size=0.5)+
  stat_compare_means(comparisons = list(c("Cluster1", "Cluster2"),c("Cluster2", "Cluster3"),c("Cluster1", "Cluster3")))+
  # stat_compare_means(comparisons = list(c("High","Low")))+
  #stat_compare_means(method="wilcox.test",label.y =4000,label.x =1.25) +
  ggtitle("")+
  theme_classic()+
  theme(legend.position = "None")+
  theme(axis.text.x=element_text(colour="black",family="Times",size=14),
        axis.text.y=element_text(family="Times",size=14,face="plain"),
        axis.title.y=element_text(family="Times",size = 14,face="plain"),
        axis.title.x=element_text(family="Times",size = 14,face="plain"),
        plot.title = element_text(family="Times",size=15,face="bold",hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ylab("StromalScore")+xlab("Subgroup")
####4.ESTIMATEScore####
ggplot(est,aes(x=factor(group1),y=ESTIMATEScore,fill=factor(group1), color=factor(group1)))+
  geom_violin(trim=FALSE)+  
  scale_fill_manual(values = c("#388E8E", "#EE4000","#aa1000"))+
  scale_color_manual(values=c("black","black","black"))+
  geom_boxplot(fill="white", width=0.2, size=0.5)+
  stat_compare_means(comparisons = list(c("Cluster1", "Cluster2"),c("Cluster2", "Cluster3"),c("Cluster1", "Cluster3")))+
  #stat_compare_means(method="wilcox.test",label.y =4000,label.x =1.25) +
  ggtitle("")+
  theme_classic()+
  theme(legend.position = "None")+
  theme(axis.text.x=element_text(colour="black",family="Times",size=14),
        axis.text.y=element_text(family="Times",size=14,face="plain"),
        axis.title.y=element_text(family="Times",size = 14,face="plain"),
        axis.title.x=element_text(family="Times",size = 14,face="plain"),
        plot.title = element_text(family="Times",size=15,face="bold",hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ylab("ESTIMATEScore")+xlab("Subgroup")
####TumorPurity####
ggplot(est,aes(x=factor(group1),y=TumorPurity,fill=factor(group1), color=factor(group1)))+
  geom_violin(trim=FALSE)+  
  scale_fill_manual(values = c("#388E8E", "#EE4000","#aa1000"))+
  scale_color_manual(values=c("black","black","black"))+
  geom_boxplot(fill="white", width=0.2, size=0.5)+
  stat_compare_means(comparisons = list(c("Cluster1", "Cluster2"),c("Cluster2", "Cluster3"),c("Cluster1", "Cluster3")))+
  #stat_compare_means(method="wilcox.test",label.y =4000,label.x =1.25) +
  ggtitle("")+
  theme_classic()+
  theme(legend.position = "None")+
  theme(axis.text.x=element_text(colour="black",family="Times",size=14),
        axis.text.y=element_text(family="Times",size=14,face="plain"),
        axis.title.y=element_text(family="Times",size = 14,face="plain"),
        axis.title.x=element_text(family="Times",size = 14,face="plain"),
        plot.title = element_text(family="Times",size=15,face="bold",hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ylab("TumorPurity")+xlab("Subgroup")