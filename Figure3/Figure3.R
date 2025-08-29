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
# setwd('/home/cmx/Others_Analysis/KeH/InputData/')
# count_files = dir("expdata/",pattern = "*.tsv$",recursive = T)

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

####三个亚型的区别####
# clinical_stage <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")
# old_patient_stage <- clinical_stage%>% as.data.frame() %>%
#   ###挑选些想要研究的信息，这里我需要的信息大概就这些
#   dplyr::select(submitter_id,
#                 age_at_index,
#                 vital_status,
#                 ajcc_pathologic_stage,
#                 days_to_death,
#                 days_to_last_follow_up,
#                 ajcc_pathologic_t,
#                 ajcc_pathologic_m,
#                 ajcc_pathologic_n) %>%
#   ###重命名得到的这些挑好的列名
#   dplyr::rename(Sample_ID=submitter_id,
#                 Age=age_at_index,
#                 Status=vital_status,
#                 Stage=ajcc_pathologic_stage,
#                 daystodeath=days_to_death,
#                 daystolastfollowup = days_to_last_follow_up,
#                 Mstage = ajcc_pathologic_m,
#                 Nstage = ajcc_pathologic_n,
#                 Tstage = ajcc_pathologic_t)
# old_patient_stage<-old_patient_stage[old_patient_stage$Sample_ID %in% clinical_OS$sample, ]
# 
# ##有些临床信息特别是分期相关的，会有些很细的分型，没有必要，所以就直接把他们简称一下方便统一
# old_patient_stage  <- old_patient_stage  %>% dplyr::mutate(Tstage = case_when(
#   grepl("T0",Tstage) ~ "T0",
#   grepl("T1",Tstage) ~ "T1",
#   grepl("T2",Tstage) ~ "T2",
#   grepl("T3",Tstage) ~ "T3",
#   grepl("T4",Tstage) ~ "T4",
#   T ~ NA_character_
# ))
# old_patient_stage  <- old_patient_stage  %>% dplyr::mutate(Mstage = case_when(
#   grepl("M0",Mstage) ~ "M0",
#   grepl("M1",Mstage) ~ "M1",
#   grepl("MX",Mstage) ~ "MX",
#   T ~ NA_character_
# ))
# 
# old_patient_stage  <- old_patient_stage  %>% dplyr::mutate(Nstage = case_when(
#   grepl("N0",Nstage) ~ "N0",
#   grepl("N1",Nstage) ~ "N1",
#   grepl("N2",Nstage) ~ "N2",
#   grepl("N3",Nstage) ~ "N3",
#   grepl("NX",Nstage) ~ "NX",
#   T ~ NA_character_
# ))
# old_patient_stage$Stage<- ifelse(old_patient_stage$Stage== "Stage I" ,"Stage IC",old_patient_stage$Stage)
# old_patient_stage  <- old_patient_stage  %>% dplyr::mutate(Stage = case_when(
#   grepl("Stage IC",Stage) ~ "Stage I",
#   grepl("Stage IA",Stage) ~ "Stage I",
#   grepl("Stage IB",Stage) ~ "Stage I",
#   grepl("Stage IIA",Stage) ~ "Stage II",
#   grepl("Stage IIB",Stage) ~ "Stage II",
#   grepl("Stage IIIA",Stage) ~ "Stage III",
#   grepl("Stage IIIB",Stage) ~ "Stage III",
#   grepl("Stage IIIC",Stage) ~ "Stage III",
#   grepl("Stage IV",Stage) ~ "Stage IV",
#   grepl("Stage X",Stage) ~ "Stage X",
#   
#   T ~ NA_character_
# ))
# 
# ###把daystodeath和lastFollowupTime合并一下
# old_patient_stage$OS.time<-ifelse(old_patient_stage$Status=='Alive',old_patient_stage$daystolastfollowup,old_patient_stage$daystodeath)
# old_patient_stage$OS<- ifelse(old_patient_stage$Status=="Alive",0,1)
# ###重复的行删掉,得到最后的临床信息
# old_patient_stage<-old_patient_stage[!duplicated(old_patient_stage),]
# old_patient_stage<-merge(old_patient_stage,group1,by.x="Sample_ID",by.y="sample")
# old_patient_stage<- old_patient_stage[order(old_patient_stage$group1, decreasing = F), ]
# 
# table(old_patient_stage$Tstage)
# table(old_patient_stage$Mstage)
# table(old_patient_stage$Nstage)
# table(old_patient_stage$Stage)
load("./outdata/old_patient_stage.Rdata")

# save(old_patient_stage,file = "./outdata/old_patient_stage.Rdata")
imm28<-read.csv("/home/cmx/RNA-splicing/indata/28immcell.csv",header = TRUE,sep = ",")
table(imm28$Cell.type)
gmt <- list()
for (i in 1:28) {
  gmt[[i]] <- imm28[which(imm28$Num == i),]$Metagene
}
names(gmt) <- c("Activated B cell","Activated CD4 T cell","Activated CD8 T cell","Activated dendritic cell",
                "CD56bright natural killer cell","CD56dim natural killer cell","Central memory CD4 T cell",
                "Central memory CD8 T cell","Effector memeory CD4 T cell","Effector memeory CD8 T cell","Eosinophil",
                "Gamma delta T cell","Immature  B cell","Immature dendritic cell","Macrophage","Mast cell","MDSC",
                "Memory B cell","Monocyte","Natural killer cell","Natural killer T cell","Neutrophil","Plasmacytoid dendritic cell",
                "Regulatory T cell","T follicular helper cell","Type 1 T helper cell","Type 17 T helper cell ","Type 2 T helper cell")
#save(gmt,file ="outdata/gmt.Rdata") #保存
TCGA_BRCA_matrix<-as.matrix(old_patient_tpm_log)
# save(TCGA_BRCA_matrix,gmt,file = "./outdata/GSVA_MY.Rdata")
# ssGSEA_matrix <- gsva(expr = TCGA_BRCA_matrix,
#                       gset.idx.list = gmt,
#                       method = 'ssgsea',
#                       kcdf = "Gaussian",
#                       abs.ranking = TRUE)##自己电脑上跑的，直接加载
# load("./indata/GSVA_out.Rdata")
imm_ssgsea<-as.data.frame(ssGSEA_matrix)
scale_imm_ssgsea<-as.matrix(t(scale(t(imm_ssgsea))))
scale_imm_ssgsea <- scale_imm_ssgsea[, order(match(colnames(scale_imm_ssgsea), old_patient_stage$Sample_ID))]
# old_patient_stage<-old_patient_stage[!old_patient_stage$group1=="Cluster1",]
# clinical_OS$OS.time<-ifelse(clinical_OS$vital_status=='Alive',clinical_OS$days_to_last_followup,clinical_OS$days_to_death)
# old_patient_stage$group1<-ifelse(old_patient_stage$group1=="Cluster2","High","Low")
list<-c("TCGA-BH-A18G","TCGA-AC-A23H","TCGA-AN-A046","TCGA-AC-A5XS","TCGA-A2-A0EV","TCGA-BH-A0HP","TCGA-E9-A1R3",
        "TCGA-C8-A1HM","TCGA-AC-A2B8","TCGA-LD-A74U","TCGA-EW-A1P5","TCGA-D8-A1J8",
        "TCGA-B6-A0RQ")
scale_imm_ssgsea <- scale_imm_ssgsea[,!colnames(scale_imm_ssgsea) %in% list ]
old_patient_stage<-old_patient_stage[!old_patient_stage$Sample_ID%in% list,]
sfd = HeatmapAnnotation(N = old_patient_stage$Nstage,
                        M = old_patient_stage$Mstage,
                        "T" = old_patient_stage$Tstage,
                        Stage = old_patient_stage$Stage,
                        Age = old_patient_stage$Age,
                        Subgroup = old_patient_stage$group1)
# scale_imm_ssgsea <- scale_imm_ssgsea[,old_patient_stage$Sample_ID]

# 绘制热图
Heatmap((scale_imm_ssgsea),
        col = colorRampPalette(c("navy","white","firebrick3"))(100),
        column_split  = factor(old_patient_stage$group1, levels = c("Cluster1","Cluster2","Cluster3")),
        name = " ",
        top_annotation = sfd,
        show_column_names = F,
        cluster_rows = F,
        cluster_columns = F,
        column_title = ""
)
#简单热图的方法#
annotation_col <- data.frame(Subtype = old_patient_stage$group1)
row.names(annotation_col) <- colnames(scale_imm_ssgsea)
#colnames(annotation_row) <- " 
#HLA_KIR[HLA_KIR > 2] <- 2
#HLA_KIR[HLA_KIR < -2] <- -2                
pheatmap::pheatmap(scale_imm_ssgsea,
                   scale = "row",
                   col = colorRampPalette(c("navy","white","firebrick3"))(200),
                   show_colnames = F,
                   cluster_rows = F,
                   cluster_cols = F,
                   #gaps_col = c(275),
                   annotation_col = annotation_col)

#加上P值#
load("./indata/GSVA_out.Rdata")
imm_ssgsea<-as.data.frame(ssGSEA_matrix)
HLA_KIR1<-imm_ssgsea
HLA_KIR1 <- as.data.frame(t(HLA_KIR1))
HLA_KIR1$Sample_ID <- rownames(HLA_KIR1)
HLA_KIR1<-merge(HLA_KIR1,old_patient_stage,by="Sample_ID")
colnames(HLA_KIR1)[41]<-"Cluster"

p_value <- data.frame(gene = rownames(data.frame(imm_ssgsea)))

HLA_KIR1$Cluster <- factor(HLA_KIR1$Cluster,levels = c("Cluster1","Cluster2","Cluster3"))

anno_pvalue <- data.frame()
for(i in p_value$gene){
  data <- dplyr::select(HLA_KIR1,gene = i,Cluster = Cluster)
  out_results <- compare_means(gene~Cluster ,data = data, method ="kruskal.test")
  out_results$.y. <- i
  anno_pvalue <- rbind(anno_pvalue,out_results)
}


####5.加上P值，复杂热图重新画图，28个免疫细胞丰度差别####
load("./indata/GSVA_out.Rdata")
imm_ssgsea<-as.data.frame(ssGSEA_matrix)
HLA_gene <- c("Activated B cell","Activated CD4 T cell","Activated CD8 T cell","Activated dendritic cell",
              "CD56bright natural killer cell","CD56dim natural killer cell","Central memory CD4 T cell",
              "Central memory CD8 T cell","Effector memeory CD4 T cell","Effector memeory CD8 T cell","Eosinophil",
              "Gamma delta T cell","Immature  B cell","Immature dendritic cell","Macrophage","Mast cell","MDSC",
              "Memory B cell","Monocyte","Natural killer cell","Natural killer T cell","Neutrophil","Plasmacytoid dendritic cell",
              "Regulatory T cell","T follicular helper cell","Type 1 T helper cell","Type 17 T helper cell ","Type 2 T helper cell")
HLA_KIR<-imm_ssgsea[HLA_gene,]
HLA_KIR<-as.matrix(t(scale(t(HLA_KIR))))

# HLA_KIR[HLA_KIR > 2] <- 2
# HLA_KIR[HLA_KIR < -2] <- -2
HLA_KIR <- HLA_KIR[, order(match(colnames(HLA_KIR), old_patient_stage$Sample_ID))]
HLA_KIR<-HLA_KIR[,!colnames(HLA_KIR)%in% list]

sfd = HeatmapAnnotation(N = old_patient_stage$Nstage,
                        M = old_patient_stage$Mstage,
                        "T" = old_patient_stage$Tstage,
                        Stage = old_patient_stage$Stage,
                        Age = old_patient_stage$Age,
                        Subgroup = old_patient_stage$group1)
row_ha2 = rowAnnotation(foo = anno_text(anno_pvalue$p.signif  , location = 0.5, just = "center"))

# 绘制热图
ComplexHeatmap::Heatmap(HLA_KIR,
                        col = colorRampPalette(c("navy","white","firebrick3"))(200),
                        column_split  = factor(old_patient_stage$group1, levels = c("Cluster1","Cluster2","Cluster3")),
                        name = " ",
                        row_names_side = "left",
                        right_annotation = row_ha2,
                        top_annotation = sfd,
                        show_column_names = F,
                        cluster_rows = F,
                        cluster_columns = F,
)

####箱线图展示#####
load("./indata/GSVA_out.Rdata")
resDD<-as.data.frame(ssGSEA_matrix)
resDD<-as.data.frame(t(resDD))
resDD$sample<-rownames(resDD)
data2<-merge(resDD,clinical_OS[,c(1,2)],by="sample")
data2 = dplyr::select(data2,30,everything())
data2 = dplyr::select(data2,2,everything())
rownames(data2) <- data2[,1]## 第一列转换为行名
#data2 <- data2[,-1] ## 删除第一列

res <- pivot_longer(data = data2,
                    cols = 3:ncol(data2),
                    names_to = "gene",
                    values_to = "value",
                    values_drop_na = TRUE)
res$value40<-res$value*80
ggboxplot(res, x = "gene", y = "value40", color = "black",outlier.size=0.2,
          outlier.colour="black",fill="group1",ylab = "Scale of fraction",
          palette =c("#CD3333", "#8DB6CD","#458B00"),)+
  theme(legend.position = "top",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60,vjust = 0.78,hjust =0.78 ,size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "#000000"),)+
  stat_compare_means(aes(group=group1),method="kruskal.test",
                     label="p.signif")
#免疫因子#
genes <- c("CD276","NRP1","CD200","CD40","VSIR","CD70","TNFSF9","CD160","CD40LG","TNFRSF18","TMIGD2",
           "CD80","MTAP","LAIR1","BTLA","CD28","TNFSF4","CD244","CD200R1","LAG3","PDCD1","TNFSF15","HAVCR2",
           "BTNL2","TNFSF18","CD274","PDCD1LG2","ICOSLG","ADORA2A","HHLA2")
TCGA_KIRC1<-old_patient_count_log[genes,]
TCGA_KIRC1<-as.matrix(t(scale(t(TCGA_KIRC1))))

TCGA_KIRC1 <- TCGA_KIRC1[, order(match(colnames(TCGA_KIRC1), old_patient_stage$Sample_ID))]

sfd = HeatmapAnnotation(Subtype = old_patient_stage$group1)


# 绘制热图
Heatmap(TCGA_KIRC1,
        col = colorRampPalette(c("navy","white","firebrick3"))(100),
        column_split  = factor(old_patient_stage$group1, levels = c("Cluster1","Cluster2","Cluster3")),
        name = " ",
        top_annotation = sfd,
        show_column_names = F,
        cluster_rows = T,
        cluster_columns = F,
        column_title = ""
)
#简单热图的方法#
annotation_col <- data.frame(Subtype = old_patient_stage$group1)
row.names(annotation_col) <- colnames(TCGA_KIRC1)
#colnames(annotation_row) <- " 

#TCGA_KIRC1[TCGA_KIRC1 > 2] <- 2
#TCGA_KIRC1[TCGA_KIRC1 < -2] <- -2
pheatmap::pheatmap(TCGA_KIRC1,
                   scale = "row",
                   col = colorRampPalette(c("navy","white","firebrick3"))(200),
                   show_colnames = F,
                   cluster_rows = T,
                   cluster_cols = F,
                   #gaps_col = c(275),
                   annotation_col = annotation_col)

####6.加上P值，复杂热图重新画图，免疫因子差别####
genes <- c("CD276","NRP1","CD200","CD40","VSIR","CD70","TNFSF9","CD160","CD40LG","TNFRSF18","TMIGD2",
           "CD80","MTAP","LAIR1","BTLA","CD28","TNFSF4","CD244","CD200R1","LAG3","PDCD1","TNFSF15","HAVCR2",
           "BTNL2","TNFSF18","CD274","PDCD1LG2","ICOSLG","ADORA2A","HHLA2")
TCGA_KIRC1<-old_patient_count_log[genes,]

p_value <- data.frame(gene = rownames(data.frame(TCGA_KIRC1)))

TCGA_KIRC1 <- as.data.frame(t(TCGA_KIRC1))
TCGA_KIRC1$Sample_ID <- rownames(TCGA_KIRC1)
TCGA_KIRC1<-merge(TCGA_KIRC1,old_patient_stage,by="Sample_ID")
colnames(TCGA_KIRC1)[42]<-"Cluster"


TCGA_KIRC1$Cluster <- factor(TCGA_KIRC1$Cluster,levels = c("Cluster1","Cluster2","Cluster3"))

anno_pvalue <- data.frame()
for(i in p_value$gene){
  data <- dplyr::select(TCGA_KIRC1,gene = i,Cluster = Cluster)
  out_results <- compare_means(gene~Cluster ,data = data, method ="kruskal.test")
  out_results$.y. <- i
  anno_pvalue <- rbind(anno_pvalue,out_results)
}


##重新画图
genes <- c("CD276","NRP1","CD200","CD40","VSIR","CD70","TNFSF9","CD160","CD40LG","TNFRSF18","TMIGD2",
           "CD80","MTAP","LAIR1","BTLA","CD28","TNFSF4","CD244","CD200R1","LAG3","PDCD1","TNFSF15","HAVCR2",
           "BTNL2","TNFSF18","CD274","PDCD1LG2","ICOSLG","ADORA2A","HHLA2")
TCGA_KIRC1<-old_patient_count_log[genes,]
TCGA_KIRC1<-as.matrix(t(scale(t(TCGA_KIRC1))))
TCGA_KIRC1 <- TCGA_KIRC1[, order(match(colnames(TCGA_KIRC1), old_patient_stage$Sample_ID))]
# TCGA_KIRC1[TCGA_KIRC1 > 2] <- 2
# TCGA_KIRC1[TCGA_KIRC1 < -2] <- -2

sfd = HeatmapAnnotation(Subtype = old_patient_stage$group1)

row_ha2 = rowAnnotation(foo = anno_text(anno_pvalue$p.signif  , location = 0.5, just = "center"))

# 绘制热图
ComplexHeatmap::Heatmap(TCGA_KIRC1,
                        col = colorRampPalette(c("navy","white","firebrick3"))(200),
                        column_split  = factor(old_patient_stage$group1, levels = c("Cluster1","Cluster2","Cluster3")),
                        name = " ",
                        row_names_side = "left",
                        right_annotation = row_ha2,
                        top_annotation = sfd,
                        show_column_names = F,
                        cluster_rows = F,
                        cluster_columns = F,
)

#HLA#
HLA_gene <- c("HLA-E","HLA-DPB2","HLA-C","HLA-DQB1","HLA-DQA1","HLA-DMA","HLA-DRB1","HLA-H","HLA-DRB5",
              "HLA-DOA","HLA-DPB1","HLA-DRA","HLA-DMB","HLA-DPA1")
HLA_KIR<-old_patient_count_log[HLA_gene,]
HLA_KIR<-as.matrix(t(scale(t(HLA_KIR))))
HLA_KIR <- HLA_KIR[, !colnames(HLA_KIR) %in% list]

#HLA_KIR[HLA_KIR > 2] <- 2
#HLA_KIR[HLA_KIR < -2] <- -2


HLA_KIR <- HLA_KIR[, order(match(colnames(HLA_KIR), old_patient_stage$Sample_ID))]
sfd = HeatmapAnnotation(Subtype = old_patient_stage$group1)

# 绘制热图
Heatmap(HLA_KIR,
        col = colorRampPalette(c("navy","white","firebrick3"))(200),
        column_split  = factor(old_patient_stage$group1, levels = c("Cluster1","Cluster2","Cluster3")),
        name = " ",
        top_annotation = sfd,
        show_column_names = F,
        cluster_rows = F,
        cluster_columns = F,
        column_title = ""
)
#简单热图的方法#
annotation_col <- data.frame(Subtype = old_patient_stage$group1)
row.names(annotation_col) <- colnames(HLA_KIR)
#colnames(annotation_row) <- " 
#HLA_KIR[HLA_KIR > 2] <- 2
#HLA_KIR[HLA_KIR < -2] <- -2                
pheatmap::pheatmap(HLA_KIR,
                   scale = "row",
                   col = colorRampPalette(c("navy","white","firebrick3"))(200),
                   show_colnames = F,
                   cluster_rows = F,
                   cluster_cols = F,
                   #gaps_col = c(275),
                   annotation_col = annotation_col)

####7.加上P值，复杂热图重新画图，HLA细胞丰度差别####
HLA_gene <- c("HLA-E","HLA-DPB2","HLA-C","HLA-DQB1","HLA-DQA1","HLA-DMA","HLA-DRB1","HLA-H","HLA-DRB5",
              "HLA-DOA","HLA-DPB1","HLA-DRA","HLA-DMB","HLA-DPA1")
HLA_KIR1<-old_patient_count_log[HLA_gene,]
HLA_KIR1 <- as.data.frame(t(HLA_KIR1))
HLA_KIR1$Sample_ID <- rownames(HLA_KIR1)
HLA_KIR1<-merge(HLA_KIR1,old_patient_stage,by="Sample_ID")
colnames(HLA_KIR1)[27]<-"Cluster"

p_value <- data.frame(gene = rownames(data.frame(HLA_KIR)))

HLA_KIR1$Cluster <- factor(HLA_KIR1$Cluster,levels = c("Cluster1","Cluster2","Cluster3"))

anno_pvalue <- data.frame()
for(i in p_value$gene){
  data <- dplyr::select(HLA_KIR1,gene = i,Cluster = Cluster)
  out_results <- compare_means(gene~Cluster ,data = data, method ="kruskal.test")
  out_results$.y. <- i
  anno_pvalue <- rbind(anno_pvalue,out_results)
}


# 重新画图 #
HLA_gene <- c("HLA-E","HLA-DPB2","HLA-C","HLA-DQB1","HLA-DQA1","HLA-DMA","HLA-DRB1","HLA-H","HLA-DRB5",
              "HLA-DOA","HLA-DPB1","HLA-DRA","HLA-DMB","HLA-DPA1")
HLA_KIR<-old_patient_count_log[HLA_gene,]
HLA_KIR<-as.matrix(t(scale(t(HLA_KIR))))

# HLA_KIR[HLA_KIR > 2] <- 2
# HLA_KIR[HLA_KIR < -2] <- -2


HLA_KIR <- HLA_KIR[, order(match(colnames(HLA_KIR), old_patient_stage$Sample_ID))]
sfd = HeatmapAnnotation(Subtype = old_patient_stage$group1)
row_ha2 = rowAnnotation(foo = anno_text(anno_pvalue$p.signif  , location = 0.5, just = "center"))
HLA_KIR <- HLA_KIR[, !colnames(HLA_KIR) %in% list]

# 绘制热图
ComplexHeatmap::Heatmap(HLA_KIR,
                        col = colorRampPalette(c("navy","white","firebrick3"))(200),
                        column_split  = factor(old_patient_stage$group1, levels = c("Cluster1","Cluster2","Cluster3")),
                        name = " ",
                        row_names_side = "left",
                        right_annotation = row_ha2,
                        top_annotation = sfd,
                        show_column_names = F,
                        cluster_rows = F,
                        cluster_columns = F,
)
####8.其他15个signature#####
D3<-read_csv("/home/cmx/Reproduce_and_become_mine/BLCA_Oncogene/InputData/3D.csv")
D3<-D3[,1:3]
gmt <- list()
for (i in 1:15) {
  gmt[[i]] <- D3[which(D3$Type == i),]$GeneSymbol
}
names(gmt) <- c("CD 8 T effector","Immune Checkpoint","Antigen processing machinery","Mismatch Repair","Nucleotide excision repair",
                "DNA damage response","DNA replication","Base excision repair","Pan-F-TBRs","EMT1","EMT2","EMT3",
                "TMEscoreA","TMEscoreB","TMEscore")
#save(gmt,file ="./outdata/gmt.Rdata") #保存
# load("/home/cmx/Reproduce_and_become_mine/BLCA_Oncogene/OutData/gmt.Rdata")
# library(GSVA)
# load("./gmt.Rdata")
# gmt1<-gmt
# load("./GSVA_MY.Rdata")
# expr_matrix<-TCGA_BRCA_matrix
# ssGSEA_matrix <- gsva(expr = TCGA_BRCA_matrix,
#                       gset.idx.list = gmt1,
#                       method = 'ssgsea',
#                       kcdf = "Gaussian",
#                       abs.ranking = TRUE)##自己电脑上跑的，直接加载
# save(ssGSEA_matrix,file = "./15SET.Rdata")

##自己电脑上跑的，直接加载
load("./indata/15SET.Rdata")
resDD<-as.data.frame(ssGSEA_matrix)
resDD<-as.data.frame(t(resDD))
resDD$sample<-rownames(resDD)
data2<-merge(resDD,clinical_OS[,c(1,2)],by="sample")
data2 = dplyr::select(data2,17,everything())
data2 = dplyr::select(data2,2,everything())
rownames(data2) <- data2[,1]## 第一列转换为行名
#data2 <- data2[,-1] ## 删除第一列
list<-c("TCGA-BH-A18G","TCGA-AC-A23H","TCGA-AN-A046","TCGA-AC-A5XS","TCGA-A2-A0EV","TCGA-BH-A0HP","TCGA-E9-A1R3",
        "TCGA-C8-A1HM","TCGA-AC-A2B8","TCGA-LD-A74U","TCGA-EW-A1P5","TCGA-D8-A1J8",
        "TCGA-B6-A0RQ")
data2 <- data2[ !rownames(data2)%in% list,]
res <- pivot_longer(data = data2,
                    cols = 3:ncol(data2),
                    names_to = "gene",
                    values_to = "value",
                    values_drop_na = TRUE)
res$value40<-res$value*80
ggboxplot(res, x = "gene", y = "value40", color = "black",outlier.size=0.2,
          outlier.colour="black",fill="group1",ylab = "Scale of fraction",
          palette =c("#CD3333", "#8DB6CD","#458B00"),)+
  theme(legend.position = "top",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 60,vjust = 0.78,hjust =0.78 ,size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "#000000"),)+
  stat_compare_means(aes(group=group1),method="kruskal.test",
                     label="p.signif")
write.table(res, file = "/home/cmx/senes_splicing/outdata/CB_data/Fig3E_TIL.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#GSVA#
#save(old_patient_count,file = "./outdata/count.Rdata")
# setwd('D:/')#设置工作路径
# load("./count.Rdata")
# library(msigdbr)
# library(GSVA)
# ## msigdbr包提取下载 一般尝试KEGG和GO做GSVA分析
# ##KEGG
# KEGG_df_all <-  msigdbr(species = "Homo sapiens", # Homo sapiens or Mus musculus
#                         category = "C2",
#                         subcategory = "CP:KEGG") 
# KEGG_df <- dplyr::select(KEGG_df_all,gs_name,gs_exact_source,gene_symbol)
# kegg_list <- split(KEGG_df$gene_symbol, KEGG_df$gs_name) ##按照gs_name给gene_symbol分组
# 
# ##GO
# GO_df_all <- msigdbr(species = "Homo sapiens",
#                      category = "C5")  
# GO_df <- dplyr::select(GO_df_all, gs_name, gene_symbol, gs_exact_source, gs_subcat)
# GO_df <- GO_df[GO_df$gs_subcat!="HPO",]
# go_list <- split(GO_df$gene_symbol, GO_df$gs_name) ##按照gs_name给gene_symbol分组
# 
# #GSVA算法需要处理logCPM, logRPKM,logTPM数据或counts数据的矩阵#
# dat <- as.matrix(old_patient_count)
# #dat <- as.matrix(log2(edgeR::cpm(counts))+1)
# #dat <- as.matrix(log2(tpm+1))
# 
# geneset <- kegg_list
# 
# gsva_mat <- gsva(expr=dat, 
#                  gset.idx.list=geneset, 
#                  kcdf="Poisson" ,#"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
#                  verbose=T, 
#                  parallel.sz = parallel::detectCores())#调用所有核
# 
# geneset <- go_list
# 
# gsva_mat_GO <- gsva(expr=dat, 
#                     gset.idx.list=geneset, 
#                     kcdf="Poisson" ,#"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
#                     verbose=T, 
#                     parallel.sz = parallel::detectCores())#调用所有核
#电脑上运行完导入到服务器
####9.KEGG画图-GSVA####
load("./indata/GSVA_KEGG_GO.Rdata")
gsva_mat <- gsva_mat[, order(match(colnames(gsva_mat), old_patient_stage$Sample_ID))]

sfd = HeatmapAnnotation(Subgroup = old_patient_stage$group1)
# 绘制热图
Heatmap(gsva_mat,
        col = colorRampPalette(c("navy","white","firebrick3"))(100),
        column_split  = factor(old_patient_stage$group1, levels = c("Cluster1","Cluster2","Cluster3")),
        name = " ",
        top_annotation = sfd,
        show_column_names = F,
        show_row_names = F,
        cluster_rows = T,
        cluster_columns = F,
        column_title = ""
)
#简单热图的方法##
annotation_col <- data.frame(Subtype = old_patient_stage$group1)
row.names(annotation_col) <- colnames(gsva_mat)
#colnames(annotation_row) <- " 
#HLA_KIR[HLA_KIR > 2] <- 2
#HLA_KIR[HLA_KIR < -2] <- -2                
pheatmap::pheatmap(gsva_mat,
                   scale = "row",
                   col = colorRampPalette(c("navy","white","firebrick3"))(200),
                   show_colnames = F,
                   show_rownames = F,
                   cluster_rows = F,
                   cluster_cols = F,
                   #gaps_col = c(275),
                   annotation_col = annotation_col)

#加上P值##
gsva_mat_p<-as.data.frame(gsva_mat)
gsva_mat_p1<-gsva_mat_p
gsva_mat_p1 <- gsva_mat_p1[, order(match(colnames(gsva_mat_p1), old_patient_stage$Sample_ID))]
gsva_mat_p1 <- as.data.frame(t(gsva_mat_p1))
gsva_mat_p1$Sample_ID <- rownames(gsva_mat_p1)
gsva_mat_p1$Cluster<-old_patient_stage$group1

p_value <- data.frame(gene = rownames(data.frame(gsva_mat_p)))

gsva_mat_p1$Cluster <- factor(gsva_mat_p1$Cluster,levels = c("Cluster1","Cluster2","Cluster3"))
anno_pvalue <- data.frame()
for(i in p_value$gene){
  data <- dplyr::select(gsva_mat_p1,gene = i,Cluster = Cluster)
  out_results <- compare_means(gene~Cluster ,data = data, method ="kruskal.test")
  out_results$.y. <- i
  anno_pvalue <- rbind(anno_pvalue,out_results)
}
anno_pvalue<-anno_pvalue[!anno_pvalue$p.signif == 'ns',] #取出ids$symbol不等于空的行
gsva_mat_p1_si<-gsva_mat_p1[,anno_pvalue$.y.]

gsva_mat_p1_si<-as.matrix(t(scale(t(gsva_mat_p1_si))))

# HLA_KIR[HLA_KIR > 2] <- 2
# HLA_KIR[HLA_KIR < -2] <- -2


# gsva_mat_GO_p1_si <- gsva_mat_GO_p1_si[, order(match(colnames(gsva_mat_GO_p1_si), old_patient_stage$Sample_ID))]

sfd = HeatmapAnnotation(Subtype = old_patient_stage$group1)
row_ha2 = rowAnnotation(foo = anno_text(anno_pvalue$p.signif  , location = 0.5, just = "center"))

# 绘制热图
ComplexHeatmap::Heatmap(t(gsva_mat_p1_si),
                        col = colorRampPalette(c("navy","white","firebrick3"))(200),
                        column_split  = factor(old_patient_stage$group1, levels = c("Cluster1","Cluster2","Cluster3")),
                        name = " ",
                        row_names_side = "left",
                        right_annotation = row_ha2,
                        top_annotation = sfd,
                        show_column_names = F,
                        cluster_rows = F,
                        cluster_columns = F)

pheatmap::pheatmap(t(gsva_mat_p1_si),
                   scale = "row",
                   col = colorRampPalette(c("navy","white","firebrick3"))(200),
                   show_colnames = F,
                   show_rownames = F,
                   cluster_rows = T,
                   cluster_cols = F,
                   #gaps_col = c(275),
                   annotation_col = annotation_col)
####10.GO画图-GAVA####
load("./indata/GSVA_KEGG_GO.Rdata")
gsva_mat_GO <- gsva_mat_GO[, order(match(colnames(gsva_mat_GO), old_patient_stage$Sample_ID))]

sfd = HeatmapAnnotation(Subgroup = old_patient_stage$group1)
# 绘制热图
Heatmap(gsva_mat_GO,
        col = colorRampPalette(c("navy","white","firebrick3"))(100),
        column_split  = factor(old_patient_stage$group1, levels = c("Cluster1","Cluster2","Cluster3")),
        name = " ",
        top_annotation = sfd,
        show_column_names = F,
        show_row_names = F,
        cluster_rows = T,
        cluster_columns = F,
        column_title = ""
)
##简单热图的方法
annotation_col <- data.frame(Subtype = old_patient_stage$group1)
row.names(annotation_col) <- colnames(gsva_mat_GO)
#colnames(annotation_row) <- " 
#HLA_KIR[HLA_KIR > 2] <- 2
#HLA_KIR[HLA_KIR < -2] <- -2                
pheatmap::pheatmap(gsva_mat_GO,
                   scale = "row",
                   col = colorRampPalette(c("navy","white","firebrick3"))(200),
                   show_colnames = F,
                   show_rownames = F,
                   cluster_rows = F,
                   cluster_cols = F,
                   #gaps_col = c(275),
                   annotation_col = annotation_col)

#加上P值#
gsva_mat_GO_p<-as.data.frame(gsva_mat_GO)
gsva_mat_GO_p1<-gsva_mat_GO_p
gsva_mat_GO_p1 <- gsva_mat_GO_p1[, order(match(colnames(gsva_mat_GO_p1), old_patient_stage$Sample_ID))]

gsva_mat_GO_p1 <- as.data.frame(t(gsva_mat_GO_p1))
gsva_mat_GO_p1$Sample_ID <- rownames(gsva_mat_GO_p1)
gsva_mat_GO_p1$Cluster<-old_patient_stage$group1

p_value <- data.frame(gene = rownames(data.frame(gsva_mat_GO_p)))

gsva_mat_GO_p1$Cluster <- factor(gsva_mat_GO_p1$Cluster,levels = c("Cluster1","Cluster2","Cluster3"))
anno_pvalue <- data.frame()
for(i in p_value$gene){
  data <- dplyr::select(gsva_mat_GO_p1,gene = i,Cluster = Cluster)
  out_results <- compare_means(gene~Cluster ,data = data, method ="kruskal.test")
  out_results$.y. <- i
  anno_pvalue <- rbind(anno_pvalue,out_results)
}
anno_pvalue<-anno_pvalue[!anno_pvalue$p.signif == 'ns',] #取出ids$symbol不等于空的行


gsva_mat_GO_p1_si<-gsva_mat_GO_p1[,anno_pvalue$.y.]
gsva_mat_GO_p1_si<-as.matrix(t(scale(t(gsva_mat_GO_p1_si))))

# HLA_KIR[HLA_KIR > 2] <- 2
# HLA_KIR[HLA_KIR < -2] <- -2


# gsva_mat_GO_p1_si <- gsva_mat_GO_p1_si[, order(match(colnames(gsva_mat_GO_p1_si), old_patient_stage$Sample_ID))]

sfd = HeatmapAnnotation(Subtype = old_patient_stage$group1)
row_ha2 = rowAnnotation(foo = anno_text(anno_pvalue$p.signif  , location = 0.5, just = "center"))

# 绘制热图
ComplexHeatmap::Heatmap(t(gsva_mat_GO_p1_si),
                        col = colorRampPalette(c("navy","white","firebrick3"))(200),
                        column_split  = factor(old_patient_stage$group1, levels = c("Cluster1","Cluster2","Cluster3")),
                        name = " ",
                        row_names_side = "left",
                        right_annotation = row_ha2,
                        top_annotation = sfd,
                        show_column_names = F,
                        cluster_rows = F,
                        cluster_columns = F)

pheatmap::pheatmap(t(gsva_mat_GO_p1_si),
                   scale = "row",
                   col = colorRampPalette(c("navy","white","firebrick3"))(200),
                   show_colnames = F,
                   show_rownames = F,
                   cluster_rows = T,
                   cluster_cols = F,
                   #gaps_col = c(275),
                   annotation_col = annotation_col)
####hallmark画图-GAVA####
load("./indata/GSVA_hallmark.Rdata")
gsva_mat_hallmark <- gsva_mat_hallmark[, order(match(colnames(gsva_mat_hallmark), old_patient_stage$Sample_ID))]

sfd = HeatmapAnnotation(Subgroup = old_patient_stage$group1)
# 绘制热图
Heatmap(gsva_mat_hallmark,
        col = colorRampPalette(c("navy","white","firebrick3"))(100),
        column_split  = factor(old_patient_stage$group1, levels = c("Cluster1","Cluster2","Cluster3")),
        name = " ",
        top_annotation = sfd,
        show_column_names = F,
        show_row_names = F,
        cluster_rows = T,
        cluster_columns = F,
        column_title = ""
)
##简单热图的方法
annotation_col <- data.frame(Subtype = old_patient_stage$group1)
row.names(annotation_col) <- colnames(gsva_mat_hallmark)
#colnames(annotation_row) <- " 
#HLA_KIR[HLA_KIR > 2] <- 2
#HLA_KIR[HLA_KIR < -2] <- -2                
pheatmap::pheatmap(gsva_mat_hallmark,
                   scale = "row",
                   col = colorRampPalette(c("navy","white","firebrick3"))(200),
                   show_colnames = F,
                   show_rownames = F,
                   cluster_rows = F,
                   cluster_cols = F,
                   #gaps_col = c(275),
                   annotation_col = annotation_col)

#加上P值#
gsva_mat_hallmark_p<-as.data.frame(gsva_mat_hallmark)
gsva_mat_hallmark_p1<-gsva_mat_hallmark_p
gsva_mat_hallmark_p1 <- gsva_mat_hallmark_p1[, order(match(colnames(gsva_mat_hallmark_p1), old_patient_stage$Sample_ID))]

gsva_mat_hallmark_p1 <- as.data.frame(t(gsva_mat_hallmark_p1))
gsva_mat_hallmark_p1$Sample_ID <- rownames(gsva_mat_hallmark_p1)
gsva_mat_hallmark_p1$Cluster<-old_patient_stage$group1

p_value <- data.frame(gene = rownames(data.frame(gsva_mat_hallmark_p)))

gsva_mat_hallmark_p1$Cluster <- factor(gsva_mat_hallmark_p1$Cluster,levels = c("Cluster1","Cluster2","Cluster3"))
anno_pvalue <- data.frame()
for(i in p_value$gene){
  data <- dplyr::select(gsva_mat_hallmark_p1,gene = i,Cluster = Cluster)
  out_results <- compare_means(gene~Cluster ,data = data, method ="kruskal.test")
  out_results$.y. <- i
  anno_pvalue <- rbind(anno_pvalue,out_results)
}
# anno_pvalue<-anno_pvalue[!anno_pvalue$p.signif == 'ns',] #取出ids$symbol不等于空的行

gsva_mat_hallmark_p1_si<-gsva_mat_hallmark_p1[,anno_pvalue$.y.]
gsva_mat_hallmark_p1_si<-as.matrix(t(scale(t(gsva_mat_hallmark_p1_si))))
gsva_mat_hallmark_p1_si<-as.matrix(t(gsva_mat_hallmark_p1_si))

# HLA_KIR[HLA_KIR > 2] <- 2
# HLA_KIR[HLA_KIR < -2] <- -2


# gsva_mat_hallmark_p1_si <- gsva_mat_hallmark_p1_si[, order(match(colnames(gsva_mat_hallmark_p1_si), old_patient_stage$Sample_ID))]

sfd = HeatmapAnnotation(Subtype = old_patient_stage$group1)
row_ha2 = rowAnnotation(foo = anno_text(anno_pvalue$p.signif  , location = 0.5, just = "center"))

# 绘制热图
ComplexHeatmap::Heatmap(gsva_mat_hallmark_p1_si,
                        col = colorRampPalette(c("navy","white","firebrick3"))(100),
                        column_split  = factor(old_patient_stage$group1, levels = c("Cluster1","Cluster2","Cluster3")),
                        name = " ",
                        row_names_side = "left",
                        right_annotation = row_ha2,
                        top_annotation = sfd,
                        show_column_names = F,
                        cluster_rows = F,
                        cluster_columns = F)
pheatmap::pheatmap(gsva_mat_hallmark_p1_si,
                   scale = "row",
                   col = colorRampPalette(c("navy","white","firebrick3"))(200),
                   show_colnames = F,
                   show_rownames = T,
                   cluster_rows = T,
                   cluster_cols = F,
                   #gaps_col = c(275),
                   annotation_col = annotation_col)
write.table(gsva_mat_hallmark_p1_si, file = "/home/cmx/senes_splicing/outdata/CB_data/Fig3HI.txt", sep = "\t",quote = FALSE)

####看Cluster3激活的通路####
# Cluster3
# 设置或导入分组
zz<-clinical_OS[,1:2]
zz$cluster3 <- ifelse(zz$group1 == "Cluster3","Cluster3","others")
group <- zz$cluster3
gsva_mat_hallmark <- gsva_mat_hallmark[,zz$sample]
table(colnames(gsva_mat_hallmark) == zz$sample)

##构造如下数据框
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(gsva_mat_hallmark)
design
##规定哪一组数据与哪一组数据比较
deg = function(exprSet,design,contrast.matrix){
  fit <- lmFit(exprSet,design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  tempOutput = topTable(fit2, coef=1, n=Inf)
  nrDEG = na.omit(tempOutput) 
  head(nrDEG)
  return(nrDEG)
}
compare <- makeContrasts(Cluster3 - others, levels=design)  
gsva_cluster3 = deg(gsva_mat_hallmark,design,compare)
gsva_cluster3 <- gsva_cluster3[order(gsva_cluster3$logFC,decreasing = T),]

#挑选激活的前几个通路
pathway <- c(rownames(gsva_cluster3[1:45,]))

names(pathway) <- c(rep("Cluster3_UP",45))

#求均值
gsva_mat_hallmark <- gsva_mat_hallmark[,zz$sample]
GSVA_data <- gsva_mat_hallmark
GSVA_data <- t(as.data.frame(GSVA_data))
table(rownames(GSVA_data) == zz$sample) 
GSVA_data <- as.data.frame(GSVA_data)
GSVA_data$cluster <- zz$group1

data_heatmap <- data.frame(Cluster1 = apply(GSVA_data[which(GSVA_data$cluster == "Cluster1"),][,-51], 2, mean),
                           Cluster2 = apply(GSVA_data[which(GSVA_data$cluster == "Cluster2"),][,-51], 2, mean),
                           Cluster3 = apply(GSVA_data[which(GSVA_data$cluster == "Cluster3"),][,-51], 2, mean))

up <- data_heatmap[pathway,]

annotation_row <- data.frame(Cluster = c(names(pathway)))
rownames(annotation_row) <- rownames(up)

pheatmap::pheatmap(up,scale = "row",cluster_cols = F,cluster_rows = F ,
                   color = colorRampPalette(c("#009FCC", "white", "#FF6666"))(100),
                   border_color = "white",
                   border_width =10,
                   cellwidth = 12,
                   cellheight = 12,
                   annotation_row = annotation_row)
####看Cluster3下调的通路####
# Cluster3
# 设置或导入分组
zz<-clinical_OS[,1:2]
zz$cluster3 <- ifelse(zz$group1 == "Cluster3","Cluster3","others")
group <- zz$cluster3
gsva_mat_hallmark <- gsva_mat_hallmark[,zz$sample]
table(colnames(gsva_mat_hallmark) == zz$sample)

##构造如下数据框
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(gsva_mat_hallmark)
design
##规定哪一组数据与哪一组数据比较
compare <- makeContrasts(Cluster3 - others, levels=design)  

gsva_cluster3 = deg(gsva_mat_hallmark,design,compare)
gsva_cluster3 <- gsva_cluster3[order(gsva_cluster3$logFC,decreasing = T),]

#挑选激活的前几个通路
pathway <- c(rownames(gsva_cluster3[46:50,]))

names(pathway) <- c(rep("Cluster3_DOWN",5))

#求均值
gsva_mat_hallmark <- gsva_mat_hallmark[,zz$sample]
GSVA_data <- gsva_mat_hallmark
GSVA_data <- t(as.data.frame(GSVA_data))
table(rownames(GSVA_data) == zz$sample) 
GSVA_data <- as.data.frame(GSVA_data)
GSVA_data$cluster <- zz$group1

data_heatmap <- data.frame(Cluster1 = apply(GSVA_data[which(GSVA_data$cluster == "Cluster1"),][,-51], 2, mean),
                           Cluster2 = apply(GSVA_data[which(GSVA_data$cluster == "Cluster2"),][,-51], 2, mean),
                           Cluster3 = apply(GSVA_data[which(GSVA_data$cluster == "Cluster3"),][,-51], 2, mean))

up <- data_heatmap[pathway,]

annotation_row <- data.frame(Cluster = c(names(pathway)))
rownames(annotation_row) <- rownames(up)

pheatmap::pheatmap(up,scale = "row",cluster_cols = F,cluster_rows = F ,
                   color = colorRampPalette(c("#009FCC", "white", "#FF6666"))(100),
                   border_color = "white",
                   border_width =10,
                   cellwidth = 12,
                   cellheight = 12,
                   annotation_row = annotation_row)

####11.画个瀑布图,看三个组突变的差别####
# 确定文件路径！
dir.path <- "./indata/maf/gdc_download_20240226_024428.515507/"

# 获取所有maf文件路径
all.maf <- list.files(path = dir.path, pattern = ".gz", 
                      full.names = T, recursive = T)

# 看看前3个
all.maf[1:3]
## [1] "G:/tcga/GDCdata/TCGA-COAD/harmonized/Simple_Nucleotide_Variation/Masked_Somatic_Mutation/007c2ae4-bbd2-42c6-ab67-bf016fbddb51/982004b5-52e1-4a69-97d3-25bdcb77b026.wxs.aliquot_ensemble_masked.maf.gz"
## [2] "G:/tcga/GDCdata/TCGA-COAD/harmonized/Simple_Nucleotide_Variation/Masked_Somatic_Mutation/010f9040-294d-4d14-a2b4-80d7a11625dd/5083b949-1bf3-4bc2-bf4f-f668f8a13792.wxs.aliquot_ensemble_masked.maf.gz"
## [3] "G:/tcga/GDCdata/TCGA-COAD/harmonized/Simple_Nucleotide_Variation/Masked_Somatic_Mutation/0148fff1-b8af-4bf0-8bcd-de1ff9f750f3/2c16cfe2-bf6d-4a39-af3a-9dfd5ada3e17.wxs.aliquot_ensemble_masked.maf.gz"
maf.list <- lapply(all.maf, data.table::fread, 
                   sep = "\t", 
                   header = T,
                   skip = 7 
)

maf.merge <- do.call(rbind,maf.list)

dim(maf.merge)
## [1] 252664    140
# 读取成功！
library(maftools)
maf1 <- read.maf(maf.merge)
## -Validating
## -Silent variants: 63597 
## -Summarizing
## --Mutiple centers found
## BCM;WUGSC;BCM;BI;BCM;WUGSC--Possible FLAGS among top ten genes:
##   TTN
##   SYNE1
##   MUC16
## -Processing clinical data
## --Missing clinical data
## -Finished in 3.970s elapsed (3.730s cpu)
plotmafSummary(maf = maf1, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
somaticInteractions(maf = maf1, top = 20)

#oncoplot for top ten mutated genes.
oncoplot(maf = maf1, top = 20)

cluster1<-clinical_OS[clinical_OS$group1 == 'Cluster1',]
cluster2<-clinical_OS[clinical_OS$group1 == 'Cluster2',]
cluster3<-clinical_OS[clinical_OS$group1 == 'Cluster3',]


#all
maf.merge <- do.call(rbind,maf.list)

BRCA <- maf.merge
maf.BRCA <- read.maf(BRCA)
maf.BRCA<-maf.BRCA@data
maf.BRCA$Tumor_Sample_Barcode <- substr(maf.BRCA$Tumor_Sample_Barcode,1,12)
maf.BRCA <- maf.BRCA[which(maf.BRCA$Tumor_Sample_Barcode %in% clinical_OS$sample),]
# clindata <- fxxx
# colnames(clindata)[1] <- "Tumor_Sample_Barcode"
mut_all<- read.maf(maf.BRCA)
plotmafSummary(maf = mut_all, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
somaticInteractions(maf = mut_all, top = 20)

oncoplot(maf = mut_all, top = 20)
#group1
maf.merge <- do.call(rbind,maf.list)

BRCA <- maf.merge
maf.BRCA <- read.maf(BRCA)
maf.BRCA<-maf.BRCA@data
maf.BRCA$Tumor_Sample_Barcode <- substr(maf.BRCA$Tumor_Sample_Barcode,1,12)
maf.BRCA <- maf.BRCA[which(maf.BRCA$Tumor_Sample_Barcode %in% cluster1$sample),]
# clindata <- fxxx
# colnames(clindata)[1] <- "Tumor_Sample_Barcode"
mut1 <- read.maf(maf.BRCA)
plotmafSummary(maf = mut1, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
somaticInteractions(maf = mut1, top = 20)

plot1<-oncoplot(maf = mut1, top = 20)


# write.table(maf.BRCA, file = "/home/cmx/senes_splicing/outdata/CB_data/Fig4A.txt",row.names = FALSE, sep="\t")

#group2
maf.merge <- do.call(rbind,maf.list)

BRCA <- maf.merge
maf.BRCA <- read.maf(BRCA)
maf.BRCA<-maf.BRCA@data
maf.BRCA$Tumor_Sample_Barcode <- substr(maf.BRCA$Tumor_Sample_Barcode,1,12)
maf.BRCA <- maf.BRCA[which(maf.BRCA$Tumor_Sample_Barcode %in% cluster2$sample),]
# clindata <- fxxx
# colnames(clindata)[1] <- "Tumor_Sample_Barcode"
mut2 <- read.maf(maf.BRCA)
plotmafSummary(maf = mut2, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
somaticInteractions(maf = mut2, top = 20)

plot2<-oncoplot(maf = mut2, top = 20)


# 使用 oncoplot 绘制这些基因的突变情况
oncoplot(maf = mut2, top = 20)
# write.table(maf.BRCA, file = "/home/cmx/senes_splicing/outdata/CB_data/Fig4B.txt",row.names = FALSE, sep="\t")

#group3
maf.merge <- do.call(rbind,maf.list)

BRCA <- maf.merge
maf.BRCA <- read.maf(BRCA)
maf.BRCA<-maf.BRCA@data
maf.BRCA$Tumor_Sample_Barcode <- substr(maf.BRCA$Tumor_Sample_Barcode,1,12)
maf.BRCA <- maf.BRCA[which(maf.BRCA$Tumor_Sample_Barcode %in% cluster3$sample),]
# clindata <- fxxx
# colnames(clindata)[1] <- "Tumor_Sample_Barcode"
mut3 <- read.maf(maf.BRCA)
plotmafSummary(maf = mut3, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
somaticInteractions(maf = mut3, top = 20)
plot3<-oncoplot(maf = mut3, top = 20)

# 使用 oncoplot 绘制这些基因的突变情况
oncoplot(maf = mut3, top = 20)
# write.table(maf.BRCA, file = "/home/cmx/senes_splicing/outdata/CB_data/Fig4C.txt",row.names = FALSE, sep="\t")
library(venn)

genes = c("TP53","TTN","PIK3CA","MUC16","KMT2C","FLG","SYNE1","USH2A","PKHD1L1","ZFHX4",
          "ANK2","DST","HMCN1","MUC17","PCLO","RELN","ATRX","LRP2","MDN1","RYR2")
# 创建列联表
contingency_table <- matrix(c(52, 43, 39, 161), nrow = 2)
colnames(contingency_table) <- c("TP53突变", "TP53野生型")
rownames(contingency_table) <- c("Cluster3", "非Cluster3")

# 执行Fisher精确检验
fisher_result <- fisher.test(contingency_table)

# 输出结果
print("列联表:")
print(contingency_table)
print("Fisher精确检验结果:")
print(fisher_result)



####12.logTMB三个分组差别####
#计算tmb值
tmb_table_wt_log = tmb(maf = maf1,logScale = F)
#查看tmb值
head(tmb_table_wt_log)
tmb_table_wt_log$Tumor_Sample_Barcode <- substr(tmb_table_wt_log$Tumor_Sample_Barcode,1,12)
tmb_table_wt_log<-as.data.frame(tmb_table_wt_log)

tmb_table_wt_log<-merge(clinical_OS[,c(1,2)],tmb_table_wt_log,by.x="sample", by.y="Tumor_Sample_Barcode" )
tmb_table_wt_log <- tmb_table_wt_log[ !tmb_table_wt_log$sample %in% list,]

ggplot(tmb_table_wt_log,aes(x=factor(group1),y=total_perMB_log,fill=factor(group1), color=factor(group1)))+
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
  ylab("lgTMB")+xlab("Subgroup")
# write.table(tmb_table_wt_log, file = "/home/cmx/senes_splicing/outdata/CB_data/Fig4F.txt",row.names = FALSE, sep="\t")

####三个Cluster的GSEA####
# Cluster1
# 设置或导入分组
zz<-clinical_OS[,1:2]
zz$cluster1 <- ifelse(zz$group1 == "Cluster1","Cluster1","others")
zz <- zz[ !zz$sample %in% list,]
group <- zz$cluster1

gsva_mat_hallmark <- gsva_mat_hallmark[,zz$sample]
table(colnames(gsva_mat_hallmark) == zz$sample)

##构造如下数据框
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
old_patient_count <- old_patient_count[, !colnames(old_patient_count) %in% list]

rownames(design) = colnames(old_patient_count)
design
##规定哪一组数据与哪一组数据比较
compare <- makeContrasts(Cluster1 - others, levels=design)  

gsea_cluster1 = deg(old_patient_count,design,compare)

diffsig <- gsea_cluster1  
diffsig$Gene<-rownames(diffsig)
#diffsig <-  diffsig[diffsig$logFC < 0,]
diffsig <- diffsig[order(diffsig$logFC,decreasing = T), ] #
#diffsig<-diffsig[-c(1:800),]
genelist_input <- diffsig[,c("Gene","logFC")]
genename <- as.character(genelist_input[,1])   #提取第一列基因名
library(ensembldb)
gene_map <-  AnnotationDbi::select(org.Hs.eg.db, keys=genename, keytype="SYMBOL", columns=c("ENTREZID"))  #将SYMBOL格式的ID换成ENTREZ格式的ID。
non_duplicates_idx <- which(duplicated(gene_map$SYMBOL) == FALSE)
gene_map <- gene_map[non_duplicates_idx, ]   #去除重复值
colnames(gene_map)[1]<-"Gene" 
#将ENTREZID与logFC对应起来，并根据logFC的值降序排列，最终生成结果如图所示。
temp<-inner_join(gene_map,genelist_input,by = "Gene")
temp<-temp[,-1]
temp<-na.omit(temp)
temp$logFC<-sort(temp$logFC,decreasing = T) 
geneList = temp[,2]
names(geneList) = as.character(temp[,1])
geneList 
#Go_gseresult <- gseGO(geneList, 'org.Hs.eg.db', keyType = "ENTREZID", ont="BP", pvalueCutoff=1)
library(clusterProfiler)
library(AnnotationDbi)
set.seed(123)  # 设置随机种子以确保结果可重复
geneList <- geneList + rnorm(length(geneList), mean = 0, sd = 1e-5)

# 按值进行降序排序
geneList <- sort(geneList, decreasing = TRUE)
KEGG_gseresult <- gseKEGG(geneList, pvalueCutoff=1)  #使用GSEA进行KEGG富集分析
#保存富集分析结果
kegg_results<-as.data.frame(KEGG_gseresult) 
myenrichplot::gseaplot2(KEGG_gseresult,c('hsa04612',"hsa04210","hsa00061"),title = 'GSEA enrichment analysis of Cluster1 vs others',pvalue_table = T)
myenrichplot::gseaplot2(KEGG_gseresult,c('hsa04612'),title = 'GSEA enrichment analysis of Cluster1 vs others',pvalue_table = T)


# Cluster2
# 设置或导入分组'
zz<-clinical_OS[,1:2]
zz$cluster1 <- ifelse(zz$group1 == "Cluster2","Cluster2","others")
zz <- zz[ !zz$sample %in% list,]
group <- zz$cluster1

gsva_mat_hallmark <- gsva_mat_hallmark[,zz$sample]
table(colnames(gsva_mat_hallmark) == zz$sample)

##构造如下数据框
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
old_patient_count <- old_patient_count[, !colnames(old_patient_count) %in% list]
rownames(design) = colnames(old_patient_count)
design
##规定哪一组数据与哪一组数据比较
compare <- makeContrasts(Cluster2 - others, levels=design)  

gsea_cluster2 = deg(old_patient_count,design,compare)

diffsig <- gsea_cluster2 
diffsig$Gene<-rownames(diffsig)
#diffsig <-  diffsig[diffsig$logFC < 0,]
diffsig <- diffsig[order(diffsig$logFC,decreasing = T), ] #
#diffsig<-diffsig[-c(1:800),]
genelist_input <- diffsig[,c("Gene","logFC")]
genename <- as.character(genelist_input[,1])   #提取第一列基因名
library(ensembldb)
gene_map <-  AnnotationDbi::select(org.Hs.eg.db, keys=genename, keytype="SYMBOL", columns=c("ENTREZID"))  #将SYMBOL格式的ID换成ENTREZ格式的ID。
non_duplicates_idx <- which(duplicated(gene_map$SYMBOL) == FALSE)
gene_map <- gene_map[non_duplicates_idx, ]   #去除重复值
colnames(gene_map)[1]<-"Gene" 
#将ENTREZID与logFC对应起来，并根据logFC的值降序排列，最终生成结果如图所示。
temp<-inner_join(gene_map,genelist_input,by = "Gene")
temp<-temp[,-1]
temp<-na.omit(temp)
temp$logFC<-sort(temp$logFC,decreasing = T) 
geneList = temp[,2]
names(geneList) = as.character(temp[,1])
geneList 
#Go_gseresult <- gseGO(geneList, 'org.Hs.eg.db', keyType = "ENTREZID", ont="BP", pvalueCutoff=1)
library(clusterProfiler)
library(AnnotationDbi)
set.seed(123)  # 设置随机种子以确保结果可重复
geneList <- geneList + rnorm(length(geneList), mean = 0, sd = 1e-5)

# 按值进行降序排序
geneList <- sort(geneList, decreasing = TRUE)
KEGG_gseresult <- gseKEGG(geneList, pvalueCutoff=1,eps =0 )  #使用GSEA进行KEGG富集分析
#保存富集分析结果
kegg_results<-as.data.frame(KEGG_gseresult) 
myenrichplot::gseaplot2(KEGG_gseresult,c('hsa00010',"hsa04066","hsa04612"),title = 'GSEA enrichment analysis of Cluster2 vs others',pvalue_table = T)
myenrichplot::gseaplot2(KEGG_gseresult,c('hsa04066'),title = 'GSEA enrichment analysis of Cluster2 vs others',pvalue_table = T)

write.table(as.data.frame(KEGG_gseresult), file = "/home/cmx/senes_splicing/outdata/CB_data/Fig3J_GSEA2.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Cluster3
# 设置或导入分组'zz<-clinical_OS[,1:2]
zz$cluster1 <- ifelse(zz$group1 == "Cluster3","Cluster3","others")
zz <- zz[ !zz$sample %in% list,]
group <- zz$cluster1

gsva_mat_hallmark <- gsva_mat_hallmark[,zz$sample]
table(colnames(gsva_mat_hallmark) == zz$sample)

##构造如下数据框
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
old_patient_count <- old_patient_count[, !colnames(old_patient_count) %in% list]
rownames(design) = colnames(old_patient_count)
design
##规定哪一组数据与哪一组数据比较
compare <- makeContrasts(Cluster3 - others, levels=design)  

gsea_cluster3 = deg(old_patient_count,design,compare)

diffsig <- gsea_cluster3
diffsig$Gene<-rownames(diffsig)
#diffsig <-  diffsig[diffsig$logFC < 0,]
diffsig <- diffsig[order(diffsig$logFC,decreasing = T), ] #
#diffsig<-diffsig[-c(1:800),]
genelist_input <- diffsig[,c("Gene","logFC")]
genename <- as.character(genelist_input[,1])   #提取第一列基因名
library(ensembldb)
gene_map <-  AnnotationDbi::select(org.Hs.eg.db, keys=genename, keytype="SYMBOL", columns=c("ENTREZID"))  #将SYMBOL格式的ID换成ENTREZ格式的ID。
non_duplicates_idx <- which(duplicated(gene_map$SYMBOL) == FALSE)
gene_map <- gene_map[non_duplicates_idx, ]   #去除重复值
colnames(gene_map)[1]<-"Gene" 
#将ENTREZID与logFC对应起来，并根据logFC的值降序排列，最终生成结果如图所示。
temp<-inner_join(gene_map,genelist_input,by = "Gene")
temp<-temp[,-1]
temp<-na.omit(temp)
temp$logFC<-sort(temp$logFC,decreasing = T) 
geneList = temp[,2]
names(geneList) = as.character(temp[,1])
geneList 
#Go_gseresult <- gseGO(geneList, 'org.Hs.eg.db', keyType = "ENTREZID", ont="BP", pvalueCutoff=1)
library(clusterProfiler)
library(AnnotationDbi)
set.seed(123)  # 设置随机种子以确保结果可重复
geneList <- geneList + rnorm(length(geneList), mean = 0, sd = 1e-5)

# 按值进行降序排序
geneList <- sort(geneList, decreasing = TRUE)
KEGG_gseresult <- gseKEGG(geneList, pvalueCutoff=1,eps =0 )  #使用GSEA进行KEGG富集分析#保存富集分析结果
kegg_results<-as.data.frame(KEGG_gseresult) 
myenrichplot::gseaplot2(KEGG_gseresult,c('hsa04512'),title = 'GSEA enrichment analysis of Cluster3 vs others',pvalue_table = T)
write.table(as.data.frame(KEGG_gseresult), file = "/home/cmx/senes_splicing/outdata/CB_data/Fig3J_GSEA3.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#write.csv2(PSI_data_top,"./outdata/PSI_data_top.csv")
#考虑加进去Reatome
