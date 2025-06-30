setwd('/home/cmx/senes_splicing/')
rm(list=ls())
library(maftools)
.libPaths("/home/cmx/R/x86_64-pc-linux-gnu-library/4.2",include.site = TRUE)
load("./outdata/step1_cluster_os.Rdata")
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