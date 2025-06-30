setwd('/home/cmx/senes_splicing/')
rm(list=ls())
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(harmony)
library(tidyverse)
library(RCurl)
library(cowplot)
library(sctransform)
library(tools)
####加载数据，就五个样本####
# 加载必要的库
#Visualize QC metrics as a violin plot
setwd("/home/cmx/senes_splicing/indata/scRNA_seq65/")
lists <- list.files(pattern = ".rds")


file <- list.files()
data <- list()
for(i in 1:length(file)){
  sc_scdata <- readRDS(file[i])
  data[[i]] <- sc_scdata
}
scdata <- merge(data[[1]],data[2:length(data)])
options(future.globals.maxSize = 400 * 1024^3)
scdata <- readRDS("./12.RDS")
# 执行SCTransform
scdata<- SCTransform(
  scdata,
  vars.to.regress = c("percent.mt", "nCount_RNA"), 
  verbose = FALSE)
# 运行PCA
scdata <- RunPCA(scdata, assay = "SCT", npcs = 50)

# 选择主成分（通过ElbowPlot判断）
ElbowPlot(scdata, ndims = 50) 
# Harmony整合
scdata <- RunHarmony(scdata, group.by.vars = "orig.ident")
# 构建邻域图并聚类 (使用Harmony embeddings)
scdata <- FindNeighbors(scdata, reduction = "harmony", dims = 1:20)%>% 
  FindClusters(resolution = 0.3)

scdata <- BuildClusterTree(scdata, dims = 1:20)  # 使用Harmony降维结果
resolutions <- c(0.1,0.2, 0.3,0.4,0.5, 0.6,0.7, 0.8,0.9, 1.0,1.1, 1.2)  # 设置分辨率范围
scdata <- FindClusters(scdata, resolution = resolutions, algorithm = 1)
clustree(scdata, prefix = "SCT_snn_res.")

# UMAP可视化 (使用Harmony embeddings)
scdata <- RunUMAP(scdata, reduction = "harmony", dims = 1:20)
# 可视化
DimPlot(scdata, reduction = "umap", group.by = "seurat_clusters",raster=FALSE,label = T) 
DimPlot(scdata, reduction = "umap", group.by = "celltype",raster=FALSE,label = F) 
DimPlot(scdata,raster=FALSE,group.by="class",label = F)

genes <- c("MS4A1","CD79A","CD19", #B cell
           "PECAM1","VWF",#Endothelial cel
           "EPCAM","KRT19","PROM1",#Epithelial cell
           "COL1A1","COL1A2","DCN", #Fibroblast
           "CD68","CSF1R","MRC1",#Myeloid cell
           "IGHG1","MZB1","SDC1",#Plasma cell
           "CD3D","CD3E","CD8A","CD4"#T cell
)
Idents(scdata) <- "seurat_clusters"
DotPlot(scdata, features =genes, group.by = "celltype") + 
  RotatedAxis()
scdata@meta.data$celltype<- case_when(
  scdata@meta.data$seurat_clusters == 14 ~ "B cell",
  scdata@meta.data$seurat_clusters == 12 ~ "Endothelial",
  scdata@meta.data$seurat_clusters %in% c(5,13) ~ "Fibroblast",
  scdata@meta.data$seurat_clusters == 3~ "Myeloid",
  scdata@meta.data$seurat_clusters == 8 ~ "Plasma cell",
  scdata@meta.data$seurat_clusters %in% c(2,4)  ~ "T cell",
  scdata@meta.data$seurat_clusters %in% c(0,1,6,7,8,9,10,11,15,16,17,18,19) ~ "Epithelial",
  # 兜底（可选）
)
scdata@meta.data$class<- case_when(
  scdata@meta.data$orig.ident %in% c("CID3948","CID4513","CID4515",
                                     "ERMH0056","sc5rJUQ024",
                                     "sc5rJUQ045","sc5rJUQ051","TNMH0114","TNSH0106") ~ "Cluster3",
  scdata@meta.data$orig.ident %in% c("CID4067","CID4290A","ERMH0040","ERMH0064","ERMH0114","HER2PM0337","sc5rJUQ043",
                                     "sc5rJUQ046","sc5rJUQ050","TNB1MH0131", "TNMH0126") ~ "Non-cluster3")

Idents(scdata) <- "celltype"
Cellratio <- as.data.frame(prop.table(table(scdata$celltype, scdata$class), margin = 2))#计算各组样本不同细胞群比例
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))