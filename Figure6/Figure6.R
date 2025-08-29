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
####加载数据，就五个样本####
# 加载必要的库
library(tools)
library(Rcpp)

save.image(file = "my_environment.RData")
load("my_environment.RData")
#Visualize QC metrics as a violin plot
setwd("/home/cmx/senes_splicing/indata/scRNA_seq65/")
lists <- list.files(pattern = ".rds")

list <- c()
for (i in seq(1,length(lists),by=1)) {
  list[[i]] <- readRDS(lists[i])
}

library(Seurat)

# 假设list是一个包含78个Seurat对象的列表
# list <- list(seurat_obj1, seurat_obj2, ..., seurat_obj78)

# 遍历每个Seurat对象并修改基因名称
for (i in seq_along(list)) {
  # 获取当前Seurat对象
  seurat_obj <- list[[i]]
  
  # 获取基因名称
  genes <- rownames(seurat_obj[["RNA"]]@counts)
  
  # 检查并替换基因名称
  if ("CCN2" %in% genes) {
    # 替换基因名称的索引
    gene_idx <- which(genes == "CCN2")
    
    # 修改counts矩阵的行名
    rownames(seurat_obj[["RNA"]]@counts)[gene_idx] <- "CTGF"
    
    # 修改data矩阵的行名（如果存在）
    if ("data" %in% slotNames(seurat_obj[["RNA"]]) &&
        nrow(seurat_obj[["RNA"]]@data) == length(genes)) {
      rownames(seurat_obj[["RNA"]]@data)[gene_idx] <- "CTGF"
    }
    
    # 修改scale.data矩阵的行名（如果存在）
    if ("scale.data" %in% slotNames(seurat_obj[["RNA"]]) &&
        nrow(seurat_obj[["RNA"]]@scale.data) == length(genes)) {
      rownames(seurat_obj[["RNA"]]@scale.data)[gene_idx] <- "CTGF"
    }
  }
  
  # 更新列表中的Seurat对象
  list[[i]] <- seurat_obj
}

setwd('/home/cmx/senes_splicing/')
scRNA_harmony <- merge(list[[1]], y = list[2:30])
setwd("/home/cmx/senes_splicing/indata/scRNA_seq65/")

####整合数据####

scRNA_harmony <- NormalizeData(scRNA_harmony) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
# saveRDS(scRNA_harmony,"./scRNA_harmony.rds")
# group.by.vars最重要，就是样本之间的整合
# system.time({scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident")})
# 
# scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:15) %>% FindClusters(resolution =0.5)
# 
# scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:15)
# 在别的服务器跑了，直接读结果
scRNA_harmony<-readRDS("/home/cmx/senes_splicing/indata/scRNA_seq65/scRNA_harmony2.rds")
immune.combined<-scRNA_harmony
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 13999
## Number of edges: 589767
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.9094
## Number of communities: 15
## Elapsed time: 4 seconds
# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "orig.ident",raster=FALSE)
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE,raster=FALSE)
p1 + p2
####找marker基因####
#================================================================================
# Allmarkers <- FindAllMarkers(immune.combined, logfc.threshold = 0.3, min.pct = 0.3, only.pos = T)
# write.csv(Allmarkers, file = './outdata/Allmarkers.csv')
#================================================================================
#Manual annotation, reference to published articles
markers <- c("CD3D",'CD2',"CD3E", # T-cell CD2、CD3D、CD3E、CD3G
             "COL1A1","COL3A1","ACTA2",# Fibroblasts
             'LYZ','TYROBP','MS4A6A' ,  # Myeloid
             'KRT19','EPCAM','KRT18', # Epithelial
             'IGHM','CD79A','MS4A1', # B-cell MS4A1、CD19、CD79B、MS4A1
             'CLDN5','FLT1','RAMP2' # Endothelial
)
DotPlot(immune.combined, features = markers, col.min = 0)+coord_flip()

(n=length(unique(immune.combined@meta.data$seurat_clusters)))
celltype=data.frame(ClusterID=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,23,24),
                    celltype='unkown')

celltype[celltype$ClusterID %in% c(0,2,15),2]='T-cell'
celltype[celltype$ClusterID %in% c(4,6,11,16,18,24),2]='Fibroblasts'  # 35待定为Fibroblasts, 根据resolution图，35和30来源于同一群
celltype[celltype$ClusterID %in% c(5,17),2]='Myeloid' 
celltype[celltype$ClusterID %in% c(1,7,8,12,13,14),2]='Epithelial'
celltype[celltype$ClusterID %in% c(9,10),2]='B_cell'
celltype[celltype$ClusterID %in% c(3,18,23),2]='Endothelial'

immune.combined@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  immune.combined@meta.data[which(immune.combined@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(immune.combined@meta.data$celltype)
phe=immune.combined@meta.data

immune.combined<-subset(immune.combined, subset = celltype != 'unkown')
DimPlot(immune.combined, reduction = 'umap', group.by = 'celltype', label = T,raster=FALSE)
DimPlot(immune.combined, reduction = 'umap', label = T,raster=FALSE)
library(scRNAtoolVis)
####美化图片####
clusterCornerAxes(object=immune.combined,
                  reduction='umap',
                  noSplit=T,
                  cornerTextSize=3.5,
                  themebg='bwCorner',
                  addCircle=F,
                  cellLabel = T,
                  #cicAlpha=0.2,
                  nbin=200,
                  clusterCol='celltype')
clusterCornerAxes(object=immune.combined,
                  reduction='umap',
                  noSplit=T,
                  cornerTextSize=3.5,
                  themebg='bwCorner',
                  addCircle=F,
                  cellLabel = T,
                  #cicAlpha=0.2,
                  nbin=200,
                  clusterCol='celltype')
library(scales)
library(ggsci)

####拿免疫细胞继续分析####
bt_cell = immune.combined[,immune.combined@meta.data$celltype %in% c("T-cell")]
sce=bt_cell
sce <- NormalizeData(sce) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
# saveRDS(sce,"./sce.rds")
# group.by.vars最重要，就是样本之间的整合
# system.time({sce <- RunHarmony(sce, group.by.vars = "orig.ident")})
# 

# 在别的服务器跑了，直接读结果
sce<-readRDS("/home/cmx/senes_splicing/indata/scRNA_seq65/sce2.rds")
sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:15) %>% FindClusters(resolution =1)

sce <- RunUMAP(sce, reduction = "harmony", dims = 1:15)
# Look at cluster IDs of the first 5 cells
head(Idents(sce), 5)
table(sce$seurat_clusters) 
DimPlot(sce, reduction = 'umap')
DimPlot(sce, reduction = 'umap',group.by = "orig.ident")

# "Mast cell"=c("CD63"),
# "Neutrophil"=c("FCGR3A","ITGAM"),
# "cDendritic cell"=c("FCER1A","CST3"),
# "pDendritic cell"=c("GZMB"),
# "Monocyte"="LYZ",
# "B cell"=c("MS4A1","CD79A"),
# "Plasma cell"=c("MZB1","IGKC","JCHAIN"),
# "Proliferative signal"=c("MKI67","TOP2A","STMN1"),
# "NK cell"=c("GNLY","NKG7","KLRD1"),
# "T cell"=c("CD3D","CD3E")


markers <- c("SELL","EEF1G","LDHB",#CD4_001_Naive:10
             "JUNB","FOS","CD40LG",#CD4_002_Activated:
             "LMNA","ANXA1","CCR7",#CD4_003_Tcm:, 
             'CCR6','LTB', "TRAC",#CD4_004_Tem:
             'ICOS','CTLA4','BATF', "MAF", "NFATC1" ,"IRF4",#CD4_005_Tfh:
             "FOXP3",'TIGIT','SAT1',#CD4_006_Treg:
             'CD8A','YBX3','CCL5',# CD8_007_Naive:
             'CD69','FOSB','HSPA1B',#CD8_008_Acivated:,,
             'NKG7','CST7','PRF1',#CD8_009_Cytotoxic:
             'CXCL13',"LAG3","PRKD2",#P2RY8CD8_010_ Exhausted:CD3+，CD8+，PD1, TIM3, 1B11, LAG3, BLIMP
             'KLRC1','KLRD1','IFNG'#CD8_011_Trm:
)

DotPlot(sce, features = markers, col.min = 0)+coord_flip()
DimPlot(sce, reduction = 'umap',group.by = "orig.ident")

(n=length(unique(sce@meta.data$seurat_clusters)))
celltype=data.frame(ClusterID=0:(n-1),
                    celltype='unknow')
celltype[celltype$ClusterID %in% c(3),2]='CD4_001_Naive'
celltype[celltype$ClusterID %in% c(0,15),2]='CD4_002_Activated'  
celltype[celltype$ClusterID %in% c(1),2]='CD4_003_Tcm' 
celltype[celltype$ClusterID %in% c(2),2]='CD4_004_Tem'
celltype[celltype$ClusterID %in% c(14),2]='CD4_005_Tfh'#10
celltype[celltype$ClusterID %in% c(12),2]='CD4_006_Treg'#10,1
celltype[celltype$ClusterID %in% c(5,6),2]='CD8_007_Naive'
celltype[celltype$ClusterID %in% c(4),2]='CD8_008_Acivated'
celltype[celltype$ClusterID %in% c(11),2]='CD8_009_Cytotoxic'
celltype[celltype$ClusterID %in% c(7),2]='CD8_010_Exhausted'#11可能
celltype[celltype$ClusterID %in% c(8,9),2]='CD8_011_Trm'
# celltype[celltype$ClusterID %in% c(8,9,15),2]='unknow'
sce@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce@meta.data$celltype)
phe=sce@meta.data

sce<-subset(sce, subset = celltype != 'unknow')

DimPlot(sce, reduction = 'umap', group.by = 'celltype', label = F)
DimPlot(sce, reduction = 'umap',  label = T)

library("scales")
library(ggsci)
mycol = hue_pal()(17)

color = c(pal_d3("category20")(20))
#pal_d3("category20b")(20),
#pal_d3("category20c")(20),
#pal_d3("category10")(10))

DimPlot(sce, reduction = "umap", group.by = "celltype",
        cols = color,pt.size =0.8,label = F,label.box = F)

# DimPlot(sce, pt.size =0.8, cols = color,reduction = 'umap',label = T)

# umap = sce@reductions$umap@cell.embeddings %>%  #坐标信息
#   as.data.frame() %>% 
#   cbind(cell_type = sce@meta.data$celltype) # 注释后的label信息 ，改为cell_type
# 
# head(umap)