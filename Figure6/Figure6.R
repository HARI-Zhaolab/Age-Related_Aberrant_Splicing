setwd('/home/cmx/senes_splicing/')
rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(magrittr)
  library(harmony)
  library(tidyverse)   # includes readr, tibble, stringr, etc.
  library(cowplot)
  library(sctransform)
  library(tools)
  library(Rcpp)
  library(limma)
})

# ============== 1) Input paths & file discovery ===========
# Keep paths in variables instead of switching wd back-and-forth.
data_dir   <- "/home/cmx/senes_splicing/indata/scRNA_seq65"
out_dir    <- "/home/cmx/senes_splicing/outdata"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

rds_files <- list.files(data_dir, pattern = "\\.rds$", full.names = TRUE)
stopifnot(length(rds_files) > 0)

# ============== 2) Read all Seurat objects safely =========
# NOTE: Do not use `list` as a variable name (it masks base::list()).
obj_list <- lapply(rds_files, readRDS)
obj_ids  <- tools::file_path_sans_ext(basename(rds_files))  # sample IDs from file names

# ============== 3) Robust gene renamer (CCN2 -> CTGF) =====
# This updates counts/data/scale.data and VariableFeatures if present.
rename_gene_in_assay <- function(obj, from = "CCN2", to = "CTGF", assay = "RNA") {
  if (!assay %in% names(obj@assays)) return(obj)
  a <- obj[[assay]]
  feats <- rownames(a)
  if (!from %in% feats) return(obj)

  idx <- which(feats == from)

  # counts
  if (nrow(a@counts) == length(feats)) {
    rownames(a@counts)[idx] <- to
  }

  # data
  if (nrow(a@data) == length(feats)) {
    rownames(a@data)[idx] <- to
  }

  # scale.data (may be empty if not scaled yet)
  if (nrow(a@scale.data) > 0) {
    rn <- rownames(a@scale.data)
    if (length(rn) > 0 && from %in% rn) {
      rownames(a@scale.data)[rownames(a@scale.data) == from] <- to
    }
  }

  # variable features
  vfs <- VariableFeatures(obj, assay = assay)
  if (length(vfs)) {
    vfs[vfs == from] <- to
    VariableFeatures(obj, assay = assay) <- unique(vfs)
  }

  # Update obj
  obj[[assay]] <- a
  return(obj)
}

obj_list <- lapply(obj_list, rename_gene_in_assay, from = "CCN2", to = "CTGF", assay = "RNA")

# ============== 4) Merge all objects ======================
# Use one call to merge with add.cell.ids to keep cell barcodes unique & traceable.
# If you truly only want a subset, slice `obj_list`/`obj_ids` accordingly.
if (length(obj_list) == 1) {
  scRNA <- obj_list[[1]]
} else {
  scRNA <- merge(
    x = obj_list[[1]],
    y = obj_list[-1],
    add.cell.ids = obj_ids,
    project = "scRNA_all"
  )
}

# ============== 5) Basic preprocessing ====================
# If you already computed Harmony/PCA/UMAP on another server and saved an RDS,
# you can skip this block and load it below.
scRNA <- NormalizeData(scRNA) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(verbose = FALSE)

# Optional: run Harmony here if you want to recompute locally.
# scRNA <- RunHarmony(scRNA, group.by.vars = "orig.ident")
# scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = 1:15) %>%
#          FindClusters(resolution = 0.5)
# scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:15)

# If you already have a processed object from another server:
harmony_rds <- file.path(data_dir, "scRNA_harmony2.rds")
if (file.exists(harmony_rds)) {
  scRNA <- readRDS(harmony_rds)
}

# Keep a copy named as in your code
immune.combined <- scRNA

# ============== 6) UMAP visualization =====================
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "orig.ident", raster = FALSE)
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE, raster = FALSE)
print(p1 + p2)

# ============== 7) Quick marker panel (reference-based) ===
markers <- c(
  "CD3D","CD2","CD3E",                     # T-cell
  "COL1A1","COL3A1","ACTA2",               # Fibroblasts
  "LYZ","TYROBP","MS4A6A",                 # Myeloid
  "KRT19","EPCAM","KRT18",                 # Epithelial
  "IGHM","CD79A","MS4A1",                  # B-cell
  "CLDN5","FLT1","RAMP2"                   # Endothelial
)
print(DotPlot(immune.combined, features = markers, col.min = 0) + coord_flip())

# ============== 8) Manual cluster annotation (safe) =======
# Build a mapping only for clusters that exist. Default to "Unknown".
clusters_present <- sort(unique(immune.combined$seurat_clusters))  # factor/char → char
clusters_present <- as.integer(as.character(clusters_present))

annot <- setNames(rep("Unknown", length(clusters_present)), clusters_present)

# Edit the lines below to suit your dataset; only intersecting cluster IDs will be applied.
annot[intersect(c(0, 2, 15), clusters_present)]        <- "T-cell"
annot[intersect(c(4, 6, 11, 16, 18, 24), clusters_present)] <- "Fibroblasts"
annot[intersect(c(5, 17), clusters_present)]           <- "Myeloid"
annot[intersect(c(1, 7, 8, 12, 13, 14), clusters_present)] <- "Epithelial"
annot[intersect(c(9, 10), clusters_present)]           <- "B_cell"
annot[intersect(c(3, 18, 23), clusters_present)]       <- "Endothelial"

# Apply annotation
immune.combined$celltype <- "Unknown"
for (cid in names(annot)) {
  immune.combined$celltype[immune.combined$seurat_clusters == as.integer(cid)] <- annot[[cid]]
}

table(immune.combined$celltype)

# Drop unknowns if desired (typo-proof compared to 'unknow/unkown')
immune.combined <- subset(immune.combined, subset = celltype != "Unknown")

# Plot by celltype
print(DimPlot(immune.combined, reduction = "umap", group.by = "celltype", label = TRUE, raster = FALSE))
print(DimPlot(immune.combined, reduction = "umap", label = TRUE, raster = FALSE))

# Optional prettified plot (scRNAtoolVis)
suppressPackageStartupMessages(library(scRNAtoolVis))
clusterCornerAxes(
  object = immune.combined,
  reduction = 'umap',
  noSplit = TRUE,
  cornerTextSize = 3.5,
  themebg = 'bwCorner',
  addCircle = FALSE,
  cellLabel = TRUE,
  nbin = 200,
  clusterCol = 'celltype'
)

# ============== 9) Subset to T cells and re-analyze =======
tcell <- subset(immune.combined, subset = celltype %in% "T-cell")

# If not processed yet, run a quick pipeline; otherwise load precomputed 'sce2.rds'
tcell <- NormalizeData(tcell) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)

tcell_harmony_rds <- file.path(data_dir, "sce2.rds")
if (file.exists(tcell_harmony_rds)) {
  sce <- readRDS(tcell_harmony_rds)
} else {
  sce <- RunHarmony(tcell, group.by.vars = "orig.ident")
  sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:15) %>% FindClusters(resolution = 1.0)
  sce <- RunUMAP(sce, reduction = "harmony", dims = 1:15)
}

# Check clusters and basic plots
print(head(Idents(sce), 5))
print(table(sce$seurat_clusters))
print(DimPlot(sce, reduction = "umap"))
print(DimPlot(sce, reduction = "umap", group.by = "orig.ident"))

# T cell state markers
t_markers <- c(
  "SELL","EEF1G","LDHB",                # CD4_001_Naive
  "JUNB","FOS","CD40LG",                # CD4_002_Activated
  "LMNA","ANXA1","CCR7",                # CD4_003_Tcm
  "CCR6","LTB","TRAC",                  # CD4_004_Tem
  "ICOS","CTLA4","BATF","MAF","NFATC1","IRF4",  # CD4_005_Tfh
  "FOXP3","TIGIT","SAT1",               # CD4_006_Treg
  "CD8A","YBX3","CCL5",                 # CD8_007_Naive
  "CD69","FOSB","HSPA1B",               # CD8_008_Activated
  "NKG7","CST7","PRF1",                 # CD8_009_Cytotoxic
  "CXCL13","LAG3","PRKD2",              # CD8_010_Exhausted
  "KLRC1","KLRD1","IFNG"                # CD8_011_Trm
)
print(DotPlot(sce, features = t_markers, col.min = 0) + coord_flip())
print(DimPlot(sce, reduction = "umap", group.by = "orig.ident"))

# ============== 10) T cell annotation (safe mapping) ======
clu_present <- sort(unique(sce$seurat_clusters)) |> as.integer()
t_annot <- setNames(rep("Unknown", length(clu_present)), clu_present)

t_annot[intersect(3, clu_present)]           <- "CD4_001_Naive"
t_annot[intersect(c(0, 15), clu_present)]    <- "CD4_002_Activated"
t_annot[intersect(1, clu_present)]           <- "CD4_003_Tcm"
t_annot[intersect(2, clu_present)]           <- "CD4_004_Tem"
t_annot[intersect(14, clu_present)]          <- "CD4_005_Tfh"
t_annot[intersect(12, clu_present)]          <- "CD4_006_Treg"
t_annot[intersect(c(5, 6), clu_present)]     <- "CD8_007_Naive"
t_annot[intersect(4, clu_present)]           <- "CD8_008_Activated"
t_annot[intersect(11, clu_present)]          <- "CD8_009_Cytotoxic"
t_annot[intersect(7, clu_present)]           <- "CD8_010_Exhausted"
t_annot[intersect(c(8, 9), clu_present)]     <- "CD8_011_Trm"

sce$celltype <- "Unknown"
for (cid in names(t_annot)) {
  sce$celltype[sce$seurat_clusters == as.integer(cid)] <- t_annot[[cid]]
}
table(sce$celltype)

# Drop unknown if desired
sce <- subset(sce, subset = celltype != "Unknown")

# Plots
library(scales); library(ggsci)
palette20 <- pal_d3("category20")(20)
print(DimPlot(sce, reduction = "umap", group.by = "celltype", cols = palette20, pt.size = 0.8, label = FALSE, label.box = FALSE))
print(DimPlot(sce, reduction = "umap", label = TRUE))

# ============== 11) Save key outputs ======================
saveRDS(immune.combined, file.path(out_dir, "immune_combined.rds"))
saveRDS(sce,              file.path(out_dir, "tcell_sce.rds"))
# Optional: save whole session at the end (not at the beginning)
# save.image(file = file.path(out_dir, "my_environment_final.RData"))
