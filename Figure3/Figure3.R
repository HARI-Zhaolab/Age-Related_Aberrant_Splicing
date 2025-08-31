## 0) Setup ----------------------------------------------------------------
rm(list = ls()); set.seed(123)
setwd("/home/cmx/senes_splicing/")
.libPaths("/home/cmx/R/x86_64-pc-linux-gnu-library/4.2", include.site = TRUE)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(readr)
  library(ggplot2); library(ggpubr)
  library(ComplexHeatmap); library(pheatmap); library(RColorBrewer)
  library(limma); library(clusterProfiler); library(AnnotationDbi); library(org.Hs.eg.db)
  library(GSVA)
  library(maftools)
})

# ---- Parameters (centralize all paths & aesthetics) ---------------------
OUTDIR               <- "./outdata"
CB_OUTDIR            <- file.path(OUTDIR, "CB_data")
INDIR                <- "./indata"
EXCLUDE_IDS          <- c("TCGA-BH-A18G","TCGA-AC-A23H","TCGA-AN-A046","TCGA-AC-A5XS",
                          "TCGA-A2-A0EV","TCGA-BH-A0HP","TCGA-E9-A1R3","TCGA-C8-A1HM",
                          "TCGA-AC-A2B8","TCGA-LD-A74U","TCGA-EW-A1P5","TCGA-D8-A1J8",
                          "TCGA-B6-A0RQ")
CLUSTER_LEVELS       <- c("Cluster1","Cluster2","Cluster3")
PALETTE_CLUSTER      <- c("#388E8E", "#EE4000","#aa1000")  # violin fills
HEAT_COL             <- colorRampPalette(c("navy","white","firebrick3"))(200)

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CB_OUTDIR, showWarnings = FALSE, recursive = TRUE)

## 1) Load upstream objects ------------------------------------------------
# From your pipeline:
load(file.path(OUTDIR, "step1_cluster_os.Rdata"))   # clinical_OS, group1/cluster labels upstream
load(file.path(OUTDIR, "step2_count_log2.Rdata"))   # old_patient_count, old_patient_count_log
load(file.path(OUTDIR, "step2_tpm_log2.Rdata"))     # old_patient_tpm, old_patient_tpm_log
load(file.path(INDIR,  "GSVA_out.Rdata"))           # ssGSEA_matrix (28 immune cell signatures)
load(file.path(OUTDIR, "old_patient_stage.Rdata"))  # old_patient_stage with Stage/T/N/M/age/group1

# Expect:
# - clinical_OS: columns sample, group1 (Cluster1/2/3), vital_status, days_to_*
# - group1 dataframe upstream merged to clinical data; here we will keep "group1" as label
#   and also create a factor "cluster"
stopifnot(all(c("sample","group1") %in% colnames(clinical_OS)))
clinical_OS <- clinical_OS %>% mutate(cluster = factor(group1, levels = CLUSTER_LEVELS))

## 2) Helper functions ----------------------------------------------------
gg_violin_box <- function(df, x, y, fill_pal = PALETTE_CLUSTER, ylab = "", xlab = "Subgroup") {
  x <- enquo(x); y <- enquo(y)
  ggplot(df, aes(x = !!x, y = !!y, fill = !!x, color = !!x)) +
    geom_violin(trim = FALSE) +
    scale_fill_manual(values = fill_pal) +
    scale_color_manual(values = c("black","black","black")) +
    geom_boxplot(fill = "white", width = 0.2, size = 0.5, outlier.size = 0.2) +
    stat_compare_means(comparisons = list(c("Cluster1","Cluster2"),
                                          c("Cluster2","Cluster3"),
                                          c("Cluster1","Cluster3"))) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_text(colour="black", size=14),
          axis.text.y = element_text(size=14),
          axis.title  = element_text(size=14)) +
    labs(x = xlab, y = ylab)
}

heatmap_with_p <- function(mat_z, col_split, top_anno = NULL, pvals_row = NULL,
                           colmap = HEAT_COL, show_rows = FALSE, show_cols = FALSE, name = " ") {
  ra <- if (!is.null(pvals_row)) rowAnnotation(p = anno_text(pvals_row, location = 0.5, just = "center")) else NULL
  ComplexHeatmap::Heatmap(mat_z,
                          col = colmap,
                          column_split = factor(col_split, levels = CLUSTER_LEVELS),
                          name = name,
                          top_annotation = top_anno,
                          right_annotation = ra,
                          show_column_names = show_cols,
                          show_row_names = show_rows,
                          cluster_rows = FALSE,
                          cluster_columns = FALSE)
}

kruskal_pvals <- function(mat_rows_by_feature, group_vec) {
  # mat_rows_by_feature: rows == features, cols == samples
  df <- as.data.frame(t(mat_rows_by_feature))
  df$Cluster <- factor(group_vec, levels = CLUSTER_LEVELS)
  do.call(rbind, lapply(colnames(df)[colnames(df)!="Cluster"], function(g) {
    tibble(.y. = g, p.signif = compare_means(reformulate("Cluster", g), data = df, method="kruskal.test")$p.signif)
  }))
}

scale_rows <- function(m) { as.matrix(t(scale(t(m)))) }

## 3) ESTIMATE scores (imported from file) --------------------------------
est <- read.table(file.path(OUTDIR, "FPKM_tumor_estimate_score.txt"),
                  sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, row.names = 1)
est <- est[,-1] %>% t() %>% as.data.frame()
rownames(est) <- gsub("\\.", "-", rownames(est))
est <- est %>% mutate(across(where(is.character), as.numeric))
est$sample <- rownames(est)

# Merge cluster labels
est <- est %>%
  inner_join(clinical_OS[,c("sample","group1")], by = "sample") %>%
  mutate(group1 = factor(group1, levels = CLUSTER_LEVELS))

# Violin plots
gg_violin_box(est, group1, ImmuneScore, ylab = "ImmuneScore")
gg_violin_box(est, group1, StromalScore, ylab = "StromalScore")
gg_violin_box(est, group1, ESTIMATEScore, ylab = "ESTIMATEScore")
gg_violin_box(est, group1, TumorPurity,  ylab = "TumorPurity")

## 4) ssGSEA (28 immune cell signatures) ---------------------------------
# ssGSEA_matrix rows = signatures, cols = samples
imm_ssgsea <- as.data.frame(ssGSEA_matrix)
imm_z      <- scale_rows(imm_ssgsea)
# align sample order to old_patient_stage
imm_z      <- imm_z[, order(match(colnames(imm_z), old_patient_stage$Sample_ID))]
# drop excluded samples
keep_idx   <- !colnames(imm_z) %in% EXCLUDE_IDS
imm_z      <- imm_z[, keep_idx]
ops        <- old_patient_stage[match(colnames(imm_z), old_patient_stage$Sample_ID), ]
top_anno   <- HeatmapAnnotation(N = ops$Nstage, M = ops$Mstage, T = ops$Tstage,
                                Stage = ops$Stage, Age = ops$Age, Subgroup = ops$group1)

# p-values (Kruskal–Wallis) across clusters
pv_imm <- kruskal_pvals(imm_ssgsea[, colnames(imm_z)], ops$group1)
heatmap_with_p(imm_z, col_split = ops$group1, top_anno = top_anno, pvals_row = pv_imm$p.signif)

# Export long-format for boxplots if needed
res_long <- t(imm_ssgsea) %>% as.data.frame() %>% tibble::rownames_to_column("sample") %>%
  inner_join(clinical_OS[,c("sample","group1")], by="sample") %>%
  pivot_longer(cols = -c(sample,group1), names_to = "gene", values_to = "value") %>%
  mutate(value40 = value * 80)
write.table(res_long, file = file.path(CB_OUTDIR, "Fig3E_TIL.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

## 5) Immune checkpoint panel --------------------------------------------
genes_icp <- c("CD276","NRP1","CD200","CD40","VSIR","CD70","TNFSF9","CD160","CD40LG",
               "TNFRSF18","TMIGD2","CD80","MTAP","LAIR1","BTLA","CD28","TNFSF4","CD244",
               "CD200R1","LAG3","PDCD1","TNFSF15","HAVCR2","BTNL2","TNFSF18","CD274",
               "PDCD1LG2","ICOSLG","ADORA2A","HHLA2")

icp_mat <- old_patient_count_log[genes_icp, , drop = FALSE] %>% scale_rows()
icp_mat <- icp_mat[, order(match(colnames(icp_mat), old_patient_stage$Sample_ID))]
icp_mat <- icp_mat[, !colnames(icp_mat) %in% EXCLUDE_IDS]
ops2    <- old_patient_stage[match(colnames(icp_mat), old_patient_stage$Sample_ID), ]
pv_icp  <- kruskal_pvals(old_patient_count_log[genes_icp, colnames(icp_mat), drop = FALSE], ops2$group1)
heatmap_with_p(icp_mat, col_split = ops2$group1,
               top_anno = HeatmapAnnotation(Subtype = ops2$group1),
               pvals_row = pv_icp$p.signif)

## 6) HLA panel -----------------------------------------------------------
genes_hla <- c("HLA-E","HLA-DPB2","HLA-C","HLA-DQB1","HLA-DQA1","HLA-DMA","HLA-DRB1","HLA-H",
               "HLA-DRB5","HLA-DOA","HLA-DPB1","HLA-DRA","HLA-DMB","HLA-DPA1")
hla_mat <- old_patient_count_log[genes_hla, , drop = FALSE] %>% scale_rows()
hla_mat <- hla_mat[, order(match(colnames(hla_mat), old_patient_stage$Sample_ID))]
hla_mat <- hla_mat[, !colnames(hla_mat) %in% EXCLUDE_IDS]
ops3    <- old_patient_stage[match(colnames(hla_mat), old_patient_stage$Sample_ID), ]
pv_hla  <- kruskal_pvals(old_patient_count_log[genes_hla, colnames(hla_mat), drop = FALSE], ops3$group1)
heatmap_with_p(hla_mat, col_split = ops3$group1,
               top_anno = HeatmapAnnotation(Subtype = ops3$group1),
               pvals_row = pv_hla$p.signif)

## 7) GSVA (KEGG / GO / Hallmark) ----------------------------------------
# Pre-computed matrices provided:
load(file.path(INDIR, "GSVA_KEGG_GO.Rdata"))     # gsva_mat, gsva_mat_GO
load(file.path(INDIR, "GSVA_hallmark.Rdata"))    # gsva_mat_hallmark

# KEGG
kegg_mat <- gsva_mat[, order(match(colnames(gsva_mat), old_patient_stage$Sample_ID))]
ops4     <- old_patient_stage[match(colnames(kegg_mat), old_patient_stage$Sample_ID), ]
p_kegg   <- kruskal_pvals(kegg_mat, ops4$group1)
kegg_sel <- rownames(kegg_mat) %in% p_kegg$.y.[p_kegg$p.signif != "ns"]
heatmap_with_p(scale_rows(kegg_mat[kegg_sel, , drop = FALSE]),
               col_split = ops4$group1,
               top_anno = HeatmapAnnotation(Subtype = ops4$group1),
               pvals_row = p_kegg$p.signif[p_kegg$.y. %in% rownames(kegg_mat)[kegg_sel]])

# GO
go_mat <- gsva_mat_GO[, order(match(colnames(gsva_mat_GO), old_patient_stage$Sample_ID))]
ops5   <- old_patient_stage[match(colnames(go_mat), old_patient_stage$Sample_ID), ]
p_go   <- kruskal_pvals(go_mat, ops5$group1)
go_sel <- rownames(go_mat) %in% p_go$.y.[p_go$p.signif != "ns"]
heatmap_with_p(scale_rows(go_mat[go_sel, , drop = FALSE]),
               col_split = ops5$group1,
               top_anno = HeatmapAnnotation(Subtype = ops5$group1),
               pvals_row = p_go$p.signif[p_go$.y. %in% rownames(go_mat)[go_sel]])

# Hallmark
hall_mat <- gsva_mat_hallmark[, order(match(colnames(gsva_mat_hallmark), old_patient_stage$Sample_ID))]
ops6     <- old_patient_stage[match(colnames(hall_mat), old_patient_stage$Sample_ID), ]
p_hall   <- kruskal_pvals(hall_mat, ops6$group1)
hall_sel <- rownames(hall_mat) %in% p_hall$.y.[p_hall$p.signif != "ns"]
h_selmat <- scale_rows(hall_mat[hall_sel, , drop = FALSE])
heatmap_with_p(h_selmat, col_split = ops6$group1,
               top_anno = HeatmapAnnotation(Subtype = ops6$group1),
               pvals_row = p_hall$p.signif[p_hall$.y. %in% rownames(hall_mat)[hall_sel]])
write.table(h_selmat, file = file.path(CB_OUTDIR, "Fig3HI.txt"), sep = "\t", quote = FALSE)

## 8) Pathways enriched in Cluster3 (Hallmark, limma contrast) -----------
zz <- clinical_OS[,c("sample","group1")] %>% mutate(cluster3 = ifelse(group1=="Cluster3","Cluster3","others"))
hall_mat2 <- hall_mat[, zz$sample]
design <- model.matrix(~0 + factor(zz$cluster3))
colnames(design) <- levels(factor(zz$cluster3)); rownames(design) <- colnames(hall_mat2)
contrast <- makeContrasts(Cluster3 - others, levels = design)
fit  <- eBayes(contrasts.fit(lmFit(hall_mat2, design), contrast))
tab  <- topTable(fit, coef = 1, n = Inf) %>% na.omit() %>% arrange(desc(logFC))

path_up   <- rownames(tab)[1:45]
path_down <- rownames(tab)[46:50]

# Mean GSVA per cluster
gsva_df <- t(hall_mat2) %>% as.data.frame() %>% mutate(cluster = zz$group1)
mean_mat <- rbind(
  Cluster1 = colMeans(gsva_df[gsva_df$cluster=="Cluster1",-ncol(gsva_df)]),
  Cluster2 = colMeans(gsva_df[gsva_df$cluster=="Cluster2",-ncol(gsva_df)]),
  Cluster3 = colMeans(gsva_df[gsva_df$cluster=="Cluster3",-ncol(gsva_df)])
)
pheatmap(mean_mat[, path_up, drop = FALSE], scale = "row", cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("#009FCC","white","#FF6666"))(100), border_color = "white")
pheatmap(mean_mat[, path_down, drop = FALSE], scale = "row", cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("#009FCC","white","#FF6666"))(100), border_color = "white")

## 9) MAF & TMB -----------------------------------------------------------
# Read MAFs (ensure path points to your MAF folder)
maf_dir  <- "./indata/maf/gdc_download_20240226_024428.515507/"
all_maf  <- list.files(maf_dir, pattern = ".gz$", full.names = TRUE, recursive = TRUE)
maf_list <- lapply(all_maf, data.table::fread, sep = "\t", header = TRUE, skip = 7)
maf_merged <- do.call(rbind, maf_list)

# Full cohort restricted to study samples
maf_obj <- read.maf(maf_merged)
maf_df  <- maf_obj@data
maf_df$Tumor_Sample_Barcode <- substr(maf_df$Tumor_Sample_Barcode, 1, 12)
maf_df  <- maf_df[maf_df$Tumor_Sample_Barcode %in% clinical_OS$sample, ]
mut_all <- read.maf(maf_df)

plotmafSummary(maf = mut_all, rmOutlier = TRUE, addStat = "median", dashboard = TRUE)
somaticInteractions(maf = mut_all, top = 20)
oncoplot(maf = mut_all, top = 20)

# TMB (log scale already produced by maftools if requested)
tmb_tbl <- tmb(maf = maf_obj, logScale = FALSE) %>% as.data.frame()
tmb_tbl$Tumor_Sample_Barcode <- substr(tmb_tbl$Tumor_Sample_Barcode,1,12)
tmb_tbl <- clinical_OS[,c("sample","group1")] %>%
  inner_join(tmb_tbl, by = c("sample" = "Tumor_Sample_Barcode")) %>%
  filter(!sample %in% EXCLUDE_IDS) %>%
  mutate(group1 = factor(group1, levels = CLUSTER_LEVELS))

gg_violin_box(tmb_tbl, group1, total_perMB_log, ylab = "lgTMB")

# Example Fisher's test for TP53 (Cluster3 vs others)
# contingency: rows = Cluster3 vs non-Cluster3; cols = TP53-mut vs WT
# Replace counts below with actual cross-tab from your MAF subset
# fisher.test(matrix(c(a,b,c,d), nrow=2))

## 10) Cluster-wise GSEA (KEGG) on counts --------------------------------
deg_limma <- function(exprSet, design, contrast) {
  fit <- eBayes(contrasts.fit(lmFit(exprSet, design), contrast))
  na.omit(topTable(fit, coef = 1, n = Inf))
}
gsea_one_cluster <- function(expr_counts, target_cluster, exclude_ids = EXCLUDE_IDS) {
  zz <- clinical_OS[,c("sample","group1")] %>%
    mutate(tag = ifelse(group1 == target_cluster, target_cluster, "others")) %>%
    filter(!sample %in% exclude_ids)
  expr_counts <- expr_counts[, zz$sample, drop = FALSE]
  design <- model.matrix(~0 + factor(zz$tag))
  colnames(design) <- levels(factor(zz$tag))
  rownames(design) <- colnames(expr_counts)
  contrast <- makeContrasts(contrasts = paste0(target_cluster,"-others"), levels = design)
  de <- deg_limma(expr_counts, design, contrast) %>% arrange(desc(logFC))
  de$Gene <- rownames(de)
  glist <- de[, c("Gene","logFC")]
  # SYMBOL -> ENTREZ
  gn <- glist$Gene
  map <- AnnotationDbi::select(org.Hs.eg.db, keys = gn, keytype = "SYMBOL", columns = "ENTREZID") %>%
    filter(!duplicated(SYMBOL))
  mat <- inner_join(map %>% rename(Gene = SYMBOL), glist, by = "Gene") %>% select(ENTREZID, logFC) %>% na.omit()
  x <- mat$logFC + rnorm(nrow(mat), 0, 1e-5)     # jitter for ties
  names(x) <- mat$ENTREZID
  x <- sort(x, decreasing = TRUE)
  gseKEGG(x, pvalueCutoff = 1, eps = 0)
}

# Run GSEA for each cluster (writes text files compatible with your CB_data)
kegg_c1 <- gsea_one_cluster(old_patient_count, "Cluster1")
write.table(as.data.frame(kegg_c1), file = file.path(CB_OUTDIR, "Fig3J_GSEA1.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
kegg_c2 <- gsea_one_cluster(old_patient_count, "Cluster2")
write.table(as.data.frame(kegg_c2), file = file.path(CB_OUTDIR, "Fig3J_GSEA2.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
kegg_c3 <- gsea_one_cluster(old_patient_count, "Cluster3")
write.table(as.data.frame(kegg_c3), file = file.path(CB_OUTDIR, "Fig3J_GSEA3.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ======================================================================
# End of script
# ======================================================================
