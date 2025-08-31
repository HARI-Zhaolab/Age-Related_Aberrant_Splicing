setwd('/home/cmx/senes_splicing/')
rm(list = ls())

# Inputs you already have:
load("/home/cmx/senes_splicing/outdata/step1_cluster_os.Rdata")   # clinical_OS (samples)
load("/home/cmx/senes_splicing/outdata/step2_count_log2.Rdata")   # old_patient_count_log etc.
load("/home/cmx/senes_splicing/outdata/old_patient_stage.Rdata")  # <-- needed for group1/Cluster

library(tidyverse)
library(ggpubr)            # compare_means if you want it
library(ComplexHeatmap)
library(circlize)

## If your predictions are in .Rdata as data.frames:
load("./outdata/CCLE_predict.Rdata")     # CCLE_predict
load("./outdata/GDSC_predict.Rdata")     # GDSC_predict
load("./outdata/GTRP_predict.Rdata")     # CTRP_predict  (FYI “GTRP” looks like a typo)
load("./outdata/CCLE_predict_metabric.Rdata")
load("./outdata/GDSC_predict_metabric.Rdata")
load("./outdata/GTRP_predict_metabric.Rdata")

## ========================= Helpers =======================

# Parse a prediction table (like the CSV oncoPredict writes) into a clean matrix:
# - expects a column with sample IDs (often named "X"); remaining columns are drugs
# - returns a numeric matrix: drugs (rows) × samples (cols)
parse_pred_matrix <- function(df, sample_col_guess = c("X", "sample", "Sample_ID")) {
  # choose sample-ID column
  cands <- intersect(sample_col_guess, colnames(df))
  stopifnot(length(cands) >= 1)
  sid_col <- cands[1]

  # move sample column to rownames; drugs left as columns
  df2 <- df %>% relocate(all_of(sid_col))
  rn <- df2[[sid_col]]
  df2 <- df2 %>% select(-all_of(sid_col))
  # coerce to numeric
  df2 <- df2 %>% mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))
  mat <- t(as.matrix(df2))
  colnames(mat) <- rn
  return(mat)
}

# Standardize by row (z-score per drug), safely handling all-NA rows
zscore_rows <- function(mat) {
  t(scale(t(mat)))
}

# Compute Kruskal–Wallis p-values across clusters for each drug
# returns a data.frame with .y. (drug), p, p.adj, p.signif
kruskal_by_row <- function(mat_drug_by_sample, cluster_df, cluster_col = "group1") {
  # Intersect samples
  common <- intersect(colnames(mat_drug_by_sample), cluster_df$Sample_ID)
  stopifnot(length(common) >= 3)  # need > 2 samples
  M <- mat_drug_by_sample[, common, drop = FALSE]
  cl <- cluster_df %>%
    filter(Sample_ID %in% common) %>%
    arrange(match(Sample_ID, common)) %>%
    pull({{ cluster_col }})
  cl <- factor(cl, levels = c("Cluster1","Cluster2","Cluster3"))

  res <- lapply(rownames(M), function(drug) {
    y <- as.numeric(M[drug, ])
    # split y by cluster and run Kruskal; guard for constant vectors
    if (all(is.na(y)) || length(unique(y[!is.na(y)])) < 2) {
      return(data.frame(.y. = drug, p = NA_real_))
    }
    kw <- try(kruskal.test(y ~ cl), silent = TRUE)
    p <- if (inherits(kw, "try-error")) NA_real_ else kw$p.value
    data.frame(.y. = drug, p = p)
  }) %>% bind_rows()

  res$p.adj <- p.adjust(res$p, method = "BH")
  # simple significance stars like ggpubr::compare_means
  cut_pts <- c(1e-4, 1e-3, 1e-2, 0.05, Inf)
  labs    <- c("****","***","**","*","ns")
  res$p.signif <- labs[findInterval(res$p, cut_pts)]
  res
}

# One-stop: make heatmap for a prediction matrix vs clusters.
# - pred_in can be a data.frame with an "X" column or already a matrix
# - filter_ns: keep only drugs with significant p-values (p.signif != "ns")
plot_pred_heatmap <- function(pred_in,
                              cluster_df,
                              cluster_col = "group1",
                              filter_ns = TRUE,
                              title = "") {
  # Convert to matrix drugs×samples
  pred_mat <- if (is.matrix(pred_in)) pred_in else parse_pred_matrix(pred_in)

  # Align to cohort first to avoid “all NA” rows/cols
  common <- intersect(colnames(pred_mat), cluster_df$Sample_ID)
  if (length(common) < 3) stop("Too few overlapping samples between predictions and cluster table.")
  pred_mat <- pred_mat[, common, drop = FALSE]

  # Row z-score (twice is unnecessary; once is enough)
  pred_z <- zscore_rows(pred_mat)

  # Stats
  ptab <- kruskal_by_row(pred_z, cluster_df, cluster_col = cluster_col)
  if (filter_ns) {
    ptab <- ptab %>% filter(!is.na(p), p.signif != "ns")
  } else {
    ptab <- ptab %>% filter(!is.na(p))
  }
  if (nrow(ptab) == 0) stop("No drugs passed the filtering criteria.")

  pred_sig <- pred_z[ptab$.y., , drop = FALSE]

  # Order columns by cluster for pretty facets
  ord <- cluster_df %>%
    filter(Sample_ID %in% colnames(pred_sig)) %>%
    arrange(factor(.data[[cluster_col]], levels = c("Cluster1","Cluster2","Cluster3")))
  pred_sig <- pred_sig[, ord$Sample_ID, drop = FALSE]

  # Annotations
  top_ann <- HeatmapAnnotation(
    Subtype = ord[[cluster_col]],
    col = list(Subtype = c(Cluster1="#4DBBD5", Cluster2="#00A087", Cluster3="#E64B35"))
  )
  row_ann <- rowAnnotation(
    `p.signif` = anno_text(ptab$p.signif, just = "center", gp = gpar(fontsize = 9))
  )

  # Heatmap colors
  col_fun <- colorRamp2(c(-4, 0, 4), c("#4DBBD5FF", "#FFFFFF", "#ED0000FF"))

  # Plot
  ht <- Heatmap(
    pred_sig,
    name = "z",
    col = col_fun,
    show_row_names = TRUE,
    show_column_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    right_annotation = row_ann,
    top_annotation = top_ann,
    column_split = factor(ord[[cluster_col]], levels = c("Cluster1","Cluster2","Cluster3")),
    column_title = title,
    row_title = paste0("Drugs (n=", nrow(pred_sig), ")")
  )
  draw(ht)
  invisible(list(zmat = pred_z, sig_mat = pred_sig, p_table = ptab, heatmap = ht))
}
