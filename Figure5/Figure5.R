## ======================= Environment =======================
setwd('/home/cmx/senes_splicing/')
rm(list = ls())

## ======================= Libraries =========================
library(estimate)
library(GDCRNATools)
library(TCGAbiolinks)
library(ComplexHeatmap)
library(GSVA)
library(Rtsne)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyverse)   # includes readr, tibble, etc.
library(Hmisc)       # rcorr
library(survival)
library(survminer)

## ======================= Data loads ========================
# Assumed objects:
# - clinical_OS (must have $sample)
# - PSI_data_top (AS_event × Sample matrix)
# - old_patient_count_log (gene × Sample matrix, log or count)
# - step3_perturbation_level.Rdata provides: PSI_data_cancer (AS_event × Sample),
#   PSI_data_normal (data.frame with column PSI_NOR for matched AS_event rows)
load("./outdata/step1_cluster_os.Rdata")
load("./outdata/step2_count_log2.Rdata")
load("./outdata/step2_tpm_log2.Rdata")
load("./outdata/old_patient_stage.Rdata")
load("./outdata/step3_perturbation_level.Rdata")

## ======================= Helpers ===========================
# Keep only intersecting columns/rows and preserve order
align_cols <- function(mat, target_ids) {
  keep <- intersect(colnames(mat), target_ids)
  mat[, keep, drop = FALSE]
}
align_rows <- function(df, target_ids) {
  keep <- intersect(rownames(df), target_ids)
  df[keep, , drop = FALSE]
}
to_numeric_df <- function(df) {
  df[] <- lapply(df, function(x) as.numeric(as.character(x)))
  df
}

# Extract ONLY the X–Y (RBP–AS) block from rcorr(cbind(X,Y))
xy_rcorr <- function(X, Y, method = "spearman") {
  Z <- cbind(X, Y)
  out <- Hmisc::rcorr(as.matrix(Z), type = method)
  p <- ncol(X); q <- ncol(Y)
  r_xy <- out$r[seq_len(p), p + seq_len(q), drop = FALSE]
  p_xy <- out$P[seq_len(p), p + seq_len(q), drop = FALSE]
  list(r = r_xy, P = p_xy)
}

# Tidy the r/P matrices into long form
melt_xy <- function(rmat, pmat, xnames, ynames) {
  stopifnot(identical(dim(rmat), dim(pmat)))
  df <- as.data.frame(as.table(rmat))
  names(df) <- c("RBP_idx", "AS_idx", "cor")
  df$p <- as.vector(pmat)
  df$RBP      <- factor(df$RBP_idx, levels = rownames(rmat), labels = xnames)
  df$AS_event <- factor(df$AS_idx,  levels = colnames(rmat), labels = ynames)
  df[, c("RBP","AS_event","cor","p")]
}

## ======================= PSI preprocessing =================
# Restrict to cohort samples available in PSI_data_top
psi_samples <- clinical_OS$sample
PSI_data_group <- PSI_data_top[, psi_samples[psi_samples %in% colnames(PSI_data_top)], drop = FALSE]
PSI_data_group <- to_numeric_df(PSI_data_group)

# Align cancer PSI to cohort and to the event universe present in PSI_data_group
PSI_data_cancer <- to_numeric_df(PSI_data_cancer)
PSI_data_cancer <- align_cols(PSI_data_cancer, psi_samples)
AS_ids <- rownames(PSI_data_group)
PSI_data_cancer <- align_rows(PSI_data_cancer, AS_ids)

# Make sure normal PSI_NOR matches the same AS events (rows)
PSI_data_normal <- PSI_data_normal[rownames(PSI_data_cancer), , drop = FALSE]

## ======================= Perturbation scoring ==============
# Tumor–normal perturbation per AS event × sample
BRCA_perturbation <- PSI_data_cancer - PSI_data_normal$PSI_NOR
BRCA_perturbation <- to_numeric_df(BRCA_perturbation)

# Keep cohort sample order
BRCA_perturbation <- align_cols(BRCA_perturbation, clinical_OS$sample)

# Event-level activity score (mean perturbation across samples; dynamic denominator)
BRCA_perturbation$acti <- rowSums(BRCA_perturbation, na.rm = TRUE) / ncol(BRCA_perturbation)
BRCA_perturbation$acti_level <- ifelse(BRCA_perturbation$acti > 0, "act", "res")
BRCA_perturbation$event <- rownames(BRCA_perturbation)

# Partition events
AS_act_events <- BRCA_perturbation$event[BRCA_perturbation$acti_level == "act"]
AS_res_events <- BRCA_perturbation$event[BRCA_perturbation$acti_level == "res"]

# Build tumor PSI matrices for each set (AS_event × Sample)
AS_acti_matrix <- align_rows(PSI_data_cancer, AS_act_events)
AS_res_matrix  <- align_rows(PSI_data_cancer, AS_res_events)

## ======================= RBP expression ====================
# Load RBP list; use unique HGNC symbols
Ontology <- readr::read_csv("/home/cmx/senes_splicing/indata/papers_Ontology.csv")
genes <- unique(Ontology$HGNC_symbol)

# old_patient_count_log: gene × Sample
TCGA_KIRC1 <- old_patient_count_log[genes, , drop = FALSE]
TCGA_KIRC1 <- na.omit(TCGA_KIRC1)
TCGA_KIRC1 <- to_numeric_df(as.data.frame(TCGA_KIRC1))

# Build sample-aligned matrices: Samples × Features
common_samples_act <- Reduce(intersect, list(colnames(TCGA_KIRC1), colnames(AS_acti_matrix), clinical_OS$sample))
common_samples_res <- Reduce(intersect, list(colnames(TCGA_KIRC1), colnames(AS_res_matrix),  clinical_OS$sample))

X_act <- t(TCGA_KIRC1[, common_samples_act, drop = FALSE])       # Samples × RBPs
Y_act <- t(AS_acti_matrix[, common_samples_act, drop = FALSE])   # Samples × AS events
X_res <- t(TCGA_KIRC1[, common_samples_res, drop = FALSE])
Y_res <- t(AS_res_matrix[,  common_samples_res, drop = FALSE])

# Drop all-NA columns (required by rcorr)
X_act <- X_act[, colSums(!is.na(X_act)) > 0, drop = FALSE]
Y_act <- Y_act[, colSums(!is.na(Y_act)) > 0, drop = FALSE]
X_res <- X_res[, colSums(!is.na(X_res)) > 0, drop = FALSE]
Y_res <- Y_res[, colSums(!is.na(Y_res)) > 0, drop = FALSE]

## ======================= Correlation (RBP–AS) ==============
# Spearman correlations, extracting ONLY RBP ↔ AS block
corr_act <- xy_rcorr(X_act, Y_act, method = "spearman")
corr_res <- xy_rcorr(X_res, Y_res, method = "spearman")

acti_df <- melt_xy(corr_act$r, corr_act$P, colnames(X_act), colnames(Y_act))
res_df  <- melt_xy(corr_res$r, corr_res$P, colnames(X_res), colnames(Y_res))

# Keep valid AS suffixes (adjust if your naming differs)
suffix_pat <- "-AA$|-AD$|-AP$|-AT$|-ES$|-ME$|-RI$"
acti_df <- acti_df[grepl(suffix_pat, acti_df$AS_event), ]
res_df  <- res_df [grepl(suffix_pat,  res_df$AS_event), ]

# Multiple-testing correction (BH FDR)
acti_df$fdr <- p.adjust(acti_df$p, method = "BH")
res_df$fdr  <- p.adjust(res_df$p,  method = "BH")

# Define activators (positive) and repressors (negative)
acti_hits <- subset(acti_df, fdr < 0.10 & cor >  0.30)
res_hits  <- subset(res_df,  fdr < 0.10 & cor < -0.30)

# Count significant events per RBP (threshold ≥ 30 events)
h_acti_30 <- acti_hits |>
  dplyr::count(RBP, name = "Freq") |>
  dplyr::filter(Freq >= 30)

h_res_30 <- res_hits |>
  dplyr::count(RBP, name = "Freq") |>
  dplyr::filter(Freq >= 30)

# Set operations for convenience
res_only    <- setdiff(h_res_30$RBP, h_acti_30$RBP)
acti_only   <- setdiff(h_acti_30$RBP, h_res_30$RBP)
res_act_RBP <- union(h_res_30$RBP, h_acti_30$RBP)
resWact_RBP <- intersect(h_res_30$RBP, h_acti_30$RBP)

## ======================= Save outputs ======================
# Keep names compatible with your downstream expectations
h_res_30_out  <- h_res_30
h_acti_30_out <- h_acti_30
save(h_res_30_out, h_acti_30_out, file = "./outdata/h_30.Rdata")

# (Optional) also save the full long tables for auditing:
# write.csv(acti_df, "./outdata/RBP_AS_correlations_activating_full.csv", row.names = FALSE)
# write.csv(res_df,  "./outdata/RBP_AS_correlations_repressing_full.csv", row.names = FALSE)
