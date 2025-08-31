setwd('/home/cmx/senes_splicing/')

## Safer .libPaths: prepend your user lib without unknown arguments
.libPaths(c("/home/cmx/R/x86_64-pc-linux-gnu-library/4.2", .libPaths()))

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyverse)
  library(DMwR2)   # for knnImputation (optional)
})

## ---------- 1) Load raw PSI table ----------
# Expect columns:
# 1:3 (event descriptors), 4:10 (meta), 11:... (samples incl. -Norm)
PSI_raw <- read.table(
  file = "/home/cmx/My_Project/indata/TCGA_cancer/BRCA/1_raw_download_data4step1/PSI_download_BRCA.txt",
  header = TRUE, check.names = FALSE, quote = ""
)

## ---------- 2) Filter events ----------
PSI_filt <- PSI_raw %>%
  # keep high-variability events
  dplyr::filter(psi_range >= 0.75, std_psi >= 0.01)

## Identify sample columns programmatically (avoid hard-coded 11:1217)
is_sample_col <- grepl("^(TCGA-|.+-Norm$)", colnames(PSI_filt))
# normalize underscores in sample IDs to hyphens for all sample columns
colnames(PSI_filt)[is_sample_col] <- gsub("_", "-", colnames(PSI_filt)[is_sample_col])

## ---------- 3) (Optional) Impute missing values on sample block ----------
# knnImputation requires numeric; coerce only the sample block
sample_block <- PSI_filt[, is_sample_col, drop = FALSE]
# safely coerce to numeric (non-numeric → NA)
sample_block[] <- lapply(sample_block, function(x) suppressWarnings(as.numeric(x)))

# If you already imputed elsewhere, skip this and read your file.
# Otherwise you can run KNN here (commented because it can be slow):
# sample_block <- knnImputation(sample_block, k = 10)

# If you prefer to reuse the imputed CSV you created earlier:
# imputed <- read.csv("./outdata/PSI_data_filled.csv", header = TRUE)
# imputed <- imputed[, -1, drop = FALSE]
# colnames(imputed) <- gsub("\\.", "-", colnames(imputed))
# Make sure column sets match
# stopifnot(identical(colnames(imputed), colnames(sample_block)))
# sample_block <- imputed

## ---------- 4) Rebuild the full table with imputed sample values ----------
PSI_filt[, is_sample_col] <- sample_block

## Create a compact event ID from the first three columns, trimming whitespace
event_id <- apply(PSI_filt[, 1:3, drop = FALSE], 1, function(x) gsub("\\s+", "", paste(x, collapse = "-")))
rownames(PSI_filt) <- event_id

## ---------- 5) Build event × sample matrix ----------
PSI_mat <- PSI_filt[, is_sample_col, drop = TRUE]
# ensure matrix + numeric
PSI_mat <- as.data.frame(PSI_mat)
PSI_mat[] <- lapply(PSI_mat, function(x) suppressWarnings(as.numeric(x)))
PSI_mat <- as.matrix(PSI_mat)

## ---------- 6) Split into normal vs tumor by column name ----------
is_norm <- grepl("-Norm$", colnames(PSI_mat))
PSI_mat_norm   <- PSI_mat[,  is_norm, drop = FALSE]
PSI_mat_tumor  <- PSI_mat[, !is_norm, drop = FALSE]

## Optional: verify no overlap and both exist
message("Events: ", nrow(PSI_mat),
        " | Normal samples: ", ncol(PSI_mat_norm),
        " | Tumor samples: ",  ncol(PSI_mat_tumor))

## ---------- 7) Clean sample IDs ----------
# For normals, remove the '-Norm' suffix so tumor/normal share the same base Sample_ID
colnames(PSI_mat_norm) <- sub("-Norm$", "", colnames(PSI_mat_norm))

## ---------- 8) Ensure consistent event rows (intersection) ----------
common_events <- intersect(rownames(PSI_mat_norm), rownames(PSI_mat_tumor))
PSI_mat_norm  <- PSI_mat_norm [common_events, , drop = FALSE]
PSI_mat_tumor <- PSI_mat_tumor[common_events, , drop = FALSE]

## ---------- 9) Final outputs in your expected data.frames ----------
# If you need the long data.frame with a Sample_ID column:
PSI_data_normal <- cbind(as.data.frame(t(PSI_mat_norm))) %>%
  tibble::rownames_to_column("Sample_ID")

PSI_data_cancer <- cbind(as.data.frame(t(PSI_mat_tumor))) %>%
  tibble::rownames_to_column("Sample_ID")

# Quick look
# head(PSI_data_normal[, 1:6]); head(PSI_data_cancer[, 1:6])
