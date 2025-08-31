## =========================================================
## Build PSI matrices by age group (young/old × normal/cancer)
## =========================================================

setwd('/home/cmx/senes_splicing/')
rm(list = ls())

# ---- packages ----
suppressPackageStartupMessages({
  library(dplyr); library(tidyverse)
})

# ---- data ----
load("./clinical_age.Rdata")         # clinical_age with bcr_patient_barcode, age_at_initial_pathologic_diagnosis, ...
# Require PSI_data_normal (rows: events, col 'Sample_ID' + normal samples) and PSI_data_cancer (rows: events, tumor samples)
# If they live in another .Rdata, load it before this block.
stopifnot(exists("PSI_data_normal"), exists("PSI_data_cancer"))

# Coerce numeric matrices (events × samples) and keep Sample_ID column for normals
to_numeric_df <- function(df) { df[] <- lapply(df, function(x) suppressWarnings(as.numeric(as.character(x)))); df }

# Split young / old cohorts
clinical_young <- clinical_age %>% filter(age_at_initial_pathologic_diagnosis <= 45)
clinical_old   <- clinical_age %>% filter(age_at_initial_pathologic_diagnosis >= 65)

# ---------- Helper to merge clinical to PSI by barcode ----------
merge_clinical_psi <- function(clin_df, psi_df, is_normal = FALSE) {
  # psi_df for normals is expected to have a Sample_ID column; tumor usually has colnames as TCGA barcodes.
  df <- psi_df
  if ("Sample_ID" %in% colnames(df)) {
    # Long-ish format (Sample_ID column present) → make sure Sample_ID carries "-Norm" if needed
    out <- clin_df %>% select(bcr_patient_barcode) %>%
      left_join(df, by = c("bcr_patient_barcode" = "Sample_ID"))
    # If Sample_ID in PSI_data_normal was stored without "-Norm" (rare), enforce suffix for consistency after transpose
    out$bcr_patient_barcode <- ifelse(!grepl("-Norm$", out$bcr_patient_barcode),
                                      paste0(out$bcr_patient_barcode, if (is_normal) "-Norm" else ""),
                                      out$bcr_patient_barcode)
    # Drop the clinical columns we don’t need (keep only barcode + PSI columns already joined)
    # (Your original removed cols 2:5; here we keep only PSI columns + barcode)
    out
  } else {
    # Wide matrix style (events × samples): select matching tumor samples directly by column names
    keep <- intersect(colnames(df), clin_df$bcr_patient_barcode)
    bind_cols(
      tibble(bcr_patient_barcode = keep),
      t(df[, keep, drop = FALSE]) %>% as.data.frame() %>% rownames_to_column(NULL)
    )
  }
}

# If your PSI_data_normal is already a matrix with rownames=events and colnames=normal barcodes ending in "-Norm",
# rebuild a data.frame with Sample_ID for a clean merge:
if (!"Sample_ID" %in% colnames(PSI_data_normal)) {
  PSI_data_normal <- PSI_data_normal %>% as.data.frame()
  PSI_data_normal$Sample_ID <- rownames(PSI_data_normal)
  PSI_data_normal <- PSI_data_normal %>% relocate(Sample_ID)
}

# ---- Merge clinical with PSI (young/old × normal/cancer) ----
Sample_young_normol <- merge_clinical_psi(clinical_young, PSI_data_normal,  is_normal = TRUE)
Sample_old_normol   <- merge_clinical_psi(clinical_old,   PSI_data_normal,  is_normal = TRUE)

Sample_young_cancer <- merge_clinical_psi(clinical_young, PSI_data_cancer,  is_normal = FALSE)
Sample_old_cancer   <- merge_clinical_psi(clinical_old,   PSI_data_cancer,  is_normal = FALSE)

# ---- Tidy to matrices: events × samples ----
to_event_by_sample <- function(df_with_barcode_first_col) {
  stopifnot("bcr_patient_barcode" %in% colnames(df_with_barcode_first_col))
  # move barcode to rownames and transpose to events × samples
  m <- df_with_barcode_first_col %>% as.data.frame()
  rownames(m) <- m$bcr_patient_barcode
  m <- m[, setdiff(colnames(m), "bcr_patient_barcode"), drop = FALSE]
  m <- t(m) %>% as.data.frame()
  # coerce numeric
  m <- to_numeric_df(m)
  m
}

M_young_norm <- to_event_by_sample(Sample_young_normol)
M_old_norm   <- to_event_by_sample(Sample_old_normol)
M_young_tum  <- to_event_by_sample(Sample_young_cancer)
M_old_tum    <- to_event_by_sample(Sample_old_cancer)

# ---- Harmonize to same event universe (row intersection) ----
common_events <- Reduce(intersect, list(rownames(M_young_norm), rownames(M_old_norm),
                                        rownames(M_young_tum),  rownames(M_old_tum)))
M_young_norm <- M_young_norm[common_events, , drop = FALSE]
M_old_norm   <- M_old_norm  [common_events, , drop = FALSE]
M_young_tum  <- M_young_tum [common_events, , drop = FALSE]
M_old_tum    <- M_old_tum   [common_events, , drop = FALSE]

# ---- Simple QC summaries ----
message("Events (rows) shared across all matrices: ", length(common_events))
message("Young-Normal samples: ",  ncol(M_young_norm))
message("Young-Tumor  samples: ",  ncol(M_young_tum))
message("Old-Normal   samples: ",  ncol(M_old_norm))
message("Old-Tumor    samples: ",  ncol(M_old_tum))

# Optional: drop columns that are all NA (if any)
drop_all_na_cols <- function(m) { m[, colSums(!is.na(m)) > 0, drop = FALSE] }
M_young_norm <- drop_all_na_cols(M_young_norm)
M_old_norm   <- drop_all_na_cols(M_old_norm)
M_young_tum  <- drop_all_na_cols(M_young_tum)
M_old_tum    <- drop_all_na_cols(M_old_tum)
