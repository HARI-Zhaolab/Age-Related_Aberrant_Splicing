rm(list = ls())
set.seed(211)
setwd("/home/cmx/senes_splicing/")

## ---- Library paths (optional) ----------------------------------------
.libPaths("/home/cmx/R/x86_64-pc-linux-gnu-library/4.2", include.site = TRUE)

## ---- Packages ---------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(survival)
  library(survminer)
  library(NMF)            # for nmf(), predict(), consensusmap()
})

## ---- Data -------------------------------------------------------------
# Required objects (produced upstream):
#   - PSI_data_top: PSI matrix (samples x events or events x samples; see below)
#   - clinical_age: clinical table with bcr_patient_barcode, vital_status,
#                   days_to_last_followup, days_to_death
#   - Optionally other objects loaded below
load("./outdata/step1_cluster_os.Rdata")
load("./outdata/step2_count_log2.Rdata")
load("./outdata/step2_tpm_log2.Rdata")
load("./indata/GSVA_out.Rdata")

## ---- Prepare PSI for NMF ----------------------------------------------
# NMF expects a non-negative matrix (features x samples).
# Assume PSI_data_top is events (rows) x samples (cols).
# If your object is samples x events, transpose it here:
if (nrow(PSI_data_top) < ncol(PSI_data_top)) {
  PSI_for_nmf <- t(PSI_data_top)
} else {
  PSI_for_nmf <- PSI_data_top
}

# Ensure numeric matrix
PSI_for_nmf <- as.matrix(apply(PSI_for_nmf, 2, as.numeric))
stopifnot(all(PSI_for_nmf >= 0, na.rm = TRUE))

## ---- Run / load NMF ---------------------------------------------------
# If you already computed NMF, load it; otherwise compute here.
# K = 3 per manuscript
# NOTE: Brunet + 500 iterations replicate your settings.
# Uncomment to run:
# res1 <- nmf(PSI_for_nmf, rank = 3, method = "brunet",
#             nrun = 1, seed = 211, .options = "v", maxIter = 500)

# If computed earlier:
# load("./outdata/nmf_res_k3.Rdata")  # -> res1

# Safety check
if (!exists("res1")) {
  stop("NMF result 'res1' not found. Run nmf() above or load a saved object.")
}

## ---- QC heatmap for consensus -----------------------------------------
# Visual check of NMF consensus
consensusmap(res1)

## ---- Assign clusters ---------------------------------------------------
group1 <- predict(res1) |> as.data.frame()
colnames(group1) <- "cluster_id"
group1$cluster <- paste0("Cluster", group1$cluster_id)
group1$sample  <- rownames(group1)
group1 <- group1 |> arrange(cluster)
table(group1$cluster)

## ---- Merge with clinical & build OS endpoints --------------------------
# clinical_age must contain: bcr_patient_barcode, vital_status,
# days_to_last_followup, days_to_death
dat <- group1 |>
  inner_join(clinical_age, by = c("sample" = "bcr_patient_barcode"))

dat <- dat |>
  mutate(
    OS.time = ifelse(vital_status == "Alive", days_to_last_followup, days_to_death),
    OS      = ifelse(vital_status == "Alive", 0, 1)
  ) |>
  mutate(OS.time = as.numeric(OS.time)) |>
  # Optional: remove a known outlier as in your original script
  filter(sample != "TCGA-B6-A0RQ") |>
  mutate(
    cluster = factor(cluster, levels = c("Cluster1","Cluster2","Cluster3"))
  )

stopifnot(!any(is.na(dat$OS.time)))

## ---- Survival models ---------------------------------------------------
fit <- survfit(Surv(OS.time, OS) ~ cluster, data = dat)
cox <- coxph(Surv(OS.time, OS) ~ cluster, data = dat)
coxSummary <- summary(cox)

# Pairwise log-rank tests
pair_lr <- survminer::pairwise_survdiff(Surv(OS.time, OS) ~ cluster, data = dat)

print(coxSummary)
print(pair_lr)

## ---- Plot KM -----------------------------------------------------------
km_plot <- ggsurvplot(
  fit = fit, data = dat,
  title = "Elderly patients",
  legend = c(0.85, 0.85),
  legend.title = "Cluster",
  palette = c("#e07854", "#51b1b7", "#CD3333"),
  ggtheme = theme_survminer(),
  risk.table = FALSE,
  axes.offset = TRUE,
  xlab = "Time (days)",
  ylab = "Overall Survival",
  break.y.by = 0.2,
  conf.int = FALSE,
  size = 1.5,
  censor = TRUE,
  pval = TRUE
)

print(km_plot)
