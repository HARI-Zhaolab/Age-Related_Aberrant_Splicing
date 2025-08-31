## ---- Setup ------------------------------------------------------------
rm(list = ls())
setwd("/home/cmx/senes_splicing/")

# Library paths (optional; keep if you need a custom path)
.libPaths("/home/cmx/R/x86_64-pc-linux-gnu-library/4.2", include.site = TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(limma)
  library(Rtsne)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(ComplexHeatmap)
})

set.seed(422)

## ---- Parameters -------------------------------------------------------
PERPLEXITY_DEFAULT <- 20
PERPLEXITY_LARGE   <- 50
LOGFC_THRESHOLD    <- log2(1.2)   
FDR_THRESHOLD      <- 0.05
TSNE_THETA         <- 0.5

## ---- Helper functions -------------------------------------------------

# Make all character columns numeric (needed after transposing data.frames)
to_numeric_df <- function(df) {
  df %>% mutate(across(where(is.character), as.numeric))
}

# t-SNE scatter with a clean ggplot theme
plot_tsne <- function(mat, groups, colors, labels, title = "") {
  tsne_out <- Rtsne(mat, pca = FALSE, check_duplicates = FALSE,
                    perplexity = PERPLEXITY_DEFAULT, theta = TSNE_THETA)
  tsne_df <- as.data.frame(tsne_out$Y)
  colnames(tsne_df) <- c("tSNE1", "tSNE2")
  group_df <- tibble(group = groups)

  ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = group_df$group)) +
    geom_point(size = 1, shape = 19) +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(color = "black", fill = "transparent"),
      legend.title = element_text(size = 14, face = "bold"),
      legend.position = "top",
      legend.text = element_text(size = 8, color = "black", face = "bold")
    ) +
    labs(x = "tSNE-1", y = "tSNE-2", color = "", title = title) +
    scale_color_manual(values = colors, labels = labels)
}

# Two-group differential analysis (tumor vs normal) with t-test as in your code.
# Note: your original comment said "Wilcoxon" but used t.test(). I keep t.test()
# for fidelity; consider wilcox.test() for non-normal PSI.
diff_as_two_group <- function(expr_mat, n_tumor, n_normal,
                              fdr_cut = FDR_THRESHOLD, logfc_cut = LOGFC_THRESHOLD) {

  stopifnot(ncol(expr_mat) == (n_tumor + n_normal))

  grade <- c(rep(1, n_normal), rep(2, n_tumor))  # 1 = normal, 2 = tumor
  data  <- as.matrix(expr_mat)

  out <- lapply(rownames(data), function(g) {
    rt <- cbind(expression = data[g, ], grade = grade) %>% as.data.frame()
    tt <- t.test(expression ~ grade, data = rt)  # original choice
    conMean   <- mean(data[g, 1:n_normal])               # normal
    treatMean <- mean(data[g, (n_normal + 1):ncol(data)])# tumor
    logFC     <- log2(treatMean) - log2(conMean)
    conMed    <- median(data[g, 1:n_normal])
    treatMed  <- median(data[g, (n_normal + 1):ncol(data)])
    diffMed   <- treatMed - conMed
    # keep only if sign(logFC) agrees with sign(diffMed), as in your code
    if (((logFC > 0) & (diffMed > 0)) || ((logFC < 0) & (diffMed < 0))) {
      tibble(
        gene   = g,
        conMean = conMean,
        treaMean = treatMean,
        logFC  = logFC,
        pValue = tt$p.value
      )
    } else {
      NULL
    }
  }) %>% bind_rows()

  out$fdr <- p.adjust(out$pValue, method = "BH")
  out <- out %>%
    mutate(logFC = as.numeric(logFC)) %>%
    filter(fdr < fdr_cut & (logFC > logfc_cut | logFC < -logfc_cut)) %>%
    mutate(
      # your original code negated logFC afterwards
      logFC = -logFC,
      type  = substr(gene, nchar(gene) - 1, nchar(gene)),
      up_down = case_when(
        logFC > 0  ~ "Up",
        logFC < 0  ~ "Down",
        TRUE       ~ "Unchanged"
      )
    )

  return(out)
}

# Heatmap helper (sample annotation by "-Norm" suffix)
heatmap_as <- function(expr_mat_scaled) {
  sample_types <- ifelse(grepl("-Norm$", colnames(expr_mat_scaled)), "Normal", "Tumor")
  annotation <- data.frame(Sample_Type = factor(sample_types, levels = c("Normal", "Tumor")))
  rownames(annotation) <- colnames(expr_mat_scaled)

  color_annotation <- list(Sample_Type = c(Normal = "#4477AA", Tumor = "#CC6677"))
  color_heatmap <- colorRampPalette(brewer.pal(9, "RdBu"))(255)

  pheatmap::pheatmap(
    expr_mat_scaled,
    scale = "row",
    clustering_row = TRUE,
    clustering_col = TRUE,
    distance = "euclidean",
    annotation_col = annotation,
    annotation_colors = color_annotation,
    annotation_names_col = FALSE,
    show_rownames = FALSE,
    show_colnames = FALSE,
    color = color_heatmap,
    na_col = "gray"
  )
}

## ---- Load data --------------------------------------------------------
# Produced in "Raw_Data_Processing"
load("./step1_cluster_os.Rdata")
# Objects used below: clinical_young, clinical_old, PSI_data_normal, PSI_data_cancer

## ---- Build four PSI matrices ------------------------------------------
# Young - normal
young_normal <- merge(clinical_young, PSI_data_normal,
                      by.x = "bcr_patient_barcode", by.y = "Sample_ID")
young_normal$bcr_patient_barcode <- ifelse(
  !grepl("-Norm", young_normal$bcr_patient_barcode),
  paste0(young_normal$bcr_patient_barcode, "-Norm"),
  young_normal$bcr_patient_barcode
)
young_normal <- young_normal[, -c(2:5)] %>% t() %>% as.data.frame()
colnames(young_normal) <- young_normal[1, ]; young_normal <- young_normal[-1,]

# Young - cancer
young_cancer <- merge(clinical_young, PSI_data_cancer,
                      by.x = "bcr_patient_barcode", by.y = "Sample_ID")
young_cancer <- young_cancer[, -c(2:5)] %>% t() %>% as.data.frame()
colnames(young_cancer) <- young_cancer[1, ]; young_cancer <- young_cancer[-1,]

# Old - normal
old_normal <- merge(clinical_old, PSI_data_normal,
                    by.x = "bcr_patient_barcode", by.y = "Sample_ID")
old_normal$bcr_patient_barcode <- ifelse(
  !grepl("-Norm", old_normal$bcr_patient_barcode),
  paste0(old_normal$bcr_patient_barcode, "-Norm"),
  old_normal$bcr_patient_barcode
)
old_normal <- old_normal[, -c(2:5)] %>% t() %>% as.data.frame()
colnames(old_normal) <- old_normal[1, ]; old_normal <- old_normal[-1,]

# Old - cancer
old_cancer <- merge(clinical_old, PSI_data_cancer,
                    by.x = "bcr_patient_barcode", by.y = "Sample_ID")
old_cancer <- old_cancer[, -c(2:5)] %>% t() %>% as.data.frame()
colnames(old_cancer) <- old_cancer[1, ]; old_cancer <- old_cancer[-1,]

## ---- t-SNE: all four groups ------------------------------------------
expr_all <- cbind(old_cancer, old_normal, young_cancer, young_normal) %>%
  t() %>% as.data.frame() %>% to_numeric_df()

groups_all <- rep(c("old_cancer", "old_normal", "young_cancer", "young_normal"),
                  c(342, 31, 186, 31))
cols_all   <- c("#191970", "#9ACD32", "#1E90FF", "#EE6AA7")

plot_tsne(expr_all, groups_all, cols_all,
          labels = c("old_cancer","old_normal","young_cancer","young_normal"),
          title = "t-SNE on all samples")

## ---- t-SNE: old cohort only (optionally drop outliers) ----------------
# Example: drop one normal and one tumor as in your script
old_normal2 <- old_normal[, -5, drop = FALSE]
old_cancer2 <- old_cancer[, -1, drop = FALSE]

expr_old <- cbind(old_cancer2, old_normal2) %>%
  t() %>% as.data.frame() %>% to_numeric_df()
groups_old <- rep(c("old_cancer","old_normal"), c(341, 30))
cols_old   <- c("#191970", "#9ACD32")
plot_tsne(expr_old, groups_old, cols_old,
          labels = c("old_cancer","old_normal"),
          title = "t-SNE (old cohort)")

## ---- t-SNE: young cohort (optionally drop a few) ----------------------
# Example: drop selected columns as in your script
young_cancer2 <- young_cancer[, -38, drop = FALSE]
young_normal2 <- young_normal[, -c(3, 12), drop = FALSE]

expr_young <- cbind(young_cancer2, young_normal2) %>%
  t() %>% as.data.frame() %>% to_numeric_df()
groups_young <- rep(c("young_cancer","young_normal"), c(185, 29))
cols_young   <- c("#1E90FF", "#EE6AA7")
plot_tsne(expr_young, groups_young, cols_young,
          labels = c("young_cancer","young_normal"),
          title = "t-SNE (young cohort)")

## ---- Differential AS: young (normal vs tumor) -------------------------
young_group <- cbind(young_normal2, young_cancer2)
res_young <- diff_as_two_group(
  expr_mat = as.matrix(young_group),
  n_tumor  = 185,
  n_normal = 29
)

## ---- Differential AS: old (normal vs tumor) ---------------------------
old_group <- cbind(old_normal2, old_cancer2)
res_old <- diff_as_two_group(
  expr_mat = as.matrix(old_group),
  n_tumor  = 341,
  n_normal = 30
)

save(res_young, res_old, file = "./outdata/outTab_YwO.Rdata")

## ---- Heatmap (old differential events) --------------------------------
old_sig_genes <- res_old$gene
old_cancer_sig  <- old_cancer[old_sig_genes, ]
old_normal_sig  <- old_normal[old_sig_genes, ]
all_exp_old     <- cbind(old_cancer_sig, old_normal_sig) %>% to_numeric_df()
heatmap_as(all_exp_old)

## ---- AS type proportions (example plot) --------------------------------
as_counts <- tibble(
  Type = c("old","young"),
  AA = c(57, 30), AD = c(89, 51), AP = c(699, 475),
  AT = c(338,232), ES = c(499,392), ME = c(10,9), RI = c(133,61)
)
as_long <- as_counts %>%
  pivot_longer(-Type, names_to = "AS_Type", values_to = "Count")
mycols <- c("#EEB422","#6495ED","#3A5FCD","#FFF68F","#8FBC8F","#228B22","#A52A2A")

ggplot(as_long, aes(Type, Count, fill = AS_Type)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0.01, 0.02))) +
  scale_fill_manual(values = mycols) +
  labs(x = "", y = "Proportion of AS events", fill = "", title = "") +
  theme_classic() +
  theme(axis.title.y = element_text(size = 14),
        legend.text  = element_text(size = 10),
        axis.text.x  = element_text(size = 12, color = "black"),
        axis.text.y  = element_text(size = 12, color = "black"),
        axis.ticks.length = unit(0.3, "cm"))

## ---- Barplot: counts (example) ----------------------------------------
totals <- tibble(Type = c("old","young"), Total = c(1825, 1250))
ggplot(totals, aes(Type, Total, fill = Type)) +
  geom_col(color = "white") +
  geom_text(aes(label = Total), vjust = -0.3, size = 4) +
  scale_fill_manual(values = c("#3A5FCD","#8FBC8F")) +
  labs(x = "", y = "Total Number of AS events", fill = "", title = "") +
  theme_classic()

## ---- Intersection/Union analyses for t-SNE -----------------------------
# Intersect significant events between old and young (tumor sets)
intersect_genes <- intersect(rownames(old_cancer[res_old$gene, , drop = FALSE]),
                             rownames(young_cancer[res_young$gene, , drop = FALSE]))

expr_inter <- cbind(
  old_cancer[intersect_genes, ], old_normal[intersect_genes, ],
  young_cancer[intersect_genes, ], young_normal[intersect_genes, ]
) %>% t() %>% as.data.frame() %>% to_numeric_df()

groups_inter <- rep(c("old_cancer","old_normal","young_cancer","young_normal"),
                    c(341, 30, 185, 29))
plot_tsne(expr_inter, groups_inter, cols_all,
          labels = c("old_cancer","old_normal","young_cancer","young_normal"),
          title = "t-SNE on intersecting differential events")

# Union of significant events (old ∪ young)
union_genes <- union(res_old$gene, res_young$gene)
expr_union <- cbind(
  old_cancer[union_genes, ], old_normal[union_genes, ],
  young_cancer[union_genes, ], young_normal[union_genes, ]
) %>% t() %>% as.data.frame() %>% to_numeric_df()

plot_tsne(expr_union, groups_inter, cols_all,
          labels = c("old_cancer","old_normal","young_cancer","young_normal"),
          title = "t-SNE on union of differential events")

## ---- Old-specific events heatmap --------------------------------------
old_only_genes <- setdiff(res_old$gene, res_young$gene)
all_exp_old_only <- cbind(old_cancer[old_only_genes, ], old_normal[old_only_genes, ])
# row-wise scaling
all_exp_old_only <- as.data.frame(t(scale(t(apply(all_exp_old_only, 2, as.numeric)))))
sample_types <- ifelse(grepl("-Norm$", colnames(all_exp_old_only)), "Normal", "Tumor")

ha <- HeatmapAnnotation(
  Sample_Type = factor(sample_types, levels = c("Tumor", "Normal")),
  col = list(Sample_Type = c(Normal = "#4477AA", Tumor = "#CC6677")),
  annotation_name_side = "left"
)

Heatmap(as.matrix(all_exp_old_only),
        top_annotation = ha,
        column_split = factor(sample_types, levels = c("Tumor","Normal")),
        show_column_names = FALSE,
        show_row_names    = FALSE,
        cluster_rows      = TRUE,
        cluster_columns   = FALSE)
