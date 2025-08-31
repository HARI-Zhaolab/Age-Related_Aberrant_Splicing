# ------------------------ setup ------------------------
setwd('/home/cmx/senes_splicing/')
rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(TCGAbiolinks)
  library(maftools)
  library(ggplot2)
  library(ggsci)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

# -------------------- helpers --------------------------
stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path, call. = FALSE)
}

# safer loader: returns the object itself (not its name)
load_object <- function(path) {
  stop_if_missing(path)
  nm <- load(path)
  if (length(nm) != 1L) stop("Expected exactly one object in ", path)
  get(nm, envir = .GlobalEnv)
}

# force Chromosome to 1:22 integers; drop sex/MT etc.
coerce_chr_1_22 <- function(x) {
  x <- as.character(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x <- ifelse(x %in% c("X","Y","M","MT"), NA, x)
  as.integer(x)
}

# TCGA primary tumor filter (robust)
is_primary_01 <- function(barcode) {
  # TCGA barcode: positions 14–15 are sample type
  substr(barcode, 14, 15) == "01"
}

# genome lengths (hg38, autosomes 1–22)
hg38_df <- {
  df <- data.frame(
    chromName   = as.character(GenomeInfoDb::seqnames(BSgenome.Hsapiens.UCSC.hg38)),
    chromlength = as.numeric(GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg38)),
    stringsAsFactors = FALSE
  )
  df <- df[match(as.character(1:22), df$chromName), ] # order 1..22
  df$chromNum <- 1:22
  df$chromlengthCumsum    <- cumsum(df$chromlength)
  df$chormStartPosFrom0   <- c(0, head(df$chromlengthCumsum, -1))
  mid <- diff(c(0, df$chromlengthCumsum)) / 2
  df$chromMidelePosFrom0  <- df$chormStartPosFrom0 + mid
  df$ypos <- rep(c(-0.3, 0.3), length.out = 22)
  df
}

# unified plotting for G-score or frequency
plot_gistic_metric <- function(scores_path, df_genome = hg38_df,
                               metric = c("G.score", "frequency"),
                               amp_del_scale = NULL,
                               y_limits = NULL,
                               title = NULL) {
  stop_if_missing(scores_path)
  metric <- match.arg(metric)
  scores <- read.table(scores_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

  # harmonize chromosome as 1..22 integers
  scores$Chromosome <- coerce_chr_1_22(scores$Chromosome)
  scores <- dplyr::filter(scores, !is.na(Chromosome), Chromosome >= 1L, Chromosome <= 22L)

  # genome-wide positions
  scores$StartPos <- scores$Start + df_genome$chormStartPosFrom0[scores$Chromosome]
  scores$EndPos   <- scores$End   + df_genome$chormStartPosFrom0[scores$Chromosome]

  # make Del negative for signed area plots
  if (metric == "G.score") {
    scores[[metric]][scores$Type == "Del"] <- -scores[[metric]][scores$Type == "Del"]
  } else if (metric == "frequency") {
    scores[[metric]][scores$Type == "Del"] <- -scores[[metric]][scores$Type == "Del"]
    if (!is.null(amp_del_scale)) scores[[metric]] <- scores[[metric]] * amp_del_scale
  }

  p <- ggplot(scores, aes(StartPos, .data[[metric]])) +
    geom_area(aes(group = Type, fill = factor(Type, levels = c("Del","Amp")))) +
    scale_fill_lancet(guide = guide_legend(reverse = TRUE), name = "") +
    geom_vline(data = df_genome, aes(xintercept = chromlengthCumsum), linetype = 2) +
    geom_text(data = dplyr::mutate(df_genome, ypos = if (metric == "frequency") rep(c(-52, 51), length.out = 22) else df_genome$ypos),
              aes(x = chromMidelePosFrom0, y = ypos, label = chromName), size = 5) +
    scale_x_continuous(expand = c(0, -1000), limits = c(0, max(df_genome$chromlengthCumsum)),
                       name = NULL, labels = NULL) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background  = element_rect(fill = "transparent", colour = NA),
      legend.position  = "right",
      axis.text.y      = element_text(color = "black", size = 12),
      axis.title.y     = element_text(color = "black", size = 13)
    ) +
    labs(y = metric, title = title %||% basename(scores_path))

  if (!is.null(y_limits)) p <- p + coord_cartesian(ylim = y_limits)
  p
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# -------------------- data I/O --------------------------
# clinical OS / clusters
stop_if_missing("./outdata/step1_cluster_os.Rdata")
load("./outdata/step1_cluster_os.Rdata")  # provides clinical_OS (assumed)

# CNV segment download object
seg_obj <- load_object("./outdata/BRCA_CNV_download.rda")
# Expect a GDCprepare() result with columns including: Sample, Chromosome, Start, End, Num_Probes, Segment_Mean
cnv <- seg_obj %>%
  as.data.frame() %>%
  dplyr::select(Sample, Chromosome, Start, End, Num_Probes, Segment_Mean)

# primary tumors only (robust to barcode formats)
cnv_primary <- cnv %>% dplyr::filter(is_primary_01(Sample))

# add patient barcode (first 12 chars) to join with clinical_OS
cnv_primary <- cnv_primary %>% mutate(cluster_barcode = substr(Sample, 1, 12))

# ---------------- per-cluster segment export -------------
export_cluster_segments <- function(df, cluster_label, out_path) {
  keep_ids <- clinical_OS %>%
    dplyr::filter(group1 == cluster_label) %>%
    dplyr::pull(sample) %>%
    unique()
  seg <- df %>%
    dplyr::filter(cluster_barcode %in% keep_ids) %>%
    dplyr::select(Sample, Chromosome, Start, End, Num_Probes, Segment_Mean)
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  write.table(seg, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(seg)
}

seg1 <- export_cluster_segments(cnv_primary, "Cluster1", "./outdata/MaskedCopyNumberSegment1.txt")
seg2 <- export_cluster_segments(cnv_primary, "Cluster2", "./outdata/MaskedCopyNumberSegment2.txt")
seg3 <- export_cluster_segments(cnv_primary, "Cluster3", "./outdata/MaskedCopyNumberSegment3.txt")

# ---------------- marker file (freqcnv == FALSE) --------
marker_path <- "./indata/snp6.na35.remap.hg38.subset.txt"
if (file.exists(marker_path)) {
  filemarker <- read.table(marker_path, header = FALSE, sep = "", stringsAsFactors = FALSE, quote = "", comment.char = "")
  # first row contains header
  colnames(filemarker) <- as.character(filemarker[1, ])
  filemarker <- filemarker[-1, ]
  if (!all(c("freqcnv","Marker_name","Chromosome","Marker_position") %in% colnames(filemarker))) {
    warning("Marker file columns not as expected; skipping marker export.")
  } else {
    marker_clean <- filemarker %>%
      dplyr::filter(freqcnv %in% c(FALSE, "FALSE")) %>%  # handle logical or character
      dplyr::select(Marker_name, Chromosome, Marker_position)
    write.table(marker_clean, file = './outdata/marker_file.txt', sep = '\t', quote = FALSE, row.names = FALSE)
  }
}

# ---------------- GISTIC plotting (clusters 1–3) --------
make_cluster_plots <- function(cluster_id, ylims_g = c(-0.3, 0.3), ylims_f = NULL, freq_scale = 70) {
  pth <- file.path("./outdata", sprintf("cluster%d_cnv", cluster_id))
  all_lesions <- file.path(pth, "all_lesions.conf_90.txt")
  amp_genes   <- file.path(pth, "amp_genes.conf_90.txt")
  del_genes   <- file.path(pth, "del_genes.conf_90.txt")
  scores_gis  <- file.path(pth, "scores.gistic")

  # quickly validate presence of key files
  for (f in c(all_lesions, amp_genes, del_genes, scores_gis)) stop_if_missing(f)

  # (Optional) quick maftools genome plot
  gobj <- readGistic(gisticAllLesionsFile = all_lesions,
                     gisticAmpGenesFile  = amp_genes,
                     gisticDelGenesFile  = del_genes,
                     gisticScoresFile    = scores_gis,
                     isTCGA = TRUE)
  gisticChromPlot(gistic = gobj, markBands = "all", ref.build = 'hg38')

  # Custom G.score area plot
  p_g <- plot_gistic_metric(scores_gis, metric = "G.score",
                            y_limits = ylims_g,
                            title = sprintf("Cluster %d — G.score (hg38)", cluster_id))
  print(p_g)

  # Custom frequency area plot (Amp positive, Del negative); scale for visual range
  p_f <- plot_gistic_metric(scores_gis, metric = "frequency",
                            amp_del_scale = freq_scale,
                            y_limits = ylims_f,
                            title = sprintf("Cluster %d — Frequency ×%g (hg38)", cluster_id, freq_scale))
  print(p_f)
}

# tune y-limits for your preferred stagge
