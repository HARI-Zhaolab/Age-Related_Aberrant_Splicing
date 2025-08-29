# CNV Landscape & GISTIC Plots by Splicing-Defined Subtypes – TCGA BRCA

This script prepares **Masked Copy Number Segment** inputs per subtype and generates **GISTIC-style genome summaries** (G-score & Frequency) for **Cluster1/2/3**. It reads TCGA-BRCA CNV segments, filters primary tumor aliquots (**01A**), subsets samples by `clinical_OS$group1`, and visualizes **amplifications vs. deletions** across the genome using both **maftools** (chromosome ideograms) and **ggplot2** (continuous genome axis).

------

### Inputs

1. **Precomputed objects**
   - `step1_cluster_os.Rdata` → provides `clinical_OS` with `sample` and `group1` (Cluster1/2/3)
   - `BRCA_CNV_download.rda` (from `TCGAbiolinks::GDCprepare`) → masked CNV segments object
2. **External resources**
   - SNP6 marker meta: `snp6.na35.remap.hg38.subset.txt`
   - GISTIC outputs per cluster :
      `all_lesions.conf_90.txt`, `amp_genes.conf_90.txt`, `del_genes.conf_90.txt`, `scores.gistic`

------

### Outputs

- **Per-cluster segment files** (primary tumor 01A only), ready for GISTIC or downstream tools:
  - `MaskedCopyNumberSegment1.txt` (Cluster1)
  - `MaskedCopyNumberSegment2.txt` (Cluster2)
  - `MaskedCopyNumberSegment3.txt` (Cluster3)
- **Marker file**: `.marker_file.txt`
- **GISTIC genome plots** via `maftools::gisticChromPlot()` (per cluster; ref.build = hg38)
- **Continuous-genome plots** (per cluster; two flavors):
  - **G-score**: Amp positive / Del negative; y-range ≈ `[-0.3, 0.3]`
  - **Frequency**: scaled by ×70; Amp positive / Del negative; y-range tuned per cluster (e.g., ~`[-50, 60]`)
- Ready-to-publish figures with **chromosome boundaries**, **labels**, **dashed separators**, and **legend** for Amp/Del.

------

### Workflow

Figure4.R

1. **Load & prepare**
2. **Per-cluster segment export (for GISTIC)**
   - Add `cluster` tag by `substr(Sample,1,12)` → TCGA barcode at patient level.
   - For each cluster **k ∈ {1,2,3}**:
3. **Filter SNP6 markers (`freqcnv = FALSE`)**
   - Read `snp6.na35.remap.hg38.subset.txt`
   - Promote first row to header, then keep rows with `freqcnv == FALSE`
   - Select columns to `Marker_name, Chromosome, Marker_position`
4. **GISTIC summaries with maftools (per cluster)**
5. **Build genome axis (hg38) for ggplot**
6. **Continuous-genome plots from `scores.gistic` (per cluster)**