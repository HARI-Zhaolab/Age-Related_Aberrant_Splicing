# **Single-Cell RNA Sequencing Analysis of Tumor Microenvironment in BRCA**

This script processes and analyzes single-cell RNA sequencing data from BRCA samples to characterize the tumor microenvironment, with particular focus on immune cell subpopulations and their functional states.

### Inputs:

1. **Single-cell data**: Multiple RDS files containing Seurat objects from BRCA samples
2. **Gene annotation**: Manual correction of gene nomenclature (CCN2 → CTGF)

### Outputs:

- Integrated single-cell dataset with harmonized batch correction
- Cell type annotations for major TME components
- T-cell subpopulation characterization
- Visualization of cellular composition and spatial distribution

### Workflow:

**1. Data Loading and Preprocessing**

- Loads 30 single-cell RDS files from BRCA samples
- Performs gene name standardization (CCN2 to CTGF)
- Merges datasets into a combined Seurat object

**2. Data Integration and Normalization**

- Applies standard Seurat workflow: Normalization, Variable Feature selection, Scaling, PCA
- Uses Harmony algorithm for batch correction across samples
- Performs clustering and UMAP visualization

**3. Cell Type Annotation**

- Identifies major cell types using canonical markers:
  - T-cells (CD3D, CD2, CD3E)
  - Fibroblasts (COL1A1, COL3A1, ACTA2)
  - Myeloid cells (LYZ, TYROBP, MS4A6A)
  - Epithelial cells (KRT19, EPCAM, KRT18)
  - B-cells (IGHM, CD79A, MS4A1)
  - Endothelial cells (CLDN5, FLT1, RAMP2)
- Filters out unknown cell types

**4. T-cell Subpopulation Analysis**

- Subsets and re-analyzes T-cells separately
- Identifies 11 distinct T-cell states:
  - CD4+ subsets: Naive, Activated, Tcm, Tem, Tfh, Treg
  - CD8+ subsets: Naive, Activated, Cytotoxic, Exhausted, Trm
- Uses specific marker genes for each subpopulation:
  - Naive T-cells: SELL, EEF1G, LDHB
  - Activated T-cells: JUNB, FOS, CD40LG
  - Cytotoxic T-cells: NKG7, CST7, PRF1
  - Exhausted T-cells: CXCL13, LAG3, PRKD2

**5. Visualization**

- Creates UMAP plots showing:
  - Sample origin (orig.ident)
  - Cell type composition
  - T-cell subpopulations
- Uses specialized visualization tools (scRNAtoolVis) for enhanced presentation
- Pseudobulk.png：A conceptual explanation of the method for Pseudobulk.
- scRNA-seq_Patient.png：Clinical information of elderly breast cancer patients in single cell data.