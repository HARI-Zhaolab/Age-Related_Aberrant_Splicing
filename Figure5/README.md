# Splicing-Regulated Gene Analysis – TCGA BRCA

This script analyzes **splicing perturbations** between tumor and normal samples in **BRCA** using **PSI (Percent Spliced In)** values. It integrates splicing-defined **activation** (acti) and **resistance** (res) events, assesses their correlation with RBP genes, and identifies potential **regulators**. **Spearman correlation** is applied to assess relationships between **splicing events** and **gene expression** data, with **adjusted p-values** for significance.

------

### Inputs

1. **Precomputed objects**
   - `step1_cluster_os.Rdata` → subtype labels (`clinical_OS$group1`)
   - `step2_count_log2.Rdata` → gene expression counts (log2 transformed)
   - `step2_tpm_log2.Rdata` → gene expression TPM (log2 transformed)
   - `old_patient_stage.Rdata` → patient clinical stages
   - `step3_perturbation_level.Rdata` → precomputed perturbation levels (acti, res)
2. **External resources**
   - `papers_Ontology.csv` → RBP gene list (HGNC symbols)

### Outputs

- **Spearman Correlation Results**:
  - `acti_All_cor`: Correlations between **activation** events and RBP gene expressions.
  - `res_All_cor`: Correlations between **resistance** events and RBP gene expressions.
- **Filtered results**:
  - **acti**-only and **res**-only RBPs: Genes exclusively linked to activation or resistance.
  - **Intersected genes**: Genes common to both activation and resistance.
- **Gene lists** (with at least 30 occurrences):
  - `h_acti_30` → RBPs active in splicing events
  - `h_res_30` → RBPs involved in resistant splicing events
- Saved data for further exploration:
  - `h_30.Rdata` (for **acti** and **res** RBP gene lists)

------

### Workflow

Fiugure5.R

1. **Load Data**
2. **Preprocess PSI Data**
3. **Correlation with Gene Expression**
4. **Filter and Process Correlation Results**
5. **Identify Regulators**
   - **Activated** splicing events:
     - Identify RBPs with **positive correlation** (`cor > 0.3`) and **at least 30 occurrences**:
        `h_acti_30`
   - **Resistant** splicing events:
     - Identify RBPs with **negative correlation** (`cor < -0.3`) and **at least 30 occurrences**:
        `h_res_30`
   - **Intersecting genes**:
     - RBPs common to both activated and resistant events:
        `resWact_RBP = intersect(h_res_30$Var1, h_acti_30$Var1)`