# Alternative Splicing Analysis – BRCA(Age-Stratified)

This script performs a comprehensive analysis of alternative splicing (AS) events across cancer and normal samples, stratified by **age groups** (young vs. old) using TCGA data. It applies **t-SNE visualization, differential splicing event detection, clustering, and event-type characterization** to explore age-specific splicing patterns.

------

### Inputs

step1_cluster_os.Rdata

1. **Splicing Data**:
   - `Sample_old_cancer`, `Sample_old_normol`, `Sample_young_cancer`, `Sample_young_normol`
      (PSI matrices representing alternative splicing events for cancer and normal samples in young/old groups)
2. **Metadata**:
   - Sample group identifiers for stratification (cancer vs. normal, young vs. old).

------

### Outputs

- **t-SNE visualizations**:
  - Clustering patterns of splicing profiles for all four sample groups.
  - Subgroup comparisons (old vs. young, cancer vs. normal).
  - Intersection/union of differentially spliced events.
  - Old-specific AS events projected in both old and young cohorts.
- **Differential Splicing Event Tables**:
  - `outTab_young`: Differential events in young cohort.
  - `outTab_old`: Differential events in old cohort.
  - Annotated with **logFC, p-values, FDR, event type, and direction (Up/Down)**.
- **Heatmaps**: Scaled clustering heatmaps of differential splicing profiles with annotation by sample type.
- **Event Distribution Plots**:
  - Proportional distribution of AS event types (AA, AD, AP, AT, ES, ME, RI) across age groups.
  - Bar plots of total and differential AS event counts (Up vs. Down in cancer).
- **UpSet Plots**: Intersection of differential splicing event types between young and old.

------

### Workflow

Figure1.R

1. **Preprocessing**
   - Combine splicing matrices and label groups.
   - Convert expression values to numeric format.
   - Assign cancer/normal labels for downstream analysis.
2. **Dimensionality Reduction (t-SNE)**
   - Apply **Rtsne** on PSI matrices (perplexity = 20–50).
   - Visualize with **ggplot2**, stratified by group.
3. **Differential Analysis**
   - Perform **Wilcoxon tests** per splicing event between cancer and normal samples (separately for young and old).
   - Compute logFC, median differences, and FDR-adjusted p-values.
   - Select significant events (FDR < 0.05 & |logFC| ≥ log2(1.2)).
4. **Annotation and Grouping**
   - Label events by type (AA, AD, AP, AT, ES, ME, RI).
   - Classify events as **Up** or **Down** in cancer.
5. **Visualization**
   - **t-SNE plots** for all groups and filtered event subsets (intersections/unions).
   - **Heatmaps** of differentially spliced events (scaled, clustered).
   - **Stacked bar plots** of AS event distributions.