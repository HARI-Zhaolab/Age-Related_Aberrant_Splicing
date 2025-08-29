# Splicing-Defined Subtype Embedding – TCGA BRCA (t-SNE)

This snippet embeds samples—previously clustered by splicing features—into a low-dimensional space and visualizes separation among **Cluster1/2/3**.

### Inputs

- `step1_cluster_os.Rdata`
  - Provides `clinical_OS` with sample IDs (`sample`) and subtype labels (`group1`).
- `PSI_data_top` (in workspace)
  - PSI matrix of selected splicing events (rows = events, columns = samples).

### What it produces

- A **t-SNE embedding** (`tsne_out$Y` → `tsne_res`) of samples using `PSI_data_top`.
- Two **ggplot t-SNE scatterplots**:
  - Points colored by **Cluster1/2/3**.
  - **90%** (and an additional **85%**) dashed normal ellipses per cluster.
  - Axes labeled “t-SNE 1/2”; legend hidden in the plotting area.

### Workflow

Figure2.R

1. **Load & libraries**
    Loads subtype labels; sets library paths; imports analysis/plotting packages.
2. **t-SNE embedding**
3. **Visualization**
   - Plot 1: 90% ellipse; custom palette `#FF0000/#1874CD/#cc33CD`.
   - Plot 2: overlays 90% and 85% ellipses; palette `#b93a38/#5c7197/#4f7f58`.