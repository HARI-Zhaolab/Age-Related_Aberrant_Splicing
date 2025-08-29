# Immune & Pathway Landscape by Splicing-Defined Subtypes – TCGA BRCA

This script characterizes tumor microenvironment, pathways, mutations, and immunologic features across **three BRCA subgroups (Cluster1/2/3)** previously defined by splicing patterns. It integrates **ESTIMATE**, **ssGSEA/GSVA**, **ComplexHeatmap** visuals, **mutation/TMB** analysis, and **GSEA (KEGG)** comparisons.

------

### Inputs

1. **Precomputed objects**
   - `step1_cluster_os.Rdata` (subtype labels in `clinical_OS` / `group1`)
   - `step2_count_log2.Rdata` (counts and log2)
   - `step2_tpm_log2.Rdata` (TPM and log2)
   - `GSVA_out.Rdata`, `GSVA_KEGG_GO.Rdata`, `GSVA_hallmark.Rdata`, `15SET.Rdata` (ssGSEA/GSVA matrices)
2. **External resources**
   - ESTIMATE outputs: `FPKM_tumor_estimate_score.txt`
   - MAF files: `maf.gz`

------

### Outputs

- **Dimension reduction**: t-SNE and PCA plots of splicing matrices (three clusters)
- **TME scores**: Violin+box plots for **ImmuneScore**, **StromalScore**, **ESTIMATEScore**, **TumorPurity**
- **Immune cell abundance (28 ssGSEA)**: Heatmaps (ComplexHeatmap & pheatmap) with **per-row p-value annotations**
- **Immunomodulator genes** (checkpoints/ligands): Heatmaps with **Kruskal–Wallis p-values**
- **HLA genes**: Heatmaps with **p-values**
- **15 pan-TME/repair/EMT signatures**: Boxplots with statistics
- **GSVA (KEGG/GO/Hallmark)**: Global heatmaps; filtered **significant pathways** with p-value strips
- **Subtype-focused hallmark summaries**: Top **up/down pathways** for **Cluster3**
- **Mutation landscape**: Oncoplots and summaries (all and per-cluster), interaction plots; **TMB** violin plot

------

### Workflow

Figure3.R

1. **Load & prepare**
   - Load subtype labels (`clinical_OS$group1`) and expression matrices (counts, TPM).
   - Align column/sample orders across matrices; remove predefined outlier samples (`list`).
2. **Splicing-defined embedding**
   - t-SNE on `PSI_data_top` → scatter plots colored by **Cluster1/2/3** with 90% & 85% ellipses.
   - PCA on the same set for comparability.
3. **ESTIMATE scores**
   - Parse `FPKM_tumor_estimate_score.txt`, coerce numeric, merge with `group1`.
   - Violin+box plots for **ImmuneScore**, **StromalScore**, **ESTIMATEScore**, **TumorPurity** with multiple comparisons via `stat_compare_means`.
4. **Immune cell composition (28 ssGSEA)**
   - From `GSVA_out.Rdata` (`ssGSEA_matrix`), **row-scale** and heatmap by cluster;
   - Compute **Kruskal–Wallis** p-values per cell type and add as right-side text annotations;
   - Also provide a **long-format** boxplot view and export table (`Fig3E_TIL.txt`).
5. **Immunomodulators & HLA panels**
   - Select curated checkpoint/co-stimulatory/co-inhibitory gene sets; scale by gene;
   - Draw **ComplexHeatmap** split by cluster, annotate with **p-values** (Kruskal–Wallis);
   - Repeat for **HLA** genes.
6. **Fifteen additional signatures**
   - Build 15 gene sets from `3D.csv`, run ssGSEA (precomputed), melt to long format, and draw cluster-stratified boxplots with **Kruskal–Wallis** tests.
7. **GSVA (KEGG / GO / Hallmark)**
   - Use precomputed **GSVA** matrices (`gsva_mat`, `gsva_mat_GO`, `gsva_mat_hallmark`), reorder by samples;
   - Global heatmaps (cluster split).
   - For KEGG & GO: compute **per-pathway Kruskal–Wallis** across clusters, **filter significant**, scale, and redraw with p-value sidebars.
   - For **Hallmark**: same pipeline; export selected matrix (`Fig3HI.txt`).
8. **Cluster-focused hallmark summaries**
   - Contrast design (e.g., `Cluster3 vs others`) on hallmark GSVA with `limma`.
   - Rank by logFC, select **top 45 up** and **5 down** pathways for **Cluster3**;
   - Compute per-cluster means and plot compact heatmaps with pathway class annotations.
9. **Mutation & TMB analysis**
   - Read and merge all MAFs; create **MAF** objects overall and **per cluster**;
   - Produce **summary dashboards**, **co-mutation interactions**, and **oncoplots (top 20)**;
   - Example **Fisher’s exact** test for **TP53** enrichment in **Cluster3**;
   - Compute **TMB** with `maftools::tmb`, log scale, and draw cluster-wise violin+box with stats.
10. **Gene-level GSEA (KEGG) per cluster**