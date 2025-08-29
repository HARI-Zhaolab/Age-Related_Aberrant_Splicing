# **Drug Sensitivity Prediction Analysis Across Splicing-Defined Subtypes**

This script performs comprehensive drug sensitivity prediction analysis using multiple pharmacological databases to identify subtype-specific therapeutic vulnerabilities in elderly BRCA patients.

### Inputs:

1. **Precomputed data**:
   - `step1_cluster_os.Rdata` - Subtype classifications and clinical data
   - `step2_count_log2.Rdata` - Normalized count data
   - Drug sensitivity prediction results from three databases:
     - GDSC_predict (Genomics of Drug Sensitivity in Cancer)
     - CCLE_predict (Cancer Cell Line Encyclopedia)
     - CTRP_predict (Cancer Therapeutics Response Portal)
   - METABRIC validation dataset predictions

### Outputs:

- Heatmaps displaying drug sensitivity patterns across three splicing-defined subtypes
- Statistical comparisons of drug sensitivity between subtypes
- Validation using external METABRIC dataset

### Workflow:

**1. Data Preparation and Integration**

- Loads precomputed drug sensitivity predictions from GDSC, CCLE, and CTRP databases
- Processes METABRIC dataset for validation (elderly patients ≥65 years)
- Merges predictions with subtype classification data

**2. Statistical Analysis**

- Performs Kruskal-Wallis tests to identify significantly different drugs across subtypes
- Filters results to retain only statistically significant drugs (p < 0.05)
- Standardizes drug sensitivity scores using Z-score normalization

**3. Visualization**

- Creates heatmaps using ComplexHeatmap package:
  - Color scale: blue (-4) to white (0) to red (+4)
  - Column splitting by subtype (Cluster1, Cluster2, Cluster3)
  - Row annotations showing statistical significance
  - Sample ordering based on subtype classification

**4. Database-Specific Analysis**

- **GDSC Analysis**: 746 drugs tested, significant differences identified
- **CCLE Analysis**: 36 drugs tested with subtype-specific patterns
- **CTRP Analysis**: 557 drugs evaluated for differential sensitivity

**5. Validation with METABRIC Dataset**

- Repeats analysis using METABRIC elderly patient cohort
- Confirms subtype-specific drug sensitivity patterns
- Provides external validation of findings