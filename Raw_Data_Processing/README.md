# Splicing Variation Analysis - TCGA BRCA

This script processes TCGA BRCA data to analyze alternative splicing events, focusing on different age groups (young and old patients). It merges splicing data (PSI values) with clinical data and filters out irrelevant or low-confidence splicing events.

### Inputs:

1. **PSI Data**: PSI_download_BRCA.txt-The raw PSI data for BRCA samples.
2. **Clinical Data**: clinical_age.Rdata-Patient clinical data (age, diagnosis, etc.).

### Outputs:

- Filtered PSI data based on age (young vs. old patients).
- Merged datasets for different groups:
  - Young normal samples
  - Young cancer samples
  - Old normal samples
  - Old cancer samples

### Workflow:

Raw_Data_Processing_Clin.R & Raw_Data_Processing_PSI.R

1. **Load Data**: Reads raw PSI data and clinical data.
2. **Filter PSI Data**: Selects splicing events with PSI ≥ 0.75 and standard deviation ≥ 0.01.
3. **Clinical Data Filtering**: Extracts patient age information and classifies them into young (≤45) and old (≥65) groups.
4. **Merge Data**: Combines clinical data with splicing data to generate subsets for young vs. old, and normal vs. cancer samples.