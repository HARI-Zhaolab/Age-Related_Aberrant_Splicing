# Raw RNA-seq Data Availability

Ganzhou cohort：All raw RNA seq data generated from this study can be accessed in the NGDC database under accession number HRA006220 and HRA009553.

FUSCC cohort：The accession number for all the data reported in this paper is NODE: OEP000155. All data can be viewed in The National Omics Data Encyclopedia (NODE) (http://www.biosino.org/node) by pasting the accession (OEP000155) into the text search box or through the URL: http://www.biosino.org/node/project/detail/OEP000155.
![External Data Usage Diagram.png](./External_Data_Usage_Diagram.png)
## **Expression Analysis**: 

- **`1.run_fastp.sh`**
  Preprocessing of raw FASTQ files using **fastp** for quality control, adapter trimming, and filtering low-quality reads.

- **`2.align.sh`**
  Alignment of clean reads to the reference genome using a short-read aligner (e.g., STAR, HISAT2, or BWA depending on configuration).

- **`3.samTobam.sh`**
  Conversion of SAM files to BAM format, followed by sorting and indexing (typically using **samtools**).

- **`4.quantification.sh`**
  Quantification of gene or transcript expression levels (e.g., using **featureCounts**, **HTSeq**, or **Salmon**).

- **`5.merge_result.sh`**
  Merging individual sample quantification results into a unified expression matrix for downstream analysis.

  ## **Splicing Analysis**: 

- rMATS，MISO,SUPPA2 splicing analysis.
