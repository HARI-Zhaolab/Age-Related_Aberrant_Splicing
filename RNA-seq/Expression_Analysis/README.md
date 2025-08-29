# RNA-seq Expression_Analysis

This directory contains a set of shell scripts and configuration files for performing **RNA-seq expression analysis**, from raw sequencing reads to quantified expression results.

------

## Workflow Overview

The pipeline is organized into sequential steps, each implemented as a shell script:

1. **`1.run_fastp.sh`**
    Preprocessing of raw FASTQ files using **fastp** for quality control, adapter trimming, and filtering low-quality reads.
2. **`2.align.sh`**
    Alignment of clean reads to the reference genome using a short-read aligner (e.g., STAR, HISAT2, or BWA depending on configuration).
3. **`3.samTobam.sh`**
    Conversion of SAM files to BAM format, followed by sorting and indexing (typically using **samtools**).
4. **`4.quantification.sh`**
    Quantification of gene or transcript expression levels (e.g., using **featureCounts**, **HTSeq**, or **Salmon**).
5. **`5.merge_result.sh`**
    Merging individual sample quantification results into a unified expression matrix for downstream analysis.

------

## Configuration

- **`environment.yaml`**
   Conda environment file specifying all required software dependencies and versions. This ensures reproducibility of the analysis.

------

## Usage

1. Prepare raw FASTQ sequencing files and place them in the designated input directory.

2. Create and activate the conda environment:

   ```python
   bash
   conda env create -f environment.yaml
   conda activate rna-seq-env
   ```

3. Run the scripts sequentially, for example:

   ```python
   bash
   bash 1.run_fastp.sh
   bash 2.align.sh
   bash 3.samTobam.sh
   bash 4.quantification.sh
   bash 5.merge_result.sh
   ```

4. The final merged expression matrix will be available in the output directory for downstream statistical or machine learning analysis.

------

## Notes

- Input and output directory paths should be set inside each script before execution.
- Scripts can be adapted to specific datasets or aligners depending on project requirements.
- Check log files after each step to ensure successful completion.