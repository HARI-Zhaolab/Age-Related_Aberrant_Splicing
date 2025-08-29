# RNA-seq Alternative Splicing Analysis with rMATS

## Overview

This pipeline processes RNA-seq data from FASTQ files to identify and visualize alternative splicing events using rMATS. The workflow includes quality control, alignment, splicing analysis, and visualization.

## Features

- Quality control and trimming of FASTQ files
- Genome alignment using HISAT2
- SAM to BAM conversion and sorting
- Alternative splicing analysis with rMATS
- Sashimi plot visualization

## Dependencies

Before running this script, ensure the following software is installed:

- **HISAT2** (v2.2.1 or higher)
- **samtools** (v1.10 or higher)
- **rMATS** (v4.1.2)
- **rmats2sashimiplot**
- **Python** (with pysam library)

## Directory Structure

```
project/
├── data/                   # Raw data directory
│   └── *.fastq.gz         # Compressed FASTQ files
├── Index/                 # Reference files directory
│   ├── Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa
│   └── Homo_sapiens.GRCh38.107.chr.gtf
├── output/               # Analysis output directory
├── tmp_output/          # Temporary files directory
└── sashimiplot/         # Visualization output directory
```

## Input File Requirements

1. **Reference files** (located in `Index/` directory):
   - `Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa`: Reference genome sequence
   - `Homo_sapiens.GRCh38.107.chr.gtf`: Genome annotation file
2. **RNA-seq data files**:
   - Compressed FASTQ files (e.g., `1_R1_001.fastq.gz`, `14_R1_001.fastq.gz`)
   - Located in the `data/` directory

## Installation

```
# Create and activate conda environment
conda create --prefix ./rmats_env
conda activate ./rmats_env

# Install dependencies
conda install -c bioconda hisat2
conda install -c bioconda samtools openssl=1.0
conda install -c conda-forge -c bioconda rmats=4.1.2
pip install rmats2sashimiplot pysam
```

## Execution

The analysis is divided into several steps:

### 1. Build HISAT2 Index

```
hisat2-build -p 4 Index/Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa Index/genome
```

### 2. Decompress and Align FASTQ Files

```
# Decompress FASTQ files
gzip -d data/1_R1_001.fastq.gz
gzip -d data/14_R1_001.fastq.gz

# Align to reference genome
hisat2 -p 4 -x Index/genome -U data/1_R1_001.fastq -S Index/1_R1_001.sam
hisat2 -p 4 -x Index/genome -U data/14_R1_001.fastq -S Index/14_R1_001.sam
```

### 3. Convert SAM to BAM and Sort

```
# Convert SAM to BAM
samtools view -bS Index/1_R1_001.sam -o Index/1_R1_001.bam
samtools view -bS Index/14_R1_001.sam -o Index/14_R1_001.bam

# Sort BAM files
samtools sort -o Index/1_R1_001.bam Index/1_R1_001.bam
samtools sort -o Index/14_R1_001.bam Index/14_R1_001.bam

# Index BAM files
samtools index Index/1_R1_001.bam
samtools index Index/14_R1_001.bam
```

### 4. rMATS Analysis

#### Preparation Step

```
# Create sample lists
echo Index/1_R1_001.bam > prep1.txt
echo Index/14_R1_001.bam > prep2.txt

# Run rMATS prep step
python rmats.py --b1 prep1.txt --gtf Index/Homo_sapiens.GRCh38.107.chr.gtf -t paired --readLength 50 --nthread 4 --od output --tmp tmp_output_prep_1 --task prep
python rmats.py --b2 prep2.txt --gtf Index/Homo_sapiens.GRCh38.107.chr.gtf -t paired --readLength 50 --nthread 4 --od output --tmp tmp_output_prep_2 --task prep

# Copy intermediate files
python cp_with_prefix.py prep_1_ tmp_output_post/ tmp_output_prep_1/*.rmats
python cp_with_prefix.py prep_2_ tmp_output_post/ tmp_output_prep_2/*.rmats
```

#### Post Analysis Step

```
# Run rMATS post step
python rmats.py --b1 prep1.txt --b2 prep2.txt --gtf Index/Homo_sapiens.GRCh38.107.chr.gtf -t paired --readLength 50 --nthread 4 --od output --tmp tmp_output_post --task post
```

### 5. Visualization with Sashimi Plots

```
rmats2sashimiplot --b1 Index/1_R1_001.bam --b2 Index/14_R1_001.bam -t SE -e output/SE.MATS.JC.txt --l1 SampleOne --l2 SampleTwo --exon_s 1 --intron_s 5 -o sashimiplot
```