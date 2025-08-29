# RNA-seq Alternative Splicing Analysis with MISO

## Overview

This pipeline processes RNA-seq data to identify and quantify alternative splicing events using MISO (Mixture of Isoforms). The workflow includes alignment, indexing, and Bayesian analysis of alternative splicing events.

## Features

- Quality control and alignment with HISAT2
- SAM to BAM conversion and sorting
- MISO annotation preparation
- Alternative splicing quantification with MISO
- Comparison of splicing events between conditions

## Requirements

### Software Dependencies

- HISAT2 (v2.2.1 or higher)
- samtools (v1.10 or higher)
- MISO (v0.5.4 or higher)
- Python with pandas, numpy, and pysam libraries

### Reference Files

- Genome FASTA file
- GTF annotation file
- MISO annotation files (can be generated from GTF)

### Input Data

Paired-end or single-end RNA-seq FASTQ files

## Installation

1. Create and activate conda environment:

```
conda create -n miso_env python=3.7
conda activate miso_env
```

1. Install dependencies:

```
conda install -c bioconda hisat2 samtools
pip install miso pandas numpy pysam
```

1. Install MISO from source (recommended):

```
git clone https://github.com/yardenlab/miso.git
cd miso
python setup.py build
python setup.py install
```

## Directory Structure

```
project/
├── ref/                  # Reference files
│   ├── genome.fa         # Genome sequence
│   └── annotation.gtf    # Gene annotation
├── data/                 # Input FASTQ files
├── index/               # HISAT2 index
├── alignments/          # BAM files
├── miso_annotations/    # MISO annotations
├── miso_results/        # MISO output
└── comparisons/         # Comparison results
```