# RNA-seq Alternative Splicing Analysis with SUPPA2

## Overview

This is an automated Bash script for processing RNA-seq data and performing alternative splicing analysis. The pipeline integrates multiple bioinformatics tools including Salmon, SUPPA2, and Trim Galore to provide a complete workflow from raw sequencing data to alternative splicing event analysis.

## Features

- Transcript quantification (using Salmon)
- Alternative splicing event identification and PSI value calculation (using SUPPA2)
- Data quality control and trimming (using Trim Galore)
- Parallel processing for improved efficiency
- Comprehensive logging

## Dependencies

Before running this script, ensure the following software is installed:

- **Salmon** (v0.14.0 or higher)
- **SUPPA2** (v2.3 or higher)
- **Trim Galore** (v0.6.0 or higher)
- **GNU Parallel**
- **AWK**
- **Python** (required for SUPPA2)

## Directory Structure

Recommended project directory structure:

text


```
project/
├── ref/                    # Reference files directory
│   ├── gencode.v48.transcripts.fa
│   └── gencode.v48.annotation.gtf
├── rawData/               # Raw data directory
│   ├── sample.txt         # Sample list file
│   └── *_RNAseq_R{1,2}.fastq.gz
├── cleanData/             # Output directory for QCed data
├── salmon_output/         # Salmon quantification results directory
└── RNAseq_AS_Analysis.sh  # Main script file
```

## Input File Requirements

1. **Reference files** (located in `ref/` directory):
   - `gencode.v48.transcripts.fa`: GENCODE v48 transcript sequences
   - `gencode.v48.annotation.gtf`: GENCODE v48 annotation file
2. **Sample list file**:
   - `rawData/sample.txt`: Text file containing all sample names, one per line
3. **Sequencing data files**:
   - Naming format: `{sample_name}_RNAseq_R1.fastq.gz` and `{sample_name}_RNAseq_R2.fastq.gz`
   - Located in the `rawData/` directory

## Execution

bash


```
# Submit to SLURM cluster
sbatch RNAseq_AS_Analysis.sh

# Or run directly (not recommended except for testing environments)
./RNAseq_AS_Analysis.sh
```

## Output Files

### 1. Salmon Index Files

- `ref/gencode.v48.transcripts.salmon.index/`: Salmon transcript index

### 2. SUPPA2 Event Files

- `ref/gencode.v48.events.*.ioe`: Various alternative splicing event files
- `ref/gencode.v48.all.events.ioe`: Merged event file

### 3. Quality-controlled Data

- `cleanData/`: Contains trimmed fastq files and QC logs

### 4. Quantification Results

- `salmon_output/`: Contains Salmon quantification results for each sample
- `iso_tpm_formatted.txt`: Merged transcript TPM expression matrix

### 5. PSI Value Results

- `project_events.psi`: PSI value matrix for alternative splicing events across all samples

## Key Parameters

- **Computational resources**: Requests 16 CPU cores and 80GB memory
- **Quality trimming**: Phred quality score threshold 28, read length retention threshold 36bp
- **Salmon parameters**: Uses strand-specific library type (ISF) with GC bias correction

## Important Notes

1. Ensure all dependency tools are properly installed and added to PATH environment variable
2. Reference file versions should match those specified in the script
3. Sample naming must be consistent with names in `sample.txt`
4. The script is designed for SLURM cluster environments but can be run on local servers (requires modification of SBATCH directives)
5. Runtime depends on data volume and computational resources

## Error Handling

The script includes error output redirection, with stderr from all tools saved to respective log files:

- `RNAseq_AS.job.err`: Main script error log
- `../cleanData/{}_trim.log`: Trim Galore QC log
- `../salmon_output/{}_salmon.log`: Salmon quantification log
- `psiPerEvent_log.txt`: SUPPA2 PSI calculation log
