#!/bin/bash

set -e  # Exit on error

# Configuration
THREADS=16
MEMORY="80G"
REF_DIR="Index"
DATA_DIR="data"
OUTPUT_DIR="output"
TMP_DIR="tmp_output"
SASHIMI_DIR="sashimiplot"

# Sample names (modify as needed)
SAMPLES=("1" "14")

# Create directories
mkdir -p ${REF_DIR} ${DATA_DIR} ${OUTPUT_DIR} ${TMP_DIR} ${SASHIMI_DIR}

echo "=== RNA-seq Alternative Splicing Analysis Pipeline ==="
echo "Start time: $(date)"
echo "======================================================"

# Step 1: Build HISAT2 index
echo "Step 1/7: Building HISAT2 index..."
hisat2-build -p ${THREADS} \
    ${REF_DIR}/Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa \
    ${REF_DIR}/genome

# Step 2: Decompress FASTQ files
echo "Step 2/7: Decompressing FASTQ files..."
for sample in "${SAMPLES[@]}"; do
    if [ -f "${DATA_DIR}/${sample}_R1_001.fastq.gz" ]; then
        gunzip -k ${DATA_DIR}/${sample}_R1_001.fastq.gz
    fi
done

# Step 3: Alignment with HISAT2
echo "Step 3/7: Aligning reads with HISAT2..."
for sample in "${SAMPLES[@]}"; do
    if [ -f "${DATA_DIR}/${sample}_R1_001.fastq" ]; then
        hisat2 -p ${THREADS} \
            -x ${REF_DIR}/genome \
            -U ${DATA_DIR}/${sample}_R1_001.fastq \
            -S ${REF_DIR}/${sample}_R1_001.sam
    fi
done

# Step 4: Convert SAM to BAM and sort
echo "Step 4/7: Converting SAM to BAM and sorting..."
for sample in "${SAMPLES[@]}"; do
    if [ -f "${REF_DIR}/${sample}_R1_001.sam" ]; then
        # Convert SAM to BAM
        samtools view -bS ${REF_DIR}/${sample}_R1_001.sam -o ${REF_DIR}/${sample}_R1_001.bam
        
        # Sort BAM
        samtools sort -o ${REF_DIR}/${sample}_R1_001.bam ${REF_DIR}/${sample}_R1_001.bam
        
        # Index BAM
        samtools index ${REF_DIR}/${sample}_R1_001.bam
        
        # Clean up SAM file
        rm ${REF_DIR}/${sample}_R1_001.sam
    fi
done

# Step 5: rMATS preparation
echo "Step 5/7: rMATS preparation step..."
# Create sample lists
rm -f prep1.txt prep2.txt
echo "${REF_DIR}/${SAMPLES[0]}_R1_001.bam" > prep1.txt
echo "${REF_DIR}/${SAMPLES[1]}_R1_001.bam" > prep2.txt

# Run rMATS prep steps
python rmats.py --b1 prep1.txt \
    --gtf ${REF_DIR}/Homo_sapiens.GRCh38.107.chr.gtf \
    -t paired --readLength 50 --nthread ${THREADS} \
    --od ${OUTPUT_DIR} --tmp ${TMP_DIR}_prep_1 --task prep

python rmats.py --b2 prep2.txt \
    --gtf ${REF_DIR}/Homo_sapiens.GRCh38.107.chr.gtf \
    -t paired --readLength 50 --nthread ${THREADS} \
    --od ${OUTPUT_DIR} --tmp ${TMP_DIR}_prep_2 --task prep

# Copy intermediate files
python cp_with_prefix.py prep_1_ ${TMP_DIR}_post/ ${TMP_DIR}_prep_1/*.rmats
python cp_with_prefix.py prep_2_ ${TMP_DIR}_post/ ${TMP_DIR}_prep_2/*.rmats

# Step 6: rMATS post analysis
echo "Step 6/7: rMATS post analysis..."
python rmats.py --b1 prep1.txt --b2 prep2.txt \
    --gtf ${REF_DIR}/Homo_sapiens.GRCh38.107.chr.gtf \
    -t paired --readLength 50 --nthread ${THREADS} \
    --od ${OUTPUT_DIR} --tmp ${TMP_DIR}_post --task post

# Step 7: Generate Sashimi plots
echo "Step 7/7: Generating Sashimi plots..."
rmats2sashimiplot --b1 ${REF_DIR}/${SAMPLES[0]}_R1_001.bam \
    --b2 ${REF_DIR}/${SAMPLES[1]}_R1_001.bam \
    -t SE -e ${OUTPUT_DIR}/SE.MATS.JC.txt \
    --l1 SampleOne --l2 SampleTwo \
    --exon_s 1 --intron_s 5 \
    -o ${SASHIMI_DIR}

# Cleanup temporary files
echo "Cleaning up temporary files..."
rm -f prep1.txt prep2.txt
rm -rf ${TMP_DIR}_prep_1 ${TMP_DIR}_prep_2 ${TMP_DIR}_post

echo "======================================================"
echo "Pipeline completed successfully!"
echo "End time: $(date)"
echo "Output files:"
echo "  - rMATS results: ${OUTPUT_DIR}/"
echo "  - Sashimi plots: ${SASHIMI_DIR}/"
echo "======================================================"
