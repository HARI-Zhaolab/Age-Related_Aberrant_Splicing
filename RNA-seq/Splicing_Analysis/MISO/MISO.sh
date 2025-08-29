#!/bin/bash

set -e  # Exit on error

# Configuration
THREADS=16
MEMORY="80G"
REF_DIR="ref"
DATA_DIR="data"
INDEX_DIR="index"
ALIGN_DIR="alignments"
MISO_ANN_DIR="miso_annotations"
MISO_RES_DIR="miso_results"
COMP_DIR="comparisons"

# Sample information (modify as needed)
SAMPLES=("sample1" "sample2")
CONDITIONS=("condition1" "condition2")  # Group samples by condition

# Create directories
mkdir -p ${REF_DIR} ${DATA_DIR} ${INDEX_DIR} ${ALIGN_DIR} 
mkdir -p ${MISO_ANN_DIR} ${MISO_RES_DIR} ${COMP_DIR}

echo "=== MISO Alternative Splicing Analysis Pipeline ==="
echo "Start time: $(date)"
echo "==================================================="

# Step 1: Build HISAT2 index
echo "Step 1/8: Building HISAT2 index..."
if [ ! -f "${INDEX_DIR}/genome.1.ht2" ]; then
    hisat2-build -p ${THREADS} \
        ${REF_DIR}/genome.fa \
        ${INDEX_DIR}/genome
else
    echo "HISAT2 index already exists, skipping..."
fi

# Step 2: Alignment with HISAT2
echo "Step 2/8: Aligning reads with HISAT2..."
for sample in "${SAMPLES[@]}"; do
    if [ -f "${DATA_DIR}/${sample}_R1.fastq.gz" ] && [ -f "${DATA_DIR}/${sample}_R2.fastq.gz" ]; then
        echo "Processing paired-end sample: ${sample}"
        hisat2 -p ${THREADS} \
            -x ${INDEX_DIR}/genome \
            -1 ${DATA_DIR}/${sample}_R1.fastq.gz \
            -2 ${DATA_DIR}/${sample}_R2.fastq.gz \
            -S ${ALIGN_DIR}/${sample}.sam
    elif [ -f "${DATA_DIR}/${sample}.fastq.gz" ]; then
        echo "Processing single-end sample: ${sample}"
        hisat2 -p ${THREADS} \
            -x ${INDEX_DIR}/genome \
            -U ${DATA_DIR}/${sample}.fastq.gz \
            -S ${ALIGN_DIR}/${sample}.sam
    else
        echo "No FASTQ files found for sample: ${sample}"
        exit 1
    fi
done

# Step 3: Convert SAM to BAM and sort
echo "Step 3/8: Converting SAM to BAM and sorting..."
for sample in "${SAMPLES[@]}"; do
    if [ -f "${ALIGN_DIR}/${sample}.sam" ]; then
        # Convert SAM to BAM
        samtools view -bS ${ALIGN_DIR}/${sample}.sam -o ${ALIGN_DIR}/${sample}.bam
        
        # Sort BAM
        samtools sort -o ${ALIGN_DIR}/${sample}.bam ${ALIGN_DIR}/${sample}.bam
        
        # Index BAM
        samtools index ${ALIGN_DIR}/${sample}.bam
        
        # Clean up SAM file
        rm ${ALIGN_DIR}/${sample}.sam
    fi
done

# Step 4: Prepare MISO annotations
echo "Step 4/8: Preparing MISO annotations..."
if [ ! -d "${MISO_ANN_DIR}/SE" ]; then
    # Extract alternative splicing events from GTF
    python /path/to/miso/utils/gtf2gff3.py \
        ${REF_DIR}/annotation.gtf \
        --output-dir ${MISO_ANN_DIR}
else
    echo "MISO annotations already exist, skipping..."
fi

# Step 5: Index MISO annotations
echo "Step 5/8: Indexing MISO annotations..."
for event_type in SE A5SS A3SS RI MXE; do
    if [ -d "${MISO_ANN_DIR}/${event_type}" ]; then
        index_gff --index ${MISO_ANN_DIR}/${event_type}.gff3 ${MISO_ANN_DIR}/${event_type}
    fi
done

# Step 6: Run MISO for each sample
echo "Step 6/8: Running MISO for each sample..."
for sample in "${SAMPLES[@]}"; do
    if [ -f "${ALIGN_DIR}/${sample}.bam" ]; then
        for event_type in SE A5SS A3SS RI MXE; do
            if [ -d "${MISO_ANN_DIR}/${event_type}" ]; then
                echo "Processing ${event_type} events for ${sample}"
                mkdir -p ${MISO_RES_DIR}/${sample}/${event_type}
                
                run_miso --run-indexed ${MISO_ANN_DIR}/${event_type} \
                    ${ALIGN_DIR}/${sample}.bam \
                    --output-dir ${MISO_RES_DIR}/${sample}/${event_type} \
                    --read-len 100 \
                    -p ${THREADS}
            fi
        done
    fi
done

# Step 7: Summarize MISO results
echo "Step 7/8: Summarizing MISO results..."
for sample in "${SAMPLES[@]}"; do
    for event_type in SE A5SS A3SS RI MXE; do
        if [ -d "${MISO_RES_DIR}/${sample}/${event_type}" ]; then
            summarize_miso --summarize-samples ${MISO_RES_DIR}/${sample}/${event_type} \
                ${MISO_RES_DIR}/${sample}/${event_type}_summary
        fi
    done
done

# Step 8: Compare conditions
echo "Step 8/8: Comparing conditions..."
for event_type in SE A5SS A3SS RI MXE; do
    # Prepare sample lists for each condition
    rm -f ${COMP_DIR}/condition1_${event_type}.txt ${COMP_DIR}/condition2_${event_type}.txt
    
    for i in "${!SAMPLES[@]}"; do
        if [ -d "${MISO_RES_DIR}/${SAMPLES[i]}/${event_type}" ]; then
            if [ "${CONDITIONS[i]}" = "condition1" ]; then
                echo "${MISO_RES_DIR}/${SAMPLES[i]}/${event_type}" >> ${COMP_DIR}/condition1_${event_type}.txt
            else
                echo "${MISO_RES_DIR}/${SAMPLES[i]}/${event_type}" >> ${COMP_DIR}/condition2_${event_type}.txt
            fi
        fi
    done
    
    # Compare conditions if we have samples in both groups
    if [ -s ${COMP_DIR}/condition1_${event_type}.txt ] && [ -s ${COMP_DIR}/condition2_${event_type}.txt ]; then
        compare_miso --compare-samples ${COMP_DIR}/condition1_${event_type}.txt \
            ${COMP_DIR}/condition2_${event_type}.txt \
            ${COMP_DIR}/${event_type}_comparison
    fi
done

echo "==================================================="
echo "Pipeline completed successfully!"
echo "End time: $(date)"
echo "Output files:"
echo "  - Alignments: ${ALIGN_DIR}/"
echo "  - MISO results: ${MISO_RES_DIR}/"
echo "  - Comparisons: ${COMP_DIR}/"
echo "==================================================="
