#!/usr/bin/env bash
#SBATCH --job-name=RNAseq_AS_Analysis
#SBATCH -c 16
#SBATCH --mem=80g
#SBATCH --partition=allnodes
#SBATCH --output=RNAseq_AS.job.out
#SBATCH --error=RNAseq_AS.job.err

############################
# Step 1. Build Salmon Index
############################
salmon index \
    -t ../ref/gencode.v48.transcripts.fa \
    -i ../ref/gencode.v48.transcripts.salmon.index

########################################
# Step 2. Generate Alternative Splicing Events (SUPPA2)
########################################
suppa.py generateEvents \
    -i ../ref/gencode.v48.annotation.gtf \
    -o ../ref/gencode.v48.events \
    -e SE SS MX RI FL \
    -f ioe

# 合并所有事件文件，生成一个总的事件文件
awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; } 1 {print}' \
    ../ref/*.ioe > ../ref/gencode.v48.all.events.ioe

############################
# Step 3. 质控和数据清理
############################
cd ../rawData

# 使用 trim_galore 对 fastq 做质控和修剪
cat sample.txt | parallel -j 16 '
    fq1=./{}_RNAseq_R1.fastq.gz;
    fq2=./{}_RNAseq_R2.fastq.gz;
    trim_galore --paired -q 28 --phred33 --length 36 --stringency 3 \
        --gzip --cores 10 \
        -o ../cleanData $fq1 $fq2 \
        >> ../cleanData/{}_trim.log 2>&1
'

############################
# Step 4. 定量表达 (Salmon quant)
############################
mkdir -p ../salmon_output
cd ../cleanData
index=../ref/gencode.v48.transcripts.salmon.index/

cat ../rawData/sample.txt | parallel -j 10 '
    salmon quant \
        -i '"$index"' \
        -l ISF --gcBias \
        -1 ./{}_RNAseq_R1_val_1.fq.gz \
        -2 ./{}_RNAseq_R2_val_2.fq.gz \
        -p 10 \
        -o ../salmon_output/{}_output \
        1>../salmon_output/{}_salmon.log 2>&1
'

############################
# Step 5. 合并 Salmon 定量结果
############################
cd ../salmon_output/

multipleFieldSelection.py \
    -i FUSCCTNBC*/quant.sf \
    -k 1 -f 4 \
    -o iso_tpm.txt

# 格式化 isoform TPM 文件（只保留转录本 ID）
awk -F"\t" 'BEGIN{OFS="\t"} {split($1, id_parts, /\|/); $1=id_parts[1]; print}' \
    iso_tpm.txt > iso_tpm_formatted.txt

############################
# Step 6. SUPPA2 计算 PSI
############################
ioe_merge_file=../ref/gencode.v48.all.events.ioe
suppa.py psiPerEvent \
    -i $ioe_merge_file \
    -e iso_tpm_formatted.txt \
    -o project_events \
    1>psiPerEvent_log.txt 2>&1

echo "==== RNA-seq Alternative Splicing Pipeline Finished ===="
