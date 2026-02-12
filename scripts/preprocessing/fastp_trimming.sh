#!/bin/bash
#SBATCH --job-name=fastp_bulk
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=fastp_%j.out
#SBATCH --error=fastp_%j.err

# -------------------------------
# Load conda environment
# -------------------------------
source ~/.bashrc
conda activate rnaseq_env

# -------------------------------
# Input and output directories
# -------------------------------
INPUT_DIR="/mnt/vstor/CAS_BIOL_BMM30/Squid/bulk_stat_gill/fastqFiles"
OUT_DIR="/mnt/vstor/CAS_BIOL_BMM30/Squid/bulk_stat_gill/fastqFiles_trimmed"

mkdir -p $OUT_DIR

# -------------------------------
# Loop over paired-end samples
# Assumes filenames: SAMPLE_R1.fastq.gz and SAMPLE_R2.fastq.gz
# -------------------------------
for R1 in $INPUT_DIR/*_R1.fastq.gz; do
    SAMPLE=$(basename $R1 _R1.fastq.gz)
    R2="$INPUT_DIR/${SAMPLE}_R2.fastq.gz"

    echo "Processing sample: $SAMPLE"

    fastp \
      -i $R1 \
      -I $R2 \
      -o $OUT_DIR/${SAMPLE}_R1.trimmed.fastq.gz \
      -O $OUT_DIR/${SAMPLE}_R2.trimmed.fastq.gz \
      -h $OUT_DIR/${SAMPLE}_fastp.html \
      -j $OUT_DIR/${SAMPLE}_fastp.json \
      -q 20 \
      -u 30 \
      -l 50 \
      --detect_adapter_for_pe \
      --thread 8
done

echo "All samples processed!"
