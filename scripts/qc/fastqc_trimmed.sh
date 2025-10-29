#!/bin/bash
#SBATCH --job-name=fastqc_trimmed
#SBATCH --output=fastqc_trimmed_%j.out
#SBATCH --error=fastqc_trimmed_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

# Load FastQC module or set PATH if installed manually
export PATH=$HOME/FastQC:$PATH   # Or use export PATH=... if not using modules

# Input and output directories
INPUT_DIR="/mnt/vstor/CAS_BIOL_BMM30/Squid/bulk_stat_gill/fastqFiles_trimmed"
OUTPUT_DIR="/mnt/vstor/CAS_BIOL_BMM30/Squid/bulk_stat_gill/fastqc_reports_trimmed"

mkdir -p "$OUTPUT_DIR"

# Run FastQC on all fastq files (gzipped or not)
fastqc "$INPUT_DIR"/*.fastq.gz \
  --outdir="$OUTPUT_DIR" \
  --threads $SLURM_CPUS_PER_TASK
