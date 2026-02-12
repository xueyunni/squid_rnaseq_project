#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --output=fastqc_%j.out
#SBATCH --error=fastqc_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

# Load FastQC module or set PATH if installed manually
export PATH=$HOME/FastQC:$PATH   # Or use export PATH=... if not using modules

# Input and output directories
INPUT_DIR=./fastqFiles
OUTPUT_DIR=./fastqc_reports

mkdir -p "$OUTPUT_DIR"

# Run FastQC on all fastq files (gzipped or not)
fastqc "$INPUT_DIR"/*.fastq.gz \
  --outdir="$OUTPUT_DIR" \
  --threads $SLURM_CPUS_PER_TASK
