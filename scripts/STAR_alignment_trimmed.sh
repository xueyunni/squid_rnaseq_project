#!/bin/bash
#SBATCH --partition=batch           # Cluster partition
#SBATCH --time=0-24:00                    # Job time: D-HH:MM
#SBATCH --cpus-per-task=8               # Number of threads
#SBATCH --mem=64G                       # Memory
#SBATCH --job-name=STAR_alignment_trimmed        # Job name
#SBATCH --output=%j_STAR_alignment_trimmed.out      # Standard output
#SBATCH --error=%j_STAR_alignment_trimmed.err       # Standard error

# Move to working directory (adjust as needed)
cd /mnt/vstor/CAS_BIOL_BMM30/Squid/bulk_stat_gill

# Load STAR module
module load STAR

# Set variables
GENOME_DIR=indices/Dpe_v2.1_index
FASTQ_DIR=fastqFiles_trimmed
OUT_DIR=STAR_align_trimmed
THREADS=8

mkdir -p $OUT_DIR

# Loop over paired-end FASTQ files
for  R1 in  ${FASTQ_DIR}/*_R1.trimmed.fastq.gz
do
    # Derive sample name
    SAMPLE=$(basename $R1 _R1.trimmed.fastq.gz)
    R2=${FASTQ_DIR}/${SAMPLE}_R2.trimmed.fastq.gz

    # Run STAR
    STAR --runThreadN $THREADS \
         --genomeDir $GENOME_DIR \
         --readFilesIn $R1 $R2 \
         --readFilesCommand zcat \
         --outFileNamePrefix ${OUT_DIR}/${SAMPLE}_ \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMunmapped Within \
         --outSAMattributes Standard
done

