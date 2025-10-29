#!/bin/bash
#SBATCH --partition=batch           # Cluster partition
#SBATCH --time=0-24:00                    # Job time: D-HH:MM
#SBATCH --cpus-per-task=8               # Number of threads
#SBATCH --mem=64G                       # Memory
#SBATCH --job-name=featureCounts_trimmed       # Job name
#SBATCH --output=%j_featureCounts_trimmed.out      # Standard output
#SBATCH --error=%j_featureCounts_trimmed.err       # Standard error

cd /mnt/vstor/CAS_BIOL_BMM30/Squid/bulk_stat_gill

# Set variables
GTF_FILE=Annotation/pita_genes_v0.3.bed_12to.gtf
OUT_DIR=featureCounts_trimmed
OUT_FILE=${OUT_DIR}/gene_id_trimmed.counts

mkdir -p $OUT_DIR

# Run featureCounts
~/software/subread-2.1.1-Linux-x86_64/bin/featureCounts -T 8 -s 2 -p \
  -a "$GTF_FILE" \
  -o "$OUT_FILE" \
  STAR_align_trimmed/*.bam
