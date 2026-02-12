#!/bin/bash
#SBATCH --partition=batch           # Cluster partition
#SBATCH --time=0-24:00                    # Job time: D-HH:MM
#SBATCH --cpus-per-task=8               # Number of threads
#SBATCH --mem=64G                       # Memory
#SBATCH --job-name=STAR_indexing        # Job name
#SBATCH --output=%j_STAR_index.out      # Standard output
#SBATCH --error=%j_STAR_index.err       # Standard error

# Move to working directory (adjust as needed)
cd /mnt/vstor/CAS_BIOL_BMM30/Squid/bulk_stat_gill

# Load STAR module
module load STAR

# Create output directory if it doesn't exist
mkdir -p indices/Dpe_v2.1_index

# Run STAR indexing
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir indices/Dpe_v2.1_index \
     --genomeFastaFiles Genome/Dpe_v2.1.fasta \
     --sjdbGTFfile Annotation/pita_genes_v0.3.bed_12to.gtf \
     --sjdbOverhang 99
