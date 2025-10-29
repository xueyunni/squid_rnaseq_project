#!/bin/bash

# ====== USER CONFIGURATION ======
GENOME_FA="genome/Dpe_v2.1.fasta"         # e.g., Homo_sapiens.GRCh38.dna.primary_assembly.fa
GTF_FILE="annotation/pita_genes_v0.3.bed_12to.gtf"          # e.g., Homo_sapiens.GRCh38.110.gtf
OUTDIR="indices"                    # Output directory for index
THREADS=8                             # Number of threads
READ_LENGTH=100                        # Read length of your FASTQ files

# ====== CREATE OUTPUT DIR ======
mkdir -p "$OUTDIR"

# ====== RUN STAR INDEXING ======
STAR \
  --runThreadN "$THREADS" \
  --runMode genomeGenerate \
  --genomeDir "$OUTDIR" \
  --genomeFastaFiles "$GENOME_FA" \
  --sjdbGTFfile "$GTF_FILE" \
  --sjdbOverhang $((READ_LENGTH - 1))
