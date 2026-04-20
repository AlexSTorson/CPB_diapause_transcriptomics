#!/bin/bash

#SBATCH --account=igb_fargo
#SBATCH --job-name="08_genome_index"
#SBATCH --output="08_genome_index.out"
#SBATCH --error="08_genome_index.err"
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --partition=short
#SBATCH --mem=64G

# --- Variables ---
REFERENCE_DIR="/project/igb_fargo/genome_references/Ldec_2.0"
GENOME_FASTA="${REFERENCE_DIR}/GCF_000500325.1_Ldec_2.0_genomic.fna"
SPLICE_SITES="${REFERENCE_DIR}/Ldec_2.0.ss"
EXONS="${REFERENCE_DIR}/Ldec_2.0.exon"
INDEX_DIR="${REFERENCE_DIR}/genome_index"
INDEX_PREFIX="${INDEX_DIR}/Ldec_2.0_index"

# --- Modules ---
module load hisat2/2.2.1

mkdir -p "$INDEX_DIR"

# --- Build the splice-aware genome index ---
hisat2-build \
    -p 20 \
    --ss "$SPLICE_SITES" \
    --exon "$EXONS" \
    "$GENOME_FASTA" \
    "$INDEX_PREFIX"
