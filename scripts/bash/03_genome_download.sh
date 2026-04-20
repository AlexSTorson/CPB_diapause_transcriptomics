#!/bin/bash

#SBATCH --account=igb_fargo
#SBATCH --job-name="03_genome_download"
#SBATCH --output="03_genome_download.out"
#SBATCH --error="03_genome_download.err"
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=short

# --- Variables ---
REFERENCE_DIR="/project/igb_fargo/genome_references/Ldec_2.0"
FTP_BASE="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/500/325/GCF_000500325.1_Ldec_2.0"

mkdir -p "$REFERENCE_DIR"

# --- Download genome FASTA, GTF annotation, and feature table from NCBI FTP ---
wget -P "$REFERENCE_DIR" "${FTP_BASE}/GCF_000500325.1_Ldec_2.0_genomic.fna.gz"
wget -P "$REFERENCE_DIR" "${FTP_BASE}/GCF_000500325.1_Ldec_2.0_genomic.gtf.gz"
wget -P "$REFERENCE_DIR" "${FTP_BASE}/GCF_000500325.1_Ldec_2.0_feature_table.txt.gz"

# --- Decompress ---
gunzip "${REFERENCE_DIR}/GCF_000500325.1_Ldec_2.0_genomic.fna.gz"
gunzip "${REFERENCE_DIR}/GCF_000500325.1_Ldec_2.0_genomic.gtf.gz"
gunzip "${REFERENCE_DIR}/GCF_000500325.1_Ldec_2.0_feature_table.txt.gz"
