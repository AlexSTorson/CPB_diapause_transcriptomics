#!/bin/bash

#SBATCH --account=igb_fargo
#SBATCH --job-name="04_fastqc_pre"
#SBATCH --output="04_fastqc_pre.out"
#SBATCH --error="04_fastqc_pre.err"
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --partition=short

# --- Variables ---
RAW_READS_DIR="/project/igb_fargo/cpb_diapause_rnaseq/raw_reads"
FASTQC_DIR="/project/igb_fargo/cpb_diapause_rnaseq/raw_reads_fastqc"
RENAME_CSV="/project/igb_fargo/cpb_diapause_rnaseq/scripts/cpb_sample_rename.csv"

# --- Modules ---
module load fastqc/0.11.9

mkdir -p "$FASTQC_DIR"

# --- Run FastQC on all samples ---
# Sample names are read from column 2 of the rename CSV, skipping the header
while IFS=',' read -r library_name sample_name; do
    [ "$library_name" = "LibraryName" ] && continue
    sample_name=$(echo "$sample_name" | tr -d '\r')

    fastqc \
        --outdir "$FASTQC_DIR" \
        --threads 2 \
        "${RAW_READS_DIR}/${sample_name}_R1.fastq.gz" \
        "${RAW_READS_DIR}/${sample_name}_R2.fastq.gz"

done < <(sed 's/\r//g' "$RENAME_CSV" | sed '1s/^\xef\xbb\xbf//')
