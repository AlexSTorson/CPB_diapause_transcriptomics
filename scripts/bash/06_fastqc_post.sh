#!/bin/bash

#SBATCH --account=igb_fargo
#SBATCH --job-name="06_fastqc_post"
#SBATCH --output="06_fastqc_post.out"
#SBATCH --error="06_fastqc_post.err"
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --partition=short

# --- Variables ---
TRIMMED_READS_DIR="/project/igb_fargo/cpb_diapause_rnaseq/trimmed_reads"
FASTQC_DIR="/project/igb_fargo/cpb_diapause_rnaseq/trimmed_reads_fastqc"
RENAME_CSV="/project/igb_fargo/cpb_diapause_rnaseq/scripts/cpb_sample_rename.csv"

# --- Modules ---
module load fastqc/0.11.9

mkdir -p "$FASTQC_DIR"

# --- Run FastQC on all trimmed samples ---
# Sample names are read from column 2 of the rename CSV, skipping the header
while IFS=',' read -r library_name sample_name; do
    [ "$library_name" = "LibraryName" ] && continue
    sample_name=$(echo "$sample_name" | tr -d '\r')

    fastqc \
        --outdir "$FASTQC_DIR" \
        --threads 2 \
        "${TRIMMED_READS_DIR}/${sample_name}_R1_T.fastq.gz" \
        "${TRIMMED_READS_DIR}/${sample_name}_R2_T.fastq.gz"

done < <(sed 's/\r//g' "$RENAME_CSV" | sed '1s/^\xef\xbb\xbf//')
