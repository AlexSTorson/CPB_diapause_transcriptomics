#!/bin/bash

#SBATCH --account=igb_fargo
#SBATCH --job-name="09_read_mapping"
#SBATCH --output="09_read_mapping.out"
#SBATCH --error="09_read_mapping.err"
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --partition=short
#SBATCH --mem=32G

# --- Variables ---
TRIMMED_READS_DIR="/project/igb_fargo/cpb_diapause_rnaseq/trimmed_reads"
SAM_DIR="/project/igb_fargo/cpb_diapause_rnaseq/alignments_sam"
INDEX_PREFIX="/project/igb_fargo/genome_references/Ldec_2.0/genome_index/Ldec_2.0_index"
RENAME_CSV="/project/igb_fargo/cpb_diapause_rnaseq/scripts/cpb_sample_rename.csv"

# --- Modules ---
module load hisat2/2.2.1

mkdir -p "$SAM_DIR"

# --- Align all samples to the genome ---
# --dta: optimizes alignments for downstream transcript assembly (required for novel discovery)
# Sample names are read from column 2 of the rename CSV, skipping the header
while IFS=',' read -r library_name sample_name; do
    [ "$library_name" = "LibraryName" ] && continue
    sample_name=$(echo "$sample_name" | tr -d '\r')

    hisat2 \
        -p 8 \
        --dta \
        --no-unal \
        --no-mixed \
        --no-discordant \
        -x "$INDEX_PREFIX" \
        -1 "${TRIMMED_READS_DIR}/${sample_name}_R1_T.fastq.gz" \
        -2 "${TRIMMED_READS_DIR}/${sample_name}_R2_T.fastq.gz" \
        -S "${SAM_DIR}/${sample_name}.sam"

done < <(sed 's/\r//g' "$RENAME_CSV" | sed '1s/^\xef\xbb\xbf//')
