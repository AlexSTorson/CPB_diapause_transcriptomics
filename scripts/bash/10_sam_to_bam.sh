#!/bin/bash

#SBATCH --account=igb_fargo
#SBATCH --job-name="10_sam_to_bam"
#SBATCH --output="10_sam_to_bam.out"
#SBATCH --error="10_sam_to_bam.err"
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --partition=short
#SBATCH --mem=32G

# --- Variables ---
SAM_DIR="/project/igb_fargo/cpb_diapause_rnaseq/alignments_sam"
BAM_DIR="/project/igb_fargo/cpb_diapause_rnaseq/alignments_bam"
RENAME_CSV="/project/igb_fargo/cpb_diapause_rnaseq/scripts/cpb_sample_rename.csv"

# --- Modules ---
module load samtools/1.17

mkdir -p "$BAM_DIR"

# --- Convert all samples from SAM to sorted, indexed BAM ---
# Sample names are read from column 2 of the rename CSV, skipping the header
while IFS=',' read -r library_name sample_name; do
    [ "$library_name" = "LibraryName" ] && continue
    sample_name=$(echo "$sample_name" | tr -d '\r')

    samtools view -@ 4 -bS "${SAM_DIR}/${sample_name}.sam" \
        | samtools sort -@ 4 -o "${BAM_DIR}/${sample_name}.bam"

    samtools index "${BAM_DIR}/${sample_name}.bam"

done < <(sed 's/\r//g' "$RENAME_CSV" | sed '1s/^\xef\xbb\xbf//')
