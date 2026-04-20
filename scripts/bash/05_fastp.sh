#!/bin/bash

#SBATCH --account=igb_fargo
#SBATCH --job-name="05_fastp"
#SBATCH --output="05_fastp.out"
#SBATCH --error="05_fastp.err"
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --partition=short

# --- Variables ---
RAW_READS_DIR="/project/igb_fargo/cpb_diapause_rnaseq/raw_reads"
TRIMMED_READS_DIR="/project/igb_fargo/cpb_diapause_rnaseq/trimmed_reads"
SUMMARIES_DIR="/project/igb_fargo/cpb_diapause_rnaseq/trimming_summaries"
RENAME_CSV="/project/igb_fargo/cpb_diapause_rnaseq/scripts/cpb_sample_rename.csv"

# --- Modules ---
module load miniconda
source activate fastpEnv

mkdir -p "$TRIMMED_READS_DIR" "$SUMMARIES_DIR"

# --- Run fastp on all samples ---
# Sample names are read from column 2 of the rename CSV, skipping the header
while IFS=',' read -r library_name sample_name; do
    [ "$library_name" = "LibraryName" ] && continue
    sample_name=$(echo "$sample_name" | tr -d '\r')

    fastp \
        -i "${RAW_READS_DIR}/${sample_name}_R1.fastq.gz" \
        -I "${RAW_READS_DIR}/${sample_name}_R2.fastq.gz" \
        -o "${TRIMMED_READS_DIR}/${sample_name}_R1_T.fastq.gz" \
        -O "${TRIMMED_READS_DIR}/${sample_name}_R2_T.fastq.gz" \
        -w 4 \
        --detect_adapter_for_pe \
        --json "${SUMMARIES_DIR}/${sample_name}_fastp.json" \
        --html "${SUMMARIES_DIR}/${sample_name}_fastp.html"

done < <(sed 's/\r//g' "$RENAME_CSV" | sed '1s/^\xef\xbb\xbf//')
