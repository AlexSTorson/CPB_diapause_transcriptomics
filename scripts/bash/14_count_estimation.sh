#!/bin/bash

#SBATCH --account=igb_fargo
#SBATCH --job-name="14_count_estimation"
#SBATCH --output="14_count_estimation.out"
#SBATCH --error="14_count_estimation.err"
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --partition=short
#SBATCH --mem=32G

# --- Variables ---
BAM_DIR="/project/igb_fargo/cpb_diapause_rnaseq/alignments_bam"
BALLGOWN_DIR="/project/igb_fargo/cpb_diapause_rnaseq/ballgown"
MERGED_GTF="/project/igb_fargo/cpb_diapause_rnaseq/assemblies/taco_assembly_merge/assembly.gtf"
RENAME_CSV="/project/igb_fargo/cpb_diapause_rnaseq/scripts/cpb_sample_rename.csv"

# --- Modules ---
module load stringtie/2.2.0

# --- Re-quantify all samples against the merged GTF ---
# -e: only estimate expression for transcripts in the merged GTF, no new assembly
# -B: write Ballgown-compatible .ctab files to the output directory
# Sample names are read from column 2 of the rename CSV, skipping the header
while IFS=',' read -r library_name sample_name; do
    [ "$library_name" = "LibraryName" ] && continue
    sample_name=$(echo "$sample_name" | tr -d '\r')

    mkdir -p "${BALLGOWN_DIR}/${sample_name}"

    stringtie \
        -p 8 \
        -e \
        -B \
        -G "$MERGED_GTF" \
        -o "${BALLGOWN_DIR}/${sample_name}/${sample_name}.gtf" \
        "${BAM_DIR}/${sample_name}.bam"

done < <(sed 's/\r//g' "$RENAME_CSV" | sed '1s/^\xef\xbb\xbf//')
