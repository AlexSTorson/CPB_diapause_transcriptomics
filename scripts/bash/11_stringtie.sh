#!/bin/bash

#SBATCH --account=igb_fargo
#SBATCH --job-name="11_stringtie"
#SBATCH --output="11_stringtie.out"
#SBATCH --error="11_stringtie.err"
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --partition=short
#SBATCH --mem=32G

# --- Variables ---
BAM_DIR="/project/igb_fargo/cpb_diapause_rnaseq/alignments_bam"
ASSEMBLIES_DIR="/project/igb_fargo/cpb_diapause_rnaseq/assemblies"
GTF="/project/igb_fargo/genome_references/Ldec_2.0/GCF_000500325.1_Ldec_2.0_genomic.gtf"
RENAME_CSV="/project/igb_fargo/cpb_diapause_rnaseq/scripts/cpb_sample_rename.csv"

# --- Modules ---
module load stringtie/2.2.0

mkdir -p "$ASSEMBLIES_DIR"

# --- Assemble transcripts for all samples ---
# -G: use reference GTF to guide assembly
# Novel transcripts beyond the reference are assembled by default (no -e flag)
# Sample names are read from column 2 of the rename CSV, skipping the header
while IFS=',' read -r library_name sample_name; do
    [ "$library_name" = "LibraryName" ] && continue
    sample_name=$(echo "$sample_name" | tr -d '\r')

    stringtie \
        -p 8 \
        -G "$GTF" \
        -o "${ASSEMBLIES_DIR}/${sample_name}.gtf" \
        "${BAM_DIR}/${sample_name}.bam"

done < <(sed 's/\r//g' "$RENAME_CSV" | sed '1s/^\xef\xbb\xbf//')
