#!/bin/bash

#SBATCH --account=igb_fargo
#SBATCH --job-name="07_extract_splice_exons"
#SBATCH --output="07_extract_splice_exons.out"
#SBATCH --error="07_extract_splice_exons.err"
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --partition=short

# --- Variables ---
REFERENCE_DIR="/project/igb_fargo/genome_references/Ldec_2.0"
GTF="${REFERENCE_DIR}/GCF_000500325.1_Ldec_2.0_genomic.gtf"
SPLICE_SITES="${REFERENCE_DIR}/Ldec_2.0.ss"
EXONS="${REFERENCE_DIR}/Ldec_2.0.exon"
SCRIPTS_DIR="/project/igb_fargo/cpb_diapause_rnaseq/scripts"

# --- Modules ---
module load python_3/3.11.1

# --- Extract splice sites and exons from the reference GTF ---
python "${SCRIPTS_DIR}/hisat2_extract_splice_sites.py" "$GTF" > "$SPLICE_SITES"
python "${SCRIPTS_DIR}/hisat2_extract_exons.py" "$GTF" > "$EXONS"
