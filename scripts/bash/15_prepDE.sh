#!/bin/bash

#SBATCH --account=igb_fargo
#SBATCH --job-name="15_prepDE"
#SBATCH --output="15_prepDE.out"
#SBATCH --error="15_prepDE.err"
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --partition=short

# --- Variables ---
BALLGOWN_DIR="/project/igb_fargo/cpb_diapause_rnaseq/ballgown"
SCRIPTS_DIR="/project/igb_fargo/cpb_diapause_rnaseq/scripts"
GENE_MATRIX="${BALLGOWN_DIR}/gene_count_matrix.csv"
TRANSCRIPT_MATRIX="${BALLGOWN_DIR}/transcript_count_matrix.csv"

# --- Modules ---
module load python_3/3.11.1

# --- Run prepDE.py ---
# prepDE.py expects the ballgown directory structure directly:
# ballgown/
#   sample_name/
#     sample_name.gtf
python "${SCRIPTS_DIR}/prepDE.py3" \
    -i "$BALLGOWN_DIR" \
    -g "$GENE_MATRIX" \
    -t "$TRANSCRIPT_MATRIX" \
    -l 150
