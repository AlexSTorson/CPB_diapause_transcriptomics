#!/bin/bash

#SBATCH --account=igb_fargo
#SBATCH --job-name="12_create_mergelist"
#SBATCH --output="12_create_mergelist.out"
#SBATCH --error="12_create_mergelist.err"
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=short

# --- Variables ---
ASSEMBLIES_DIR="/project/igb_fargo/cpb_diapause_rnaseq/assemblies"

# --- Collect paths to all per-sample GTFs into a single list file ---
# StringTie --merge will read this file in the next step
find "$ASSEMBLIES_DIR" \
    -name "*.gtf" \
    -not -name "stringtie_merged.gtf" \
    > "${ASSEMBLIES_DIR}/mergelist.txt"
