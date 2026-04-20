#!/bin/bash

#SBATCH --account=igb_fargo
#SBATCH --job-name="13_taco_merge"
#SBATCH --output="13_taco_merge.out"
#SBATCH --error="13_taco_merge.err"
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --partition=short
#SBATCH --mem=120G

# --- Variables ---
ASSEMBLIES_DIR="/project/igb_fargo/cpb_diapause_rnaseq/assemblies"
TACO_OUTPUT_DIR="/project/igb_fargo/cpb_diapause_rnaseq/assemblies/taco_assembly_merge"
MERGELIST="${ASSEMBLIES_DIR}/mergelist.txt"
TACO_EXEC="/project/igb_fargo/programs/taco_run"

# --- Run TACO ---
# TACO requires the output directory to not already exist
rm -rf "$TACO_OUTPUT_DIR"

"$TACO_EXEC" \
    -p 20 \
    -o "$TACO_OUTPUT_DIR" \
    "$MERGELIST"
