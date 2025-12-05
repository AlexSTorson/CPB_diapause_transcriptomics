#!/bin/bash
# SLURM job configuration
#SBATCH --account=igb_fargo
#SBATCH --time=00:05:00                  # Walltime limit (HH:MM:SS)
#SBATCH --nodes=1                        # Number of nodes
#SBATCH --ntasks=1                       # Single task is sufficient
#SBATCH --partition=atlas                # Partition type
#SBATCH --job-name="create_log_dirs"     # Descriptive job name
#SBATCH --output="01.create_log_dirs_logs/01.create_log_dirs_%j.out"  # Standard output log file
#SBATCH --error="01.create_log_dirs_logs/01.create_log_dirs_%j.err"   # Standard error log file
#SBATCH --mem=1G                         # Minimal memory requirement

# --- Load the Pipeline Configuration ---
# Get the directory where this script is located
script_dir="$(dirname "$(readlink -f "$0")")"
config_file="/project/igb_fargo/cpb_diapause_rnaseq/scripts/00.pipeline_config.sh"

# Verify the config file exists
if [ ! -f "$config_file" ]; then
    echo "ERROR: Configuration file not found at $config_file"
    exit 1
fi

# Load the configuration
source "$config_file"

# Script to create all log directories needed for the RNA-seq pipeline
# Run this script before submitting any jobs to ensure log directories exist

# The base directory is the current directory
base_dir="$(pwd)"

# Array of all log directories needed for the pipeline
log_dirs=(
    "01.create_log_dirs_logs"
    "02.data_download_logs"
    "02.data_download_logs/sample_logs"
    "03.genome_references_download_logs"
    "04.fastqc_pre_trimming_logs"
    "04.fastqc_pre_trimming_logs/sample_logs"
    "05.fastp_quality_trim_logs"
    "05.fastp_quality_trim_logs/sample_logs"
    "06.fastqc_post_trimming_logs"
    "06.fastqc_post_trimming_logs/sample_logs"
    "07.exon_and_splice_sites_logs"
    "08.genome_index_logs"
    "09.read_mapping_logs"
    "09.read_mapping_logs/sample_logs"
    "10.sam_to_bam_logs"
    "10.sam_to_bam_logs/sample_logs"
    "11.stringtie_logs"
    "11.stringtie_logs/sample_logs"
    "12.create_mergelist_logs"
    "13.stringtie_merge_logs"
    "13.taco_merge_logs"
    "14.gffcompare_logs"
    "15.count_estimation_logs"
    "15.count_estimation_logs/sample_logs"
    "16.prepDE_logs"
    "17.multiqc_logs"
)

echo "===== Creating log directories for RNA-seq pipeline ====="
echo "Base directory: $base_dir"
echo

# Create each directory
created=0
already_existed=0
for dir in "${log_dirs[@]}"; do
    full_path="${base_dir}/${dir}"
    
    if [ -d "$full_path" ]; then
        echo "Directory already exists: $dir"
        already_existed=$((already_existed + 1))
    else
        mkdir -p "$full_path"
        if [ $? -eq 0 ]; then
            echo "Created directory: $dir"
            created=$((created + 1))
        else
            echo "ERROR: Failed to create directory: $dir"
        fi
    fi
done

echo
echo "===== Summary ====="
echo "Total directories processed: ${#log_dirs[@]}"
echo "Directories created: $created"
echo "Directories already existing: $already_existed"
echo
echo "Log directories are now ready for pipeline execution."
