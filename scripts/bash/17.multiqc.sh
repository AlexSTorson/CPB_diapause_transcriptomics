#!/bin/bash

#SBATCH --job-name=multiqc               # Job name
#SBATCH --output="17.multiqc_logs/17.multiqc_%j.out"
#SBATCH --error="17.multiqc_logs/17.multiqc_%j.err"
#SBATCH --time=02:00:00                  # Time limit (hh:mm:ss)
#SBATCH --partition=short                # Partition/queue
#SBATCH --ntasks=1                       # Number of tasks
#SBATCH --cpus-per-task=4                # Number of CPU cores per task
#SBATCH --mem=16G                        # Request memory

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

# --- Standard log setup (MOVED UP) ---
current_dir="$(pwd)"
logs_dir="${current_dir}/17.multiqc_logs"

# Define log files with consistent naming convention
error_log="${logs_dir}/17.multiqc_errors.log"
progress_log="${logs_dir}/17.multiqc_progress.log"

# Clear logs at the start of the script
> "$error_log"
> "$progress_log"

# --- Centralized Variables ---
base_dir="$PROJECT_DIR"
raw_fastqc_dir="${base_dir}/raw_reads_fastqc_results"
trimmed_fastqc_dir="${base_dir}/trimmed_reads_fastqc_results"
fastp_dir="${base_dir}/trimming_summaries"
mapping_logs="${base_dir}/scripts/09.read_mapping_logs/sample_logs"  # Path to HISAT logs

# Create temporary directory for processed logs
temp_dir="${logs_dir}/temp_logs"
processed_hisat_dir="${temp_dir}/hisat_metrics"

multiqc_output_dir="${base_dir}/multiqc_report"

# --- Logging Functions ---
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$progress_log"
}

log_error() {
    echo "[ERROR] [$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$error_log"
}

# Record start time for performance tracking
start_time=$(date +%s)

# --- Step 1: Load Necessary Modules ---
log_message "Loading required modules"
$MODULE_MULTIQC

# --- Step 2: Validate Input Directories ---
log_message "Validating input directories"

# Array of directories to check
input_dirs=(
    "$raw_fastqc_dir"
    "$trimmed_fastqc_dir"
    "$fastp_dir"
    "$mapping_logs"
)

missing_dirs=0
for dir in "${input_dirs[@]}"; do
    if [ ! -d "$dir" ]; then
        log_error "Directory not found: $dir"
        missing_dirs=$((missing_dirs + 1))
    fi
done

if [ $missing_dirs -gt 0 ]; then
    log_error "Found $missing_dirs missing input directories"
    log_message "Will continue with available directories"
fi

# --- Step 3: Create Output and Temporary Directories ---
log_message "Creating MultiQC output directory: $multiqc_output_dir"
mkdir -p "$multiqc_output_dir"
if [ $? -ne 0 ]; then
    log_error "Failed to create output directory: $multiqc_output_dir"
    exit 1
fi

# Create temporary directory for processed logs
log_message "Creating temporary directory for processed logs: $temp_dir"
mkdir -p "$processed_hisat_dir"
if [ $? -ne 0 ]; then
    log_error "Failed to create temporary directory: $processed_hisat_dir"
    exit 1
fi

# --- Step 4: Process HISAT2 logs to avoid redundancy ---
log_message "Processing HISAT2 logs to avoid redundancy"

# Find and copy only *_metrics.txt files (not the full output files)
metrics_count=0
while IFS= read -r file; do
    if [ -f "$file" ]; then
        sample_name=$(basename "$file" _metrics.txt)
        cp "$file" "${processed_hisat_dir}/${sample_name}_hisat2.txt"
        metrics_count=$((metrics_count + 1))
    fi
done < <(find "$mapping_logs" -name "*_metrics.txt" -type f)

log_message "Processed $metrics_count HISAT2 metrics files"

# --- Step 5: Create MultiQC Configuration File ---
log_message "Creating MultiQC configuration file"
multiqc_config="${logs_dir}/multiqc_config.yaml"

cat > "$multiqc_config" <<EOL
# MultiQC Configuration for RNA-seq Pipeline
title: "RNA-seq Pipeline QC Report"
subtitle: "Quality control metrics for the RNA-seq pipeline"
intro_text: "This report summarizes quality control metrics from raw reads through alignment."

# Always put HISAT2 reports at the end
report_section_order:
  fastqc_pre_trimming:
    order: 10
  fastp:
    order: 20
  fastqc_post_trimming:
    order: 30
  hisat2:
    order: 900

# Custom module order
module_order:
  - fastqc:
      name: "Pre-Trimming FastQC"
      path_filters:
        - "*/raw_reads_fastqc_results/*"
  - fastp
  - fastqc:
      name: "Post-Trimming FastQC"
      path_filters:
        - "*/trimmed_reads_fastqc_results/*"
  - hisat2

# Custom sample naming
sample_names_replace_dash: False
sample_names_replace: 
  "_fastqc.zip": ""
  "_1_fastqc.zip": "_R1"
  "_2_fastqc.zip": "_R2"
  "_R1_T_fastqc.zip": "_R1_trimmed"
  "_R2_T_fastqc.zip": "_R2_trimmed"

# File name cleaning
fn_clean_exts:
  - ".fastq.gz"
  - ".fq.gz"
  - "_fastqc.zip"
  - ".bam"
  - ".sam"
  - "_metrics.txt"
  - "_hisat2.txt"
  - "_full_output.txt"

# Remove redundant modules by explicitly excluding them
hisat2/full_output:
  fn_re: '.*_full_output\.txt$'
  exclude: true

# Make sure report is in the correct order
show_analysis_paths: False
show_analysis_time: False
EOL

# --- Step 6: Create File Lists for Each Component ---
log_message "Creating file lists for each component"

# 1. Pre-Trimming FastQC
raw_fastqc_list="${logs_dir}/raw_fastqc_list.txt"
find "$raw_fastqc_dir" -name "*_fastqc.zip" -type f > "$raw_fastqc_list"
raw_fastqc_count=$(wc -l < "$raw_fastqc_list")

# 2. Fastp Summaries
fastp_list="${logs_dir}/fastp_list.txt"
find "$fastp_dir" -name "*.json" -type f > "$fastp_list"
fastp_count=$(wc -l < "$fastp_list")

# 3. Post-Trimming FastQC
trimmed_fastqc_list="${logs_dir}/trimmed_fastqc_list.txt"
find "$trimmed_fastqc_dir" -name "*_fastqc.zip" -type f > "$trimmed_fastqc_list"
trimmed_fastqc_count=$(wc -l < "$trimmed_fastqc_list")

# 4. HISAT2 Logs (processed)
hisat_list="${logs_dir}/hisat_list.txt"
find "$processed_hisat_dir" -name "*_hisat2.txt" -type f > "$hisat_list"
hisat_count=$(wc -l < "$hisat_list")

# Log the counts
log_message "Found $raw_fastqc_count pre-trimming FastQC files"
log_message "Found $fastp_count fastp summary files"
log_message "Found $trimmed_fastqc_count post-trimming FastQC files"
log_message "Found $hisat_count processed HISAT2 log files"

# --- Step 7: Run MultiQC with Custom Configuration ---
log_message "Running MultiQC with custom configuration"

multiqc_cmd="multiqc -c $multiqc_config \
             -o $multiqc_output_dir \
             -n multiqc_report \
             -f \
             $raw_fastqc_dir $fastp_dir $trimmed_fastqc_dir $processed_hisat_dir"

log_message "Command: $multiqc_cmd"
$multiqc_cmd &> "${logs_dir}/multiqc_output.log"

# Check if MultiQC ran successfully
if [ $? -ne 0 ]; then
    log_error "MultiQC failed. See ${logs_dir}/multiqc_output.log for details."
    
    # Try running without the config file as a fallback
    log_message "Trying to run MultiQC without custom configuration"
    
    multiqc -o "$multiqc_output_dir" -n "multiqc_report_fallback" -f \
            $raw_fastqc_dir $fastp_dir $trimmed_fastqc_dir $processed_hisat_dir \
            &> "${logs_dir}/multiqc_fallback.log"
    
    if [ $? -ne 0 ]; then
        log_error "Fallback MultiQC also failed. Check ${logs_dir}/multiqc_fallback.log"
        exit 1
    else
        log_message "Fallback MultiQC report created successfully (but without custom ordering)"
    fi
else
    log_message "MultiQC report created successfully"
fi

# --- Step 8: Verify Output Files ---
log_message "Verifying MultiQC output files"

# Check for primary report
primary_html="${multiqc_output_dir}/multiqc_report.html"
fallback_html="${multiqc_output_dir}/multiqc_report_fallback.html"

if [ -f "$primary_html" ]; then
    html_size=$(stat -c%s "$primary_html")
    html_size_mb=$(echo "scale=2; $html_size / 1048576" | bc)
    log_message "MultiQC report size: ${html_size_mb} MB"
    report_path="$primary_html"
elif [ -f "$fallback_html" ]; then
    html_size=$(stat -c%s "$fallback_html")
    html_size_mb=$(echo "scale=2; $html_size / 1048576" | bc)
    log_message "Fallback MultiQC report size: ${html_size_mb} MB"
    report_path="$fallback_html"
else
    log_error "No MultiQC report was created"
    exit 1
fi

# Count data files
data_dir="${multiqc_output_dir}/multiqc_data"
if [ -d "$data_dir" ]; then
    data_files=$(find "$data_dir" -type f | wc -l)
    log_message "MultiQC data directory contains $data_files files"
else
    log_error "MultiQC data directory not found: $data_dir"
fi

# Record end time and calculate duration
end_time=$(date +%s)
duration=$((end_time - start_time))
duration_formatted=$(printf '%dm:%ds' $((duration/60)) $((duration%60)))

# --- Step 9: Clean up temporary files ---
log_message "Cleaning up temporary files"
if [ -d "$temp_dir" ]; then
    rm -rf "$temp_dir"
    log_message "Removed temporary directory: $temp_dir"
fi

# --- Step 10: Generate Summary Report ---
summary_file="${logs_dir}/17.multiqc_summary.txt"
{
    echo "===== MultiQC Summary ====="
    echo "Date: $(date)"
    echo "Execution time: $duration_formatted"
    echo
    echo "Input Files:"
    echo "  Pre-Trimming FastQC files: $raw_fastqc_count"
    echo "  Fastp summary files: $fastp_count"
    echo "  Post-Trimming FastQC files: $trimmed_fastqc_count"
    echo "  HISAT2 log files: $hisat_count"
    echo
    echo "Output:"
    echo "  MultiQC output directory: $multiqc_output_dir"
    echo "  Report: $(basename "$report_path")"
    echo "  Report size: ${html_size_mb} MB"
    echo "  Data files: $data_files"
    echo
    echo "======================================="
} > "$summary_file"

# --- Step 11: Log Completion ---
log_message "MultiQC completed successfully in $duration_formatted"
log_message "Report is available at: $report_path"
cat "$summary_file"
exit 0
