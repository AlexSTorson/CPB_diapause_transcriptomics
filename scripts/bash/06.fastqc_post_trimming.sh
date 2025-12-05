#!/bin/bash

#SBATCH --time=8:00:00                  # Walltime limit (HH:MM:SS)
#SBATCH --nodes=1                        # Number of nodes
#SBATCH --ntasks-per-node=40             # 40 threads per node (20 cores x 2 threads/core)
#SBATCH --partition=short                # Partition type (e.g., short, normal, long)
#SBATCH --job-name="fastqc_post_trim"    # Job name
#SBATCH --output="06.fastqc_post_trimming_logs/06.fastqc_post_trimming_%j.out"
#SBATCH --error="06.fastqc_post_trimming_logs/06.fastqc_post_trimming_%j.err"

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

# --- Centralized Variables ---
base_dir="$PROJECT_DIR"
trimmed_reads_dir="${base_dir}/trimmed_reads"
fastqc_results_dir="${base_dir}/trimmed_reads_fastqc_results"

# Standard log setup
scripts_dir="$(pwd)"
logs_dir="${scripts_dir}/06.fastqc_post_trimming_logs"
sample_logs_dir="${logs_dir}/sample_logs"

# Define log files with consistent naming convention
error_log="${logs_dir}/06.fastqc_post_trimming_errors.log"
progress_log="${logs_dir}/06.fastqc_post_trimming_progress.log"
parallel_joblog="${logs_dir}/06.fastqc_post_trimming_parallel_joblog.txt"

# Performance tuning
threads_per_fastqc=2  # FastQC doesn't benefit from many threads, use fewer per job
max_parallel_jobs=20  # Run multiple FastQC instances in parallel (40/2=20)

# Create required directories
mkdir -p "$logs_dir" "$sample_logs_dir"

# Clear logs at the start of the script
> "$error_log"
> "$progress_log"

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
$MODULE_FASTQC
$MODULE_PARALLEL

# --- Step 2: Create Output Directory (if it does not exist) ---
log_message "Creating output directory: $fastqc_results_dir"
mkdir -p "$fastqc_results_dir"
if [ $? -ne 0 ]; then
    log_error "Failed to create output directory: $fastqc_results_dir"
    exit 1
fi

# --- Step 3: Find all trimmed FASTQ files ---
log_message "Locating trimmed FASTQ files"

# Add debug logging for file finding
log_message "Searching for trimmed FASTQ files in directory: $trimmed_reads_dir"

# List all files in the directory
log_message "All files in the directory:"
ls -l "$trimmed_reads_dir"

# Expanded file search patterns for trimmed reads
log_message "Searching for files with these patterns:"
log_message "1. *_R1_T.fastq.gz"
log_message "2. *_1_T.fastq.gz"
log_message "3. *_R1_T.fq.gz"
log_message "4. *_1_T.fq.gz"
log_message "5. *_R2_T.fastq.gz"
log_message "6. *_2_T.fastq.gz"
log_message "7. *_R2_T.fq.gz"
log_message "8. *_2_T.fq.gz"

# Find both R1 and R2 trimmed files
fastq_files=()
while IFS= read -r -d '' file; do
    fastq_files+=("$file")
    log_message "Found trimmed FASTQ file: $file"
done < <(find "$trimmed_reads_dir" -type f \( -name "*_R1_T.fastq.gz" -o -name "*_1_T.fastq.gz" -o -name "*_R1_T.fq.gz" -o -name "*_1_T.fq.gz" -o -name "*_R2_T.fastq.gz" -o -name "*_2_T.fastq.gz" -o -name "*_R2_T.fq.gz" -o -name "*_2_T.fq.gz" \) -print0)

# Check if any files were found
if [ ${#fastq_files[@]} -eq 0 ]; then
    log_error "No trimmed FASTQ files found in directory: $trimmed_reads_dir"
    log_message "Checking subdirectories..."
    
    # Search in subdirectories
    while IFS= read -r -d '' dir; do
        log_message "Examining subdirectory: $dir"
        subdir_fastq_files=()
        while IFS= read -r -d '' file; do
            subdir_fastq_files+=("$file")
            log_message "Found file in subdirectory: $file"
        done < <(find "$dir" -type f \( -name "*_R1_T.fastq.gz" -o -name "*_1_T.fastq.gz" -o -name "*_R1_T.fq.gz" -o -name "*_1_T.fq.gz" -o -name "*_R2_T.fastq.gz" -o -name "*_2_T.fastq.gz" -o -name "*_R2_T.fq.gz" -o -name "*_2_T.fq.gz" \) -print0)
        
        if [ ${#subdir_fastq_files[@]} -gt 0 ]; then
            fastq_files=("${subdir_fastq_files[@]}")
            trimmed_reads_dir="$dir"
            break
        fi
    done < <(find "$trimmed_reads_dir" -mindepth 1 -maxdepth 1 -type d -print0)
fi

# Final check
if [ ${#fastq_files[@]} -eq 0 ]; then
    log_error "No trimmed FASTQ files found in directory or its subdirectories: $trimmed_reads_dir"
    exit 1
fi

log_message "Found ${#fastq_files[@]} trimmed FASTQ files to process"

# --- Step 4: Define FastQC function for parallel execution ---
run_fastqc() {
    local fastq_file="$1"
    local sample_name=$(basename "$fastq_file")
    local sample_log="${sample_logs_dir}/${sample_name}.log"
    
    echo "Processing: $sample_name" > "$sample_log"
    
    # Run FastQC with specified number of threads
    fastqc -t "$threads_per_fastqc" "$fastq_file" --outdir="$fastqc_results_dir" &>> "$sample_log"
    
    local exit_code=$?
    if [ $exit_code -ne 0 ]; then
        echo "FastQC failed for file: $sample_name (Exit code: $exit_code)" >> "$error_log"
        return 1
    else
        echo "FastQC completed for file: $sample_name" >> "$progress_log"
        return 0
    fi
}

export -f run_fastqc
export fastqc_results_dir error_log progress_log threads_per_fastqc sample_logs_dir

# --- Step 5: Run FastQC in parallel ---
log_message "Starting parallel FastQC processing with $max_parallel_jobs jobs"

parallel -u --jobs "$max_parallel_jobs" --no-notice --joblog "$parallel_joblog" run_fastqc ::: "${fastq_files[@]}"

# --- Step 6: Validate Output Files ---
log_message "Validating FastQC output files"
expected_html_files=${#fastq_files[@]}
expected_zip_files=${#fastq_files[@]}

actual_html_files=$(find "$fastqc_results_dir" -name "*.html" | wc -l)
actual_zip_files=$(find "$fastqc_results_dir" -name "*.zip" | wc -l)

log_message "Expected HTML reports: $expected_html_files, Found: $actual_html_files"
log_message "Expected ZIP archives: $expected_zip_files, Found: $actual_zip_files"

# Record end time and calculate duration
end_time=$(date +%s)
duration=$((end_time - start_time))
duration_formatted=$(printf '%dh:%dm:%ds' $((duration/3600)) $((duration%3600/60)) $((duration%60)))

# --- Step 7: Generate Summary Report ---
log_message "Generating summary report"

# Count successful and failed jobs - Use a safer approach
total_files=${#fastq_files[@]}
# Ensure we have a valid integer for failed_jobs
failed_jobs=0
if [ -f "$parallel_joblog" ]; then
    failed_count=$(grep -c 'exit code: [^0]' "$parallel_joblog" 2>/dev/null || echo "0")
    # Make sure it's a valid integer before using it
    if [[ "$failed_count" =~ ^[0-9]+$ ]]; then
        failed_jobs="$failed_count"
    fi
fi
successful_jobs=$((total_files - failed_jobs))

# Create a summary report
summary_file="${logs_dir}/06.fastqc_post_trimming_summary.txt"
{
    echo "===== FastQC Post-Trimming Summary ====="
    echo "Date: $(date)"
    echo "Execution time: $duration_formatted"
    echo
    echo "Processing Statistics:"
    echo "  Total FASTQ files processed: $total_files"
    echo "  Successfully processed: $successful_jobs"
    echo "  Failed: $failed_jobs"
    echo
    echo "Output Files:"
    echo "  Expected HTML reports: $expected_html_files"
    echo "  Actual HTML reports: $actual_html_files"
    echo "  Expected ZIP archives: $expected_zip_files"
    echo "  Actual ZIP archives: $actual_zip_files"
    echo
    if [ "$failed_jobs" -gt 0 ]; then
        echo "Failed files:"
        grep 'FastQC failed for file' "$error_log" | awk '{print "  - " $6}'
    fi
    echo
    echo "FastQC results location: $fastqc_results_dir"
    echo "========================================"
} > "$summary_file"

# --- Step 8: Log Completion ---
if [ "$failed_jobs" -gt 0 ] || [ "$actual_html_files" -lt "$expected_html_files" ] || [ "$actual_zip_files" -lt "$expected_zip_files" ]; then
    log_error "FastQC post-trimming completed with errors. See $error_log for details."
    cat "$summary_file"
    exit 1
else
    log_message "FastQC post-trimming completed successfully for all $total_files files in $duration_formatted"
    cat "$summary_file"
    exit 0
fi
