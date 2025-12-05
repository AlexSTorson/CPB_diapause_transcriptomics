#!/bin/bash

#SBATCH --time=8:00:00                  # Walltime limit (HH:MM:SS)
#SBATCH --nodes=1                        # Number of nodes
#SBATCH --ntasks-per-node=40             # 40 threads per node (20 cores x 2 threads/core)
#SBATCH --partition=short                # Partition type (e.g., short, normal, long)
#SBATCH --job-name="fastqc_pre_trim"     # Job name
#SBATCH --output="04.fastqc_pre_trimming_logs/04.fastqc_pre_trimming_%j.out"
#SBATCH --error="04.fastqc_pre_trimming_logs/04.fastqc_pre_trimming_%j.err"

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
raw_reads_dir="${base_dir}/raw_reads"
fastqc_results_dir="${base_dir}/raw_reads_fastqc_results"

# Standard log setup
scripts_dir="$(pwd)"
logs_dir="${scripts_dir}/04.fastqc_pre_trimming_logs"
sample_logs_dir="${logs_dir}/sample_logs"

# Define log files with consistent naming convention
error_log="${logs_dir}/04.fastqc_pre_trimming_errors.log"
progress_log="${logs_dir}/04.fastqc_pre_trimming_progress.log"
parallel_joblog="${logs_dir}/04.fastqc_pre_trimming_parallel_joblog.txt"

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

# --- Step 3: Find all FASTQ files ---
log_message "Locating FASTQ files"

# Add debug logging for file finding
log_message "Searching for FASTQ files in directory: $raw_reads_dir"

# List all files in the directory
log_message "All files in the directory:"
ls -l "$raw_reads_dir"

# Expanded file search patterns for raw reads
log_message "Searching for files with these patterns:"
log_message "1. *_R1.fastq.gz"
log_message "2. *_1.fastq.gz"
log_message "3. *_R1.fq.gz"
log_message "4. *_1.fq.gz"
log_message "5. *_R2.fastq.gz"
log_message "6. *_2.fastq.gz"
log_message "7. *_R2.fq.gz"
log_message "8. *_2.fq.gz"

# Find both R1 and R2 files
fastq_files=()
while IFS= read -r -d '' file; do
    fastq_files+=("$file")
    log_message "Found FASTQ file: $file"
done < <(find "$raw_reads_dir" -type f \( -name "*_R1.fastq.gz" -o -name "*_1.fastq.gz" -o -name "*_R1.fq.gz" -o -name "*_1.fq.gz" -o -name "*_R2.fastq.gz" -o -name "*_2.fastq.gz" -o -name "*_R2.fq.gz" -o -name "*_2.fq.gz" \) -print0)

# Check if any files were found
if [ ${#fastq_files[@]} -eq 0 ]; then
    log_error "No FASTQ files found in directory: $raw_reads_dir"
    log_message "Checking subdirectories..."
    
    # Search in subdirectories
    while IFS= read -r -d '' dir; do
        log_message "Examining subdirectory: $dir"
        subdir_fastq_files=()
        while IFS= read -r -d '' file; do
            subdir_fastq_files+=("$file")
            log_message "Found file in subdirectory: $file"
        done < <(find "$dir" -type f \( -name "*_R1.fastq.gz" -o -name "*_1.fastq.gz" -o -name "*_R1.fq.gz" -o -name "*_1.fq.gz" -o -name "*_R2.fastq.gz" -o -name "*_2.fastq.gz" -o -name "*_R2.fq.gz" -o -name "*_2.fq.gz" \) -print0)
        
        if [ ${#subdir_fastq_files[@]} -gt 0 ]; then
            fastq_files=("${subdir_fastq_files[@]}")
            raw_reads_dir="$dir"
            break
        fi
    done < <(find "$raw_reads_dir" -mindepth 1 -maxdepth 1 -type d -print0)
fi

# Final check
if [ ${#fastq_files[@]} -eq 0 ]; then
    log_error "No FASTQ files found in directory or its subdirectories: $raw_reads_dir"
    exit 1
fi

log_message "Found ${#fastq_files[@]} FASTQ files to process"

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

# Record end time and calculate duration
end_time=$(date +%s)
duration=$((end_time - start_time))
duration_formatted=$(printf '%dh:%dm:%ds' $((duration/3600)) $((duration%3600/60)) $((duration%60)))

# --- Step 6: Generate Summary Report ---
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
summary_file="${logs_dir}/04.fastqc_pre_trimming_summary.txt"
{
    echo "===== FastQC Pre-Trimming Summary ====="
    echo "Date: $(date)"
    echo "Execution time: $duration_formatted"
    echo "Total FASTQ files processed: $total_files"
    echo "Successfully processed: $successful_jobs"
    echo "Failed: $failed_jobs"
    
    if [ "$failed_jobs" -gt 0 ]; then
        echo "Failed files:"
        grep 'FastQC failed for file' "$error_log" | awk '{print "  - " $6}'
    fi
    
    echo "FastQC results location: $fastqc_results_dir"
    echo "========================================"
} > "$summary_file"

# --- Step 7: Log Completion ---
if [ "$failed_jobs" -gt 0 ]; then
    log_error "FastQC completed with $failed_jobs errors. See $error_log for details."
    cat "$summary_file"
    exit 1
else
    log_message "FastQC completed successfully for all $total_files FASTQ files in $duration_formatted"
    cat "$summary_file"
    exit 0
fi
