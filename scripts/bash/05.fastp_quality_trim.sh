#!/bin/bash

#SBATCH --time=48:00:00                  # Walltime limit (HH:MM:SS)
#SBATCH --nodes=1                        # Number of nodes
#SBATCH --ntasks-per-node=40             # 40 threads per node
#SBATCH --partition=short                # Partition type
#SBATCH --job-name="fastp_quality_trim"  # Descriptive job name
#SBATCH --output="05.fastp_quality_trim_logs/05.fastp_quality_trim_%j.out"
#SBATCH --error="05.fastp_quality_trim_logs/05.fastp_quality_trim_%j.err"

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
trimmed_reads_dir="${base_dir}/trimmed_reads"
summaries_dir="${base_dir}/trimming_summaries"
conda_env="$CONDA_ENV_FASTP"  # Name of your conda environment for fastp

# Standard log setup
scripts_dir="$(pwd)"
logs_dir="${scripts_dir}/05.fastp_quality_trim_logs"
sample_logs_dir="${logs_dir}/sample_logs"

# Define log files with consistent naming convention
error_log="${logs_dir}/05.fastp_quality_trim_errors.log"
progress_log="${logs_dir}/05.fastp_quality_trim_progress.log"
parallel_joblog="${logs_dir}/05.fastp_quality_trim_parallel_joblog.txt"

# Performance tuning
threads_per_job=4  # Threads per fastp job
max_parallel_jobs=10  # Number of parallel jobs (40/4=10)

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

# --- Step 1: Load Required Modules ---
log_message "Loading required modules"
$MODULE_MINICONDA || { log_error "Failed to load Miniconda module."; exit 1; }
$MODULE_PARALLEL || { log_error "Failed to load Parallel module."; exit 1; }

# --- Step 2: Activate Conda Environment ---
log_message "Activating conda environment: $conda_env"
source activate "$conda_env" || { log_error "Failed to activate conda environment: $conda_env"; exit 1; }

# --- Step 3: Create Output Directories ---
log_message "Creating output directories"
mkdir -p "$trimmed_reads_dir" "$summaries_dir"
if [ $? -ne 0 ]; then
    log_error "Failed to create output directories."
    exit 1
fi

# --- Step 4: Find Raw Reads for Processing ---
log_message "Locating raw read files for processing"

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

r1_files=()
while IFS= read -r -d '' file; do
    r1_files+=("$file")
    log_message "Found R1 file: $file"
done < <(find "$raw_reads_dir" -type f \( -name "*_R1.fastq.gz" -o -name "*_1.fastq.gz" -o -name "*_R1.fq.gz" -o -name "*_1.fq.gz" \) -print0)

# Check if any files were found
if [ ${#r1_files[@]} -eq 0 ]; then
    log_error "No FASTQ files found in directory: $raw_reads_dir"
    log_message "Checking subdirectories..."
    
    # Search in subdirectories
    while IFS= read -r -d '' dir; do
        log_message "Examining subdirectory: $dir"
        subdir_r1_files=()
        while IFS= read -r -d '' file; do
            subdir_r1_files+=("$file")
            log_message "Found file in subdirectory: $file"
        done < <(find "$dir" -type f \( -name "*_R1.fastq.gz" -o -name "*_1.fastq.gz" -o -name "*_R1.fq.gz" -o -name "*_1.fq.gz" \) -print0)
        
        if [ ${#subdir_r1_files[@]} -gt 0 ]; then
            r1_files=("${subdir_r1_files[@]}")
            raw_reads_dir="$dir"
            break
        fi
    done < <(find "$raw_reads_dir" -mindepth 1 -maxdepth 1 -type d -print0)
fi

# Final check
if [ ${#r1_files[@]} -eq 0 ]; then
    log_error "No FASTQ files found in directory or its subdirectories: $raw_reads_dir"
    exit 1
fi

log_message "Found ${#r1_files[@]} R1 files to process"

# --- Step 5: Define Fastp Function for Parallel Execution ---
run_fastp() {
    local r1_file="$1"
    local sample_name=$(basename "$r1_file" | sed -E 's/(_R1)?(_1)?\.(fastq|fq)\.gz//')
    
    # Find corresponding R2 file
    local r2_file=""
    for pattern in "_R2.fastq.gz" "_2.fastq.gz" "_R2.fq.gz" "_2.fq.gz"; do
        potential_r2="${r1_file/_R1.fastq.gz/$pattern}"
        potential_r2="${potential_r2/_1.fastq.gz/$pattern}"
        potential_r2="${potential_r2/_R1.fq.gz/$pattern}"
        potential_r2="${potential_r2/_1.fq.gz/$pattern}"
        
        if [ -f "$potential_r2" ]; then
            r2_file="$potential_r2"
            break
        fi
    done
    
    # Validate R2 file
    if [ -z "$r2_file" ]; then
        log_error "No R2 file found for sample: $sample_name"
        return 1
    fi
    
    local sample_log="${sample_logs_dir}/${sample_name}.log"
    
    # Define output files
    local output_r1="${trimmed_reads_dir}/${sample_name}_R1_T.fastq.gz"
    local output_r2="${trimmed_reads_dir}/${sample_name}_R2_T.fastq.gz"
    local json_report="${summaries_dir}/${sample_name}_fastp.json"
    local html_report="${summaries_dir}/${sample_name}_fastp.html"
    
    echo "Processing: $sample_name" > "$sample_log"
    echo "R1 file: $r1_file" >> "$sample_log"
    echo "R2 file: $r2_file" >> "$sample_log"
    
    # Run fastp with optimal settings
    fastp \
        -i "$r1_file" -I "$r2_file" \
        -o "$output_r1" -O "$output_r2" \
        -w "$threads_per_job" \
        --json "$json_report" --html "$html_report" \
        --detect_adapter_for_pe \
        --cut_front --cut_tail \
        --cut_window_size 4 --cut_mean_quality 20 \
        --length_required 36 \
        --correction \
        --verbose >> "$sample_log" 2>&1
    
    local exit_code=$?
    if [ $exit_code -ne 0 ]; then
        echo "fastp failed for sample: $sample_name (Exit code: $exit_code)" >> "$error_log"
        return 1
    else
        echo "fastp successfully processed sample: $sample_name" >> "$progress_log"
        return 0
    fi
}

export -f run_fastp
export trimmed_reads_dir summaries_dir threads_per_job error_log progress_log sample_logs_dir

# --- Step 6: Process Files in Parallel ---
log_message "Starting parallel fastp processing with $max_parallel_jobs jobs"
parallel -u --jobs "$max_parallel_jobs" --no-notice --joblog "$parallel_joblog" run_fastp ::: "${r1_files[@]}"

# --- Step 7: Validate Trimmed Files ---
log_message "Validating trimmed output files"
total_samples=${#r1_files[@]}
successful_samples=0
failed_samples=0

for r1_file in "${r1_files[@]}"; do
    # Extract sample name
    sample_name=$(basename "$r1_file" | sed -E 's/(_R1)?(_1)?\.(fastq|fq)\.gz//')
    
    # Check output files
    r1_out="${trimmed_reads_dir}/${sample_name}_R1_T.fastq.gz"
    r2_out="${trimmed_reads_dir}/${sample_name}_R2_T.fastq.gz"
    
    if [ -f "$r1_out" ] && [ -f "$r2_out" ]; then
        # Check if files are valid gzip files
        gzip -t "$r1_out" && gzip -t "$r2_out"
        if [ $? -eq 0 ]; then
            successful_samples=$((successful_samples + 1))
        else
            log_error "Corrupted output files for sample: $sample_name"
            failed_samples=$((failed_samples + 1))
        fi
    else
        log_error "Missing output files for sample: $sample_name"
        failed_samples=$((failed_samples + 1))
    fi
done

# Record end time and calculate duration
end_time=$(date +%s)
duration=$((end_time - start_time))
duration_formatted=$(printf '%dh:%dm:%ds' $((duration/3600)) $((duration%3600/60)) $((duration%60)))

# --- Step 8: Generate Summary Report ---
log_message "Generating summary report"
failed_jobs=$(grep -c 'exit code: [^0]' "$parallel_joblog" 2>/dev/null || echo "0")

summary_file="${logs_dir}/05.fastp_quality_trim_summary.txt"
{
    echo "===== Fastp Quality Trimming Summary ====="
    echo "Date: $(date)"
    echo "Execution time: $duration_formatted"
    echo
    echo "Sample Processing:"
    echo "  Total samples: $total_samples"
    echo "  Successfully processed: $successful_samples"
    echo "  Failed: $failed_samples"
    echo "  Failed jobs (from parallel): $failed_jobs"
    echo
    echo "Output Directories:"
    echo "  Trimmed reads: $trimmed_reads_dir"
    echo "  Summary reports: $summaries_dir"
    echo
    if [ "$failed_samples" -gt 0 ]; then
        echo "Failed samples:"
        grep 'fastp failed for sample' "$error_log" | sed 's/.*fastp failed for sample: /  - /'
        echo
    fi
    echo "========================================"
} > "$summary_file"

# --- Step 9: Deactivate Conda Environment ---
log_message "Deactivating conda environment"
conda deactivate

# --- Step 10: Log Completion ---
if [ "$failed_samples" -gt 0 ]; then
    log_error "Fastp quality trimming completed with $failed_samples errors. See $error_log for details."
    cat "$summary_file"
    exit 1
else
    log_message "Fastp quality trimming completed successfully for all $total_samples samples"
    cat "$summary_file"
    exit 0
fi
