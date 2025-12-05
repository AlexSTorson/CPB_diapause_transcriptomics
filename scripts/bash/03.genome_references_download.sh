#!/bin/bash

#SBATCH --time=1:00:00                  # Walltime limit (HH:MM:SS)
#SBATCH --nodes=1                        # Number of nodes
#SBATCH --ntasks-per-node=1             # Single thread
#SBATCH --partition=short                # Partition type
#SBATCH --job-name="genome_ref_download" # Descriptive job name
#SBATCH --output="03.genome_references_download_logs/03.genome_references_download_%j.out"
#SBATCH --error="03.genome_references_download_logs/03.genome_references_download_%j.err"

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
output_dir="$REFERENCE_DIR"
ftp_base="$GENOME_FTP_BASE"

# Standard log setup
scripts_dir="$(pwd)"
logs_dir="${scripts_dir}/03.genome_references_download_logs"

# Define log files with consistent naming convention
error_log="${logs_dir}/03.genome_references_download_errors.log"
progress_log="${logs_dir}/03.genome_references_download_progress.log"
download_log="${logs_dir}/03.genome_references_download_downloads.log"
wget_log="${logs_dir}/03.genome_references_download_wget.log"  # New log for wget output

# Specific files to download
refseq_genome="$GENOME_FTP_FASTA"
gtf_file="$GENOME_FTP_GTF"
feature_file="$GENOME_FTP_FEATURE"

# Clear logs at the start of the script
> "$error_log"
> "$progress_log"
> "$download_log"
> "$wget_log"

# --- Logging Functions ---
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$progress_log"
}

log_error() {
    echo "[ERROR] [$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$error_log"
}

# Record start time for performance tracking
start_time=$(date +%s)

# --- Step 1: Create Output Directory ---
log_message "Creating output directory: $output_dir"
mkdir -p "$output_dir"
if [ $? -ne 0 ]; then
    log_error "Failed to create output directory: $output_dir"
    exit 1
fi

# --- Step 2: Download Files with Retry Logic ---
download_with_retry() {
    local url="$1"
    local output_dir="$2"
    local max_retries=3
    local retry_count=0
    local success=false
    local filename=$(basename "$url")
    local output_file="${output_dir}/${filename}"
    
    # Check if file already exists
    if [ -f "$output_file" ]; then
        log_message "File already exists: $output_file. Will be overwritten."
        rm -f "$output_file"
    fi
    
    log_message "Starting download of $filename"
    
    while [ $retry_count -lt $max_retries ] && [ "$success" = "false" ]; do
        # Redirect wget's stderr to the wget log instead of the error log
        wget --progress=dot:giga -P "$output_dir" "$url" >> "$wget_log" 2>&1
        
        if [ $? -eq 0 ]; then
            success=true
            log_message "Successfully downloaded $filename"
            echo "$(date): $filename - Download successful" >> "$download_log"
        else
            retry_count=$((retry_count + 1))
            log_error "Attempt $retry_count failed. Retrying download of $filename..."
            # Clean up potentially partial/corrupted file before retry
            rm -f "${output_dir}/${filename}"
            sleep 5
        fi
    done
    
    if [ "$success" = "false" ]; then
        log_error "Failed to download $filename after $max_retries attempts"
        echo "$(date): $filename - Download FAILED after $max_retries attempts" >> "$download_log"
        return 1
    fi
    
    return 0
}

# --- Step 3: Download Files ---
log_message "Starting genome reference downloads"

download_failed=0
for file_url in "$refseq_genome" "$gtf_file" "$feature_file"; do
    filename=$(basename "$file_url")
    target_file="${output_dir}/${filename}"
    
    # Check if file already exists
    if [ -f "$target_file" ]; then
        log_message "File already exists: $target_file. Will be overwritten."
        rm -f "$target_file"  # Remove existing file to ensure clean download
    fi
    
    download_with_retry "$file_url" "$output_dir"
    if [ $? -ne 0 ]; then
        log_error "Download process failed for: $file_url"
        download_failed=1
    fi
done

if [ $download_failed -eq 1 ]; then
    log_error "One or more downloads failed. See $error_log for details."
    exit 1
fi

# --- Step 3.5: Unzip Downloaded Files ---
log_message "Unzipping downloaded genome files"

# Function to unzip a file with error handling
unzip_file() {
    local compressed_file="$1"
    local output_file="${compressed_file%.gz}"
    
    # Check if the output file already exists
    if [ -f "$output_file" ]; then
        log_message "Uncompressed file already exists: $output_file. Will be overwritten."
        rm -f "$output_file"
    fi
    
    log_message "Unzipping: $compressed_file"
    gunzip -c "$compressed_file" > "$output_file"
    
    if [ $? -ne 0 ]; then
        log_error "Failed to unzip file: $compressed_file"
        return 1
    else
        log_message "Successfully unzipped to: $output_file"
        # Remove the compressed file to save space (optional)
        # rm "$compressed_file"
        return 0
    fi
}

# List of files to unzip
compressed_files=(
    "${output_dir}/${GENOME_FASTA_FILENAME}.gz"
    "${output_dir}/${GENOME_GTF_FILENAME}.gz"
    "${output_dir}/${GENOME_FEATURE_FILENAME}.gz"
)

# Unzip each file
unzip_failures=0
for file in "${compressed_files[@]}"; do
    if [ -f "$file" ]; then
        unzip_file "$file"
        if [ $? -ne 0 ]; then
            unzip_failures=$((unzip_failures + 1))
        fi
    else
        log_error "Compressed file not found for unzipping: $file"
        unzip_failures=$((unzip_failures + 1))
    fi
done

if [ $unzip_failures -gt 0 ]; then
    log_error "Some files failed to unzip: $unzip_failures failures"
    # You can choose to exit here if unzipping is critical
    # exit 1
else
    log_message "All files successfully unzipped"
fi

# --- Step 4: Verify Downloads ---
log_message "Verifying downloaded files"

# Check for uncompressed files (after unzipping)
expected_files=(
    "${output_dir}/${GENOME_FASTA_FILENAME}"
    "${output_dir}/${GENOME_GTF_FILENAME}"
    "${output_dir}/${GENOME_FEATURE_FILENAME}"
)

missing_files=0
for file in "${expected_files[@]}"; do
    if [ ! -f "$file" ]; then
        log_error "Missing expected file: $file"
        missing_files=$((missing_files + 1))
    else
        # Get file size for verification
        file_size=$(stat -c%s "$file")
        file_size_mb=$(echo "scale=2; $file_size / 1048576" | bc)
        log_message "Verified file: $(basename "$file") (${file_size_mb} MB)"
    fi
done

if [ $missing_files -gt 0 ]; then
    log_error "Download verification failed: $missing_files files are missing"
    exit 1
fi

# Record end time and calculate duration
end_time=$(date +%s)
duration=$((end_time - start_time))
duration_formatted=$(printf '%dh:%dm:%ds' $((duration/3600)) $((duration%3600/60)) $((duration%60)))

# --- Step 5: Generate Summary Report ---
summary_file="${logs_dir}/03.genome_references_download_summary.txt"
{
    echo "===== Genome References Download Summary ====="
    echo "Date: $(date)"
    echo "Execution time: $duration_formatted"
    echo
    echo "Downloaded Files:"
    for file in "${expected_files[@]}"; do
        if [ -f "$file" ]; then
            file_size=$(stat -c%s "$file")
            file_size_mb=$(echo "scale=2; $file_size / 1048576" | bc)
            echo "  - $(basename "$file"): ${file_size_mb} MB"
        fi
    done
    echo
    echo "Files downloaded to: $output_dir"
    echo "========================================"
} > "$summary_file"

# --- Step 6: Log Completion ---
if [ -s "$error_log" ]; then
    # Check if there are real errors or just wget output
    real_errors=$(grep -c "\[ERROR\]" "$error_log")
    if [ "$real_errors" -gt 0 ]; then
        log_error "Genome references download completed with errors. Check $error_log for details."
        cat "$summary_file"
        exit 1
    else
        log_message "Genome references download completed successfully in $duration_formatted"
        cat "$summary_file"
        exit 0
    fi
else
    log_message "Genome references download completed successfully in $duration_formatted"
    cat "$summary_file"
    exit 0
fi
