#!/bin/bash

#SBATCH --time=1:00:00                  # Walltime limit (HH:MM:SS)
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks-per-node=8             # Reduced threads - this is mostly I/O bound
#SBATCH --partition=short               # Partition type (e.g., short, normal, long)
#SBATCH --job-name="exon_splice_extract" # Job name
#SBATCH --output="07.exon_and_splice_sites_logs/07.exon_and_splice_sites_%j.out"
#SBATCH --error="07.exon_and_splice_sites_logs/07.exon_and_splice_sites_%j.err"

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
base_dir="$REFERENCE_DIR"
gtf_file="$REFERENCE_GTF"
output_dir="$REFERENCE_DIR"

# Standard log setup
current_dir="$(pwd)"
logs_dir="${current_dir}/07.exon_and_splice_sites_logs"

# Define log files with consistent naming convention
error_log="${logs_dir}/07.exon_and_splice_sites_errors.log"
progress_log="${logs_dir}/07.exon_and_splice_sites_progress.log"

# File paths for extracted splice sites and exons
splice_sites_file="${SPLICE_SITES_FILE}"
exons_file="${EXONS_FILE}"

# Create log directory
mkdir -p "$logs_dir"

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
$MODULE_PYTHON

# --- Step 2: Check and Unzip GTF File if Needed ---
if [[ "$gtf_file" == *.gz ]]; then
    log_message "Detected compressed GTF file. Unzipping..."
    
    # Create uncompressed filename
    uncompressed_file="${gtf_file%.gz}"
    
    # Check if file already exists (from previous run)
    if [ -f "$uncompressed_file" ]; then
        log_message "Using existing uncompressed file: $uncompressed_file"
    else
        gunzip -c "$gtf_file" > "$uncompressed_file"
        if [ $? -ne 0 ]; then
            log_error "Failed to unzip GTF file: $gtf_file"
            exit 1
        else
            log_message "Unzipping completed: $uncompressed_file"
        fi
    fi
    
    gtf_file="$uncompressed_file"  # Update the GTF file path to the uncompressed version
fi

# --- Step 3: Validate Input GTF File ---
if [ ! -f "$gtf_file" ]; then
    log_error "Input GTF file not found: $gtf_file"
    exit 1
fi

log_message "Using GTF file: $gtf_file"

# Get GTF file size for logging
gtf_size=$(stat -c%s "$gtf_file")
gtf_size_mb=$(echo "scale=2; $gtf_size / 1048576" | bc)
log_message "GTF file size: ${gtf_size_mb} MB"

# --- Step 4: Check Required Python Scripts ---
scripts_dir="${PROJECT_DIR}/scripts"
hisat2_extract_splice_sites="${scripts_dir}/hisat2_extract_splice_sites.py"
hisat2_extract_exons="${scripts_dir}/hisat2_extract_exons.py"

log_message "Checking Python scripts in directory: $scripts_dir"
log_message "Splice sites script path: $hisat2_extract_splice_sites"
log_message "Exons script path: $hisat2_extract_exons"

for script in "$hisat2_extract_splice_sites" "$hisat2_extract_exons"; do
    if [ ! -f "$script" ]; then
        log_error "Required Python script not found: $script"
        # List contents of scripts directory for debugging
        log_message "Contents of scripts directory:"
        ls -l "$scripts_dir"
        exit 1
    fi
done

log_message "Found required Python scripts"

# --- Step 5: Extract Splice Sites ---
log_message "Extracting splice sites from GTF file..."
python "$hisat2_extract_splice_sites" "$gtf_file" > "$splice_sites_file"
if [ $? -ne 0 ]; then
    log_error "Failed to extract splice sites from GTF file: $gtf_file"
    exit 1
else
    # Count the number of splice sites extracted
    splice_count=$(wc -l < "$splice_sites_file")
    splice_size=$(stat -c%s "$splice_sites_file")
    splice_size_kb=$(echo "scale=2; $splice_size / 1024" | bc)
    log_message "Successfully extracted $splice_count splice sites to: $splice_sites_file (${splice_size_kb} KB)"
fi

# --- Step 6: Extract Exons ---
log_message "Extracting exons from GTF file..."
python "$hisat2_extract_exons" "$gtf_file" > "$exons_file"
if [ $? -ne 0 ]; then
    log_error "Failed to extract exons from GTF file: $gtf_file"
    exit 1
else
    # Count the number of exons extracted
    exon_count=$(wc -l < "$exons_file")
    exon_size=$(stat -c%s "$exons_file")
    exon_size_kb=$(echo "scale=2; $exon_size / 1024" | bc)
    log_message "Successfully extracted $exon_count exons to: $exons_file (${exon_size_kb} KB)"
fi

# --- Step 7: Validate Output Files ---
validation_errors=0

if [ ! -s "$splice_sites_file" ]; then
    log_error "Splice sites file is empty: $splice_sites_file"
    validation_errors=$((validation_errors + 1))
fi

if [ ! -s "$exons_file" ]; then
    log_error "Exons file is empty: $exons_file"
    validation_errors=$((validation_errors + 1))
fi

if [ $validation_errors -gt 0 ]; then
    log_error "Validation found $validation_errors errors with output files"
else
    log_message "Output file validation successful"
fi

# Record end time and calculate duration
end_time=$(date +%s)
duration=$((end_time - start_time))
duration_formatted=$(printf '%dm:%ds' $((duration/60)) $((duration%60)))

# --- Step 8: Generate Summary Report ---
summary_file="${logs_dir}/07.exon_and_splice_sites_summary.txt"
{
    echo "===== Exon and Splice Site Extraction Summary ====="
    echo "Date: $(date)"
    echo "Execution time: $duration_formatted"
    echo
    echo "Input:"
    echo "  GTF file: $gtf_file"
    echo "  GTF file size: ${gtf_size_mb} MB"
    echo
    echo "Output:"
    echo "  Splice sites file: $splice_sites_file"
    echo "  Number of splice sites: $splice_count"
    echo "  Splice sites file size: ${splice_size_kb} KB"
    echo
    echo "  Exons file: $exons_file"
    echo "  Number of exons: $exon_count"
    echo "  Exons file size: ${exon_size_kb} KB"
    echo
    if [ $validation_errors -gt 0 ]; then
        echo "Validation Errors: $validation_errors"
        echo "Check $error_log for details"
    else
        echo "Validation: All output files passed validation"
    fi
    echo "=================================================="
} > "$summary_file"

# --- Step 9: Log Completion ---
if [ -s "$error_log" ]; then
    log_error "Exon and splice site extraction completed with errors. Check $error_log for details."
    cat "$summary_file"
    exit 1
else
    log_message "Exon and splice site extraction completed successfully in $duration_formatted"
    cat "$summary_file"
    exit 0
fi
