#!/bin/bash

#SBATCH --job-name="genome_index"       # Job name for identification
#SBATCH --output="08.genome_index_logs/08.genome_index_%j.out"
#SBATCH --error="08.genome_index_logs/08.genome_index_%j.err"
#SBATCH --time=24:00:00                 # Reduced walltime based on expected performance
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks-per-node=40            # Number of tasks per node
#SBATCH --partition=short               # Partition to use
#SBATCH --mem=0                         # Request all available memory

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
scripts_dir="$PROJECT_DIR/scripts"
genome_references_dir="$REFERENCE_DIR"
genome_index_dir="$GENOME_INDEX_DIR"
ss_file="$SPLICE_SITES_FILE"
exon_file="$EXONS_FILE"
genome_file_base="${REFERENCE_DIR}/${GENOME_FASTA_FILENAME%.*}"  # Remove extension
output_prefix="$GENOME_INDEX"

# Standard log setup
current_dir="$(pwd)"
logs_dir="${current_dir}/08.genome_index_logs"

# Define log files with consistent naming convention
error_log="${logs_dir}/08.genome_index_errors.log"
progress_log="${logs_dir}/08.genome_index_progress.log"
hisat2_log="${logs_dir}/08.genome_index_hisat2.log"

# Performance tuning
hisat_threads="$THREADS_PER_HISAT"  # Use threads from config

# Clear logs at the start of the script
> "$error_log"
> "$progress_log"
> "$hisat2_log"

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
$MODULE_HISAT2

# --- Step 2: Create Required Directories ---
log_message "Creating output directory: $genome_index_dir"
mkdir -p "$genome_index_dir"
if [ $? -ne 0 ]; then
    log_error "Failed to create output directory: $genome_index_dir"
    exit 1
fi

# --- Step 3: Handle Compressed or Uncompressed Genome Files ---
log_message "Checking genome file"
if [[ -f "${genome_file_base}.fna.gz" ]]; then
    log_message "Compressed genome file detected (${genome_file_base}.fna.gz). Unzipping..."
    
    # Check if uncompressed file already exists
    if [[ -f "${genome_file_base}.fna" ]]; then
        log_message "Uncompressed genome file already exists: ${genome_file_base}.fna"
    else
        log_message "Uncompressing genome file..."
        gunzip -c "${genome_file_base}.fna.gz" > "${genome_file_base}.fna"
        if [ $? -ne 0 ]; then
            log_error "Failed to unzip genome file: ${genome_file_base}.fna.gz"
            exit 1
        else
            log_message "Unzipped genome file: ${genome_file_base}.fna"
        fi
    fi
    genome_file="${genome_file_base}.fna"  # Use the uncompressed file
elif [[ -f "${genome_file_base}.fna" ]]; then
    log_message "Uncompressed genome file detected: ${genome_file_base}.fna"
    genome_file="${genome_file_base}.fna"
else
    log_error "Genome file not found in either compressed (.gz) or uncompressed format: ${genome_file_base}"
    exit 1
fi

# --- Step 4: Validate Input Files ---
log_message "Validating input files"

if [ ! -f "$ss_file" ]; then
    log_error "Splice sites file not found: $ss_file"
    exit 1
fi

if [ ! -f "$exon_file" ]; then
    log_error "Exons file not found: $exon_file"
    exit 1
fi

if [ ! -f "$genome_file" ]; then
    log_error "Genome file not found: $genome_file"
    exit 1
fi

# Get file sizes and counts for reporting
ss_count=$(wc -l < "$ss_file")
exon_count=$(wc -l < "$exon_file")
genome_size=$(stat -c%s "$genome_file")
genome_size_gb=$(echo "scale=2; $genome_size / 1073741824" | bc)

log_message "Genome file size: $genome_size_gb GB"
log_message "Splice sites file contains $ss_count entries"
log_message "Exons file contains $exon_count entries"

# --- Step 5: Estimate Memory Requirements ---
# HISAT2 needs approximately 9-10 GB per 1 billion bases
genome_bases=$(grep -v ">" "$genome_file" | tr -d '\n' | wc -c)
est_memory_gb=$(echo "scale=1; ${genome_bases}/100000000" | bc)
log_message "Estimated genome size: $(numfmt --to=iec-i --suffix=B $genome_bases) bases"
log_message "Estimated memory requirement: approximately ${est_memory_gb} GB"

# --- Step 6: Run HISAT2 Genome Indexing ---
log_message "Starting HISAT2 genome indexing with ${hisat_threads} threads"
log_message "Command: hisat2-build --ss \"$ss_file\" --exon \"$exon_file\" -p ${hisat_threads} \"$genome_file\" \"$output_prefix\""

# Run the indexing
hisat2-build --ss "$ss_file" --exon "$exon_file" -p ${hisat_threads} "$genome_file" "$output_prefix" 2> "$hisat2_log"

# Check if indexing was successful
if [ $? -ne 0 ]; then
    log_error "HISAT2 genome indexing failed. See $hisat2_log for details."
    exit 1
else
    log_message "HISAT2 genome indexing completed successfully."
fi

# --- Step 7: Verify Index Files ---
log_message "Verifying index files"
index_files=("${output_prefix}."*.ht2)
if [ ${#index_files[@]} -lt 8 ]; then
    log_error "Expected 8 index files but found ${#index_files[@]}"
    exit 1
else
    log_message "Found ${#index_files[@]} index files"
    # Log the size of each index file
    total_size=0
    for file in "${index_files[@]}"; do
        file_size=$(stat -c%s "$file")
        total_size=$((total_size + file_size))
        log_message "  $(basename "$file"): $(numfmt --to=iec-i --suffix=B $file_size)"
    done
    log_message "Total index size: $(numfmt --to=iec-i --suffix=B $total_size)"
fi

# Record end time and calculate duration
end_time=$(date +%s)
duration=$((end_time - start_time))
duration_formatted=$(printf '%dh:%dm:%ds' $((duration/3600)) $((duration%3600/60)) $((duration%60)))

# --- Step 8: Generate Summary Report ---
summary_file="${logs_dir}/08.genome_index_summary.txt"
{
    echo "===== HISAT2 Genome Indexing Summary ====="
    echo "Date: $(date)"
    echo "Execution time: $duration_formatted"
    echo
    echo "Input Files:"
    echo "  Genome file: $genome_file"
    echo "  Genome size: $genome_size_gb GB"
    echo "  Splice sites file: $ss_file with $ss_count sites"
    echo "  Exons file: $exon_file with $exon_count exons"
    echo
    echo "Output:"
    echo "  Index prefix: $output_prefix"
    echo "  Index files created: ${#index_files[@]}"
    echo "  Total index size: $(numfmt --to=iec-i --suffix=B $total_size)"
    echo
    if [ -s "$hisat2_log" ]; then
        echo "HISAT2 Summary (from log file):"
        grep "Total time for call to driver" "$hisat2_log" 2>/dev/null || echo "  No timing information found"
        grep "Memory usage" "$hisat2_log" 2>/dev/null || echo "  No memory usage information found"
    fi
    echo "=================================================="
} > "$summary_file"

# --- Step 9: Log Completion ---
if [ -s "$error_log" ]; then
    log_error "Genome indexing completed with errors. Check $error_log for details."
    cat "$summary_file"
    exit 1
else
    log_message "Genome indexing completed successfully in $duration_formatted"
    cat "$summary_file"
    exit 0
fi
