#!/bin/bash

#SBATCH --time=01:00:00                 # Walltime limit (HH:MM:SS)
#SBATCH --nodes=1                       # Single node
#SBATCH --ntasks=1                      # Single task
#SBATCH --partition=short               # Partition type
#SBATCH --job-name="create_mergelist"   # Descriptive job name
#SBATCH --output="12.create_mergelist_logs/12.create_mergelist_%j.out"
#SBATCH --error="12.create_mergelist_logs/12.create_mergelist_%j.err"

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
assemblies_dir="${base_dir}/assemblies"      # Directory containing GTF files
mergelist_file="${assemblies_dir}/mergelist.txt"   # Path to output mergelist file

# Define the taco assembly merge directory to exclude
taco_assembly_dir="${base_dir}/assemblies/taco_assembly_merge"

# Standard log setup
current_dir="$(pwd)"
logs_dir="${current_dir}/12.create_mergelist_logs"

# Define log files with consistent naming convention
error_log="${logs_dir}/12.create_mergelist_errors.log"
progress_log="${logs_dir}/12.create_mergelist_progress.log"

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

# --- Step 1: Validate Assemblies Directory ---
log_message "Validating assemblies directory: $assemblies_dir"

if [ ! -d "$assemblies_dir" ]; then
    log_error "Assemblies directory not found: $assemblies_dir"
    exit 1
fi

# Check for taco output directory and warn if it exists
if [ -d "$taco_assembly_dir" ]; then
    log_message "WARNING: taco_assembly_merge directory already exists. This directory will be excluded from mergelist."
    # Optionally, you could remove it here if appropriate:
    # log_message "Removing existing taco_assembly_merge directory"
    # rm -rf "$taco_assembly_dir"
fi

# --- Step 2: Find GTF Files ---
log_message "Locating GTF files"

gtf_files=()
while IFS= read -r -d '' file; do
    # Skip any potential merged GTF from previous runs, the mergelist file itself,
    # and any files in the taco_assembly_merge directory
    if [[ "$(basename "$file")" != "stringtie_merged.gtf" && 
          "$(basename "$file")" != "mergelist.txt" && 
          "$file" != "$taco_assembly_dir"* ]]; then
        gtf_files+=("$file")
    fi
done < <(find "$assemblies_dir" -type f -name "*.gtf" -print0)

# Check if we found any GTF files
if [ ${#gtf_files[@]} -eq 0 ]; then
    log_error "No GTF files found in directory: $assemblies_dir"
    exit 1
fi

log_message "Found ${#gtf_files[@]} GTF files to include in mergelist"

# --- Step 3: Validate GTF Files ---
log_message "Validating GTF files"
valid_gtfs=0
invalid_gtfs=0
total_file_size=0
list_of_invalid_files=""

for gtf_file in "${gtf_files[@]}"; do
    if [ ! -s "$gtf_file" ]; then
        log_error "Empty GTF file: $gtf_file"
        invalid_gtfs=$((invalid_gtfs + 1))
        list_of_invalid_files="${list_of_invalid_files}$(basename "$gtf_file")\n"
    else
        # Check if file contains valid GTF format (at least has transcript_id entries)
        if ! grep -q "transcript_id" "$gtf_file"; then
            log_error "Invalid GTF format in file: $gtf_file"
            invalid_gtfs=$((invalid_gtfs + 1))
            list_of_invalid_files="${list_of_invalid_files}$(basename "$gtf_file")\n"
        else
            valid_gtfs=$((valid_gtfs + 1))
            # Add file size to total
            file_size=$(stat -c%s "$gtf_file")
            total_file_size=$((total_file_size + file_size))
        fi
    fi
done

if [ $valid_gtfs -eq 0 ]; then
    log_error "No valid GTF files found to create mergelist"
    exit 1
fi

log_message "Validated $valid_gtfs GTF files (skipped $invalid_gtfs invalid files)"
total_file_size_mb=$(echo "scale=2; $total_file_size / 1048576" | bc)
log_message "Total size of valid GTF files: ${total_file_size_mb} MB"

# --- Step 4: Create Mergelist File ---
log_message "Creating mergelist file: $mergelist_file"

# Clear file if it exists
> "$mergelist_file"  # Clear the file if it exists

# Sort files by sample name for consistency and only include valid GTFs
for gtf_file in "${gtf_files[@]}"; do
    # Only add valid GTFs to the mergelist
    if [ -s "$gtf_file" ] && grep -q "transcript_id" "$gtf_file"; then
        echo "$gtf_file" >> "$mergelist_file"
    fi
done

# Sort the mergelist file
sort "$mergelist_file" -o "$mergelist_file"

if [ $? -ne 0 ]; then
    log_error "Failed to create mergelist file: $mergelist_file"
    exit 1
fi

# Verify the mergelist file was created successfully
if [ ! -s "$mergelist_file" ]; then
    log_error "Mergelist file is empty: $mergelist_file"
    exit 1
fi

lines_in_mergelist=$(wc -l < "$mergelist_file")
log_message "Mergelist file contains $lines_in_mergelist GTF files"

# Record end time and calculate duration
end_time=$(date +%s)
duration=$((end_time - start_time))
duration_formatted=$(printf '%dm:%ds' $((duration/60)) $((duration%60)))

# --- Step 5: Generate Summary Report ---
summary_file="${logs_dir}/12.create_mergelist_summary.txt"
{
    echo "===== Create Mergelist Summary ====="
    echo "Date: $(date)"
    echo "Execution time: $duration_formatted"
    echo
    echo "GTF File Statistics:"
    echo "  GTFs found: ${#gtf_files[@]}"
    echo "  Valid GTF files: $valid_gtfs"
    echo "  Invalid GTF files: $invalid_gtfs"
    echo "  GTF files in mergelist: $lines_in_mergelist"
    echo "  Total size of valid GTF files: ${total_file_size_mb} MB"
    echo
    echo "Output:"
    echo "  Mergelist file: $mergelist_file"
    echo
    if [ $invalid_gtfs -gt 0 ]; then
        echo "Invalid GTF files:"
        echo -e "$list_of_invalid_files"
    fi
    echo "======================================="
} > "$summary_file"

# --- Step 6: Log Completion ---
if [ $invalid_gtfs -gt 0 ]; then
    log_error "Mergelist created with warnings. Skipped $invalid_gtfs invalid GTF files."
    cat "$summary_file"
    # Don't exit with error as the mergelist was still created successfully
    exit 0
else
    log_message "Mergelist created successfully with $lines_in_mergelist GTF files in $duration_formatted"
    cat "$summary_file"
    exit 0
fi
