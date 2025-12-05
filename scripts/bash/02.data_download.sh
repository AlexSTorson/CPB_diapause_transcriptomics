#!/bin/bash

#SBATCH --time=12:00:00                  # Increased walltime to account for rate limiting
#SBATCH --nodes=1                        # Number of nodes
#SBATCH --ntasks-per-node=12             # Reduced number of tasks for better API management
#SBATCH --mem=64G                        # Memory request
#SBATCH --partition=short                # Partition type
#SBATCH --job-name="data_download"       # Descriptive job name
#SBATCH --output="02.data_download_logs/02.data_download_%j.out"
#SBATCH --error="02.data_download_logs/02.data_download_%j.err"

# --- Load the Pipeline Configuration ---
# --- Load the Pipeline Configuration with absolute path ---
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



#!/bin/bash

#!/bin/bash
# SLURM job configuration


# --- Centralized Variables ---
base_dir="$PROJECT_DIR"
raw_reads_dir="${base_dir}/raw_reads"
sra_dir="${base_dir}/sra_files"
sra_cache="${base_dir}/sra_cache"

# Create temp directory on local node storage for better I/O performance
temp_dir="${TMPDIR:-/tmp}/sra_download_${SLURM_JOB_ID}"

# Save logs in a dedicated log directory with consistent naming
scripts_dir="$(pwd)"  # Use the current working directory
logs_dir="${scripts_dir}/02.data_download_logs"

# Define log files with consistent naming convention
error_log="${logs_dir}/02.data_download_errors.log"
progress_log="${logs_dir}/02.data_download_progress.log"
parallel_log="${logs_dir}/02.data_download_parallel.log"
rename_log="${logs_dir}/02.data_download_rename.log"
validation_log="${logs_dir}/02.data_download_validation.log"
eutils_log="${logs_dir}/02.data_download_eutils.log"

mapping_file="${sra_dir}/sample_srr_list.csv"
srr_list="${sra_dir}/SRR_accessions.txt"
jobs=12  # Reduced number of parallel jobs

# BioProject ID - Replace with your actual BioProject ID
bioproject="$BIOPROJECT"

# E-utilities API key - Set as an environment variable
# Replace with your actual API key or remove if not needed
export NCBI_API_KEY="$NCBI_API_KEY"

# Create a user-specific configuration directory in your project space
export NCBI_HOME="${scripts_dir}/02.data_download_ncbi_home"
mkdir -p "$NCBI_HOME"

# Tell SRA toolkit to use this directory for user settings and cache
export HOME="$NCBI_HOME"  # Override HOME for this script only
export NCBI_SETTINGS="${NCBI_HOME}/.ncbi"
mkdir -p "$NCBI_SETTINGS"

# Rate limiting parameters
SLEEP_TIME=1   # Time between E-utilities API calls in seconds
MAX_RETRIES=5  # Maximum number of retries for E-utilities API calls

# Clear logs at the start of the script
> "$error_log"
> "$progress_log"
> "$eutils_log"

# --- Step 1: Load Necessary Modules ---
$MODULE_EDIRECT                  # Tools for accessing NCBI databases
$MODULE_SRA         # SRA Toolkit for downloading and processing SRA files
$MODULE_PARALLEL                  # GNU Parallel for multi-threaded execution
$MODULE_PIGZ                     # Parallel gzip compression

# --- Step 2: Create Output Directories (if they do not exist) ---
for dir in "$raw_reads_dir" "$sra_dir" "$sra_cache" "$temp_dir" "$NCBI_SETTINGS"; do
    if [ ! -d "$dir" ]; then
        mkdir -p "$dir"
        if [ $? -ne 0 ]; then
            echo "Failed to create directory: $dir" >> "$error_log"
            exit 1  # Exit script if directory creation fails
        fi
    fi
done

# --- Helper Functions ---
# Function to make E-utilities API calls with rate limiting and retries
function call_eutils() {
    local cmd="$1"
    local retries=0
    local success=false
    
    # Log the command to the eutils log
    echo "Executing: $cmd" >> "$eutils_log"
    
    while [ $retries -lt $MAX_RETRIES ] && [ "$success" = "false" ]; do
        result=$(eval "$cmd")
        exit_code=$?
        
        if [ $exit_code -eq 0 ] && [ -n "$result" ]; then
            success=true
            echo "$result"
        else
            retries=$((retries + 1))
            echo "Attempt $retries failed. Waiting before retry..." >> "$eutils_log"
            sleep $((SLEEP_TIME * retries))  # Progressive backoff
        fi
    done
    
    if [ "$success" = "false" ]; then
        echo "Failed after $MAX_RETRIES attempts: $cmd" >> "$error_log"
        return 1
    fi
    
    # Ensure rate limiting between successful calls
    sleep $SLEEP_TIME
    return 0
}

# --- Step 3: Generate a List of All SRA Files in the BioProject ---
echo "Retrieving SRR accessions for BioProject $bioproject..."

# Use the helper function for E-utilities calls
runinfo=$(call_eutils "esearch -db sra -query \"$bioproject\"" | call_eutils "efetch -format runinfo")
if [ $? -ne 0 ]; then
    echo "Error retrieving runinfo for BioProject: $bioproject" >> "$error_log"
    exit 1
fi

# Create SRR list - explicitly filter out header rows
echo "$runinfo" | cut -d "," -f 1 | grep -v "^Run$" | grep -v "^SRR$" | grep -v "^$" > "$srr_list"
echo "Found $(wc -l < "$srr_list") SRR accessions."

# --- Step 4: Generate Mapping File ---
echo "Creating sample to SRR mapping file..."
echo "$runinfo" | cut -d "," -f 1,12 | \
awk '
    BEGIN { FS=OFS="," }
    NR==1 && $1 != "SRR" { print "SRR,LibraryName"; next }
    $1 ~ /^(SRR|Run)$/ { next }  # Skip header rows
    { print }
' > "$mapping_file"

if [ $? -ne 0 ] || [ ! -s "$mapping_file" ]; then
    echo "Error generating mapping file (SRR and Library Name)." >> "$error_log"
    exit 1
fi
echo "Mapping file successfully created at $mapping_file."

# --- Step 5: Download and Convert SRA to FASTQ Format with Better Management ---
echo "Starting SRA downloads and FASTQ conversion..."

# Create a download function with better error handling
download_and_convert() {
    local srr=$1
    
    # Skip header rows and empty values
    if [[ "$srr" == "Run" || "$srr" == "SRR" || -z "$srr" ]]; then
        echo "Skipping header or empty value: '$srr'"
        return 0
    fi
    
    # Validate SRR accession format
    if ! [[ "$srr" =~ ^SRR[0-9]+$ ]]; then
        echo "Invalid SRR accession format: '$srr' - skipping" >> "$error_log"
        return 0
    fi
    
    local retries=0
    local max_dl_retries=3
    local success=false
    local sra_file="${temp_dir}/${srr}.sra"
    
    echo "Processing $srr..."
    
    while [ $retries -lt $max_dl_retries ] && [ "$success" = "false" ]; do
        # Download directly using prefetch with explicit output file
        prefetch -O "$temp_dir" "$srr" 2> "${temp_dir}/${srr}_prefetch.log"
        
        # Check if the SRA file exists (location can vary depending on SRA toolkit version)
        if [ -f "$sra_file" ] || [ -f "${temp_dir}/${srr}/${srr}.sra" ]; then
            # Find the actual SRA file
            local actual_sra=""
            if [ -f "$sra_file" ]; then
                actual_sra="$sra_file"
            else
                actual_sra="${temp_dir}/${srr}/${srr}.sra"
            fi
            
            # Convert SRA to FASTQ
            fasterq-dump --outdir "$temp_dir" --threads 2 --split-files "$actual_sra" 2> "${temp_dir}/${srr}_fasterq.log"
            
            if [ $? -eq 0 ] && [ -f "${temp_dir}/${srr}_1.fastq" ]; then
                # Compress FASTQ files
                pigz -p 2 "${temp_dir}/${srr}"*.fastq
                
                # Move to final destination
                mv "${temp_dir}/${srr}"*.fastq.gz "$raw_reads_dir/"
                success=true
                echo "Successfully processed $srr"
            else
                retries=$((retries + 1))
                echo "Attempt $retries failed for $srr at fasterq-dump stage. Waiting before retry..." >> "$error_log"
                cat "${temp_dir}/${srr}_fasterq.log" >> "$error_log"
                sleep $((10 * retries))  # Progressive backoff
            fi
        else
            retries=$((retries + 1))
            echo "Prefetch attempt $retries failed for $srr. Waiting before retry..." >> "$error_log"
            cat "${temp_dir}/${srr}_prefetch.log" >> "$error_log"
            sleep $((10 * retries))  # Progressive backoff
        fi
    done
    
    if [ "$success" = "false" ]; then
        echo "Failed to download $srr after $max_dl_retries attempts" >> "$error_log"
        return 1
    fi
    
    # Clean up SRA files to save space
    rm -f "$sra_file"
    rm -rf "${temp_dir}/${srr}"
    
    return 0
}

export -f download_and_convert
export temp_dir raw_reads_dir error_log

# Export NCBI settings to make them available to parallel processes
export HOME NCBI_HOME NCBI_SETTINGS NCBI_API_KEY

# Use GNU parallel but with a lower job count and export all needed variables
echo "Starting parallel downloads with $jobs concurrent jobs"
parallel --joblog "$parallel_log" --jobs "$jobs" \
    download_and_convert {} :::: "$srr_list"

# Check for errors in parallel job execution
failed_jobs=0
while IFS= read -r line; do
    # Skip header line
    if [[ "$line" == "Seq"* ]]; then continue; fi
    
    # Extract exit value (8th field)
    exitval=$(echo "$line" | awk '{print $7}')
    if [ "$exitval" -ne 0 ]; then
        command=$(echo "$line" | awk '{print $NF}')
        echo "Job failed: $command (Exit code: $exitval)" >> "$error_log"
        failed_jobs=$((failed_jobs + 1))
    fi
done < "$parallel_log"

if [ "$failed_jobs" -gt 0 ]; then
    echo "$failed_jobs jobs failed during download. Check $parallel_log for details." >> "$error_log"
    # Continue execution to process successful downloads
else
    echo "All download jobs completed successfully." >> "$progress_log"
fi

# --- Step 6: Rename FASTQ Files Using Mapping File ---
if [ ! -f "$mapping_file" ]; then
    echo "Mapping file not found: $mapping_file" >> "$error_log"
    exit 1
fi

# Rename log already defined at the top of the script
> "$rename_log"

echo "Renaming files based on mapping file..."
while IFS=',' read -r srr sample_id; do
    # Skip empty lines or header
    [ -z "$srr" ] || [ "$srr" = "SRR" ] || [ "$srr" = "Run" ] && continue
    
    # Find FASTQ files corresponding to the SRR number
    r1_file="${raw_reads_dir}/${srr}_1.fastq.gz"
    r2_file="${raw_reads_dir}/${srr}_2.fastq.gz"

    # Rename FASTQ files to use sample ID and R1/R2
    if [ -f "$r1_file" ]; then
        mv "$r1_file" "${raw_reads_dir}/${sample_id}_R1.fastq.gz"
        if [ $? -ne 0 ]; then
            echo "Failed to rename file: $r1_file" >> "$error_log"
        else
            echo "Renamed $r1_file to ${raw_reads_dir}/${sample_id}_R1.fastq.gz" >> "$rename_log"
        fi
    else
        echo "Warning: Expected file $r1_file not found" >> "$error_log"
    fi
    
    if [ -f "$r2_file" ]; then
        mv "$r2_file" "${raw_reads_dir}/${sample_id}_R2.fastq.gz"
        if [ $? -ne 0 ]; then
            echo "Failed to rename file: $r2_file" >> "$error_log"
        else
            echo "Renamed $r2_file to ${raw_reads_dir}/${sample_id}_R2.fastq.gz" >> "$rename_log"
        fi
    else
        echo "Warning: Expected file $r2_file not found" >> "$error_log"
    fi
done < "$mapping_file"

# --- Step 7: Validate Downloads ---
echo "Validating downloads..."
# Validation log already defined at the top of the script
> "$validation_log"

# Count expected and actual files - Fixed to correctly count samples
expected_samples=$(grep -v "^SRR" "$mapping_file" | grep -v "^Run" | grep -v "^$" | wc -l)
expected_files=$((expected_samples * 2))  # R1 and R2 for each sample
actual_files=$(find "$raw_reads_dir" -name "*_R[12].fastq.gz" | wc -l)

echo "Expected files: $expected_files" >> "$validation_log"
echo "Found files: $actual_files" >> "$validation_log"

if [ "$actual_files" -lt "$expected_files" ]; then
    echo "WARNING: Not all expected files were downloaded and renamed." >> "$validation_log"
    echo "Check $error_log for details on failed downloads." >> "$validation_log"
else
    echo "All expected files were successfully downloaded and renamed." >> "$validation_log"
fi

# --- Step 8: Clean up temporary directory and NCBI home directories ---
echo "Cleaning up temporary directories..."
rm -rf "$temp_dir"

# Clean up NCBI home directories now that they're no longer needed
echo "Cleaning up NCBI home directories..."
rm -rf "$NCBI_HOME"
echo "NCBI home directories successfully removed."

# --- Step 9: Log Completion ---
# Only report errors if there are real download failures
real_errors=0
if [ -s "$error_log" ]; then
    # Check for real errors, excluding expected warnings
    real_errors=$(grep -v "Invalid SRR accession format" "$error_log" | grep -v "jobs failed during download" | wc -l)
fi

if [ "$real_errors" -gt 0 ] || [ "$actual_files" -lt "$expected_files" ]; then
    echo "Data download completed with errors. Check $error_log for details."
else
    echo "Data download completed successfully at $(date)"
fi

echo "Summary:"
echo "- Files downloaded to: $raw_reads_dir"
echo "- All logs saved to: $logs_dir"
echo "- Download log: $parallel_log" 
echo "- Error log: $error_log"
echo "- Rename operations: $rename_log"
echo "- Validation results: $validation_log"
