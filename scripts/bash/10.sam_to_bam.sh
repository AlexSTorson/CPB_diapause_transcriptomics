#!/bin/bash

#SBATCH --job-name="sam_to_bam"
#SBATCH --output="10.sam_to_bam_logs/10.sam_to_bam_%j.out"
#SBATCH --error="10.sam_to_bam_logs/10.sam_to_bam_%j.err"
#SBATCH --time=24:00:00                        # Reduced walltime based on expected performance
#SBATCH --nodes=2                               # Two nodes
#SBATCH --ntasks-per-node=40                   # Number of tasks per node
#SBATCH --partition=short                      
#SBATCH --mem=0                                # Request all available memory

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



#!/bin/bash

#!/bin/bash
# SLURM job submission settings

# --- Centralized Variables ---
base_dir="$PROJECT_DIR"
sam_dir="${base_dir}/alignments_sam"
bam_dir="${base_dir}/alignments_bam"
indexed_bam_dir="${base_dir}/alignments_bam_indexed"  # New directory for indexed BAM files

# Standard log setup
current_dir="$(pwd)"
logs_dir="${current_dir}/10.sam_to_bam_logs"

# Define log files with consistent naming convention
error_log="${logs_dir}/10.sam_to_bam_errors.log"
progress_log="${logs_dir}/10.sam_to_bam_progress.log"
parallel_joblog="${logs_dir}/10.sam_to_bam_parallel_joblog.txt"

# Performance tuning
threads_per_job=4        # Threads per samtools sort job
max_parallel_jobs=15     # Adjusted for memory constraints across two nodes
memory_per_job="8G"      # Memory per sorting job

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
$MODULE_SAMTOOLS
$MODULE_PARALLEL

# --- Step 2: Create Output Directories ---
log_message "Creating output directories"
for dir in "$bam_dir" "$indexed_bam_dir"; do
    mkdir -p "$dir"
    if [ $? -ne 0 ]; then
        log_error "Failed to create directory: $dir"
        exit 1
    fi
done

# --- Step 3: Validate Input Files ---
log_message "Validating input SAM files"
sam_files=()
while IFS= read -r -d '' file; do
    sam_files+=("$file")
done < <(find "$sam_dir" -type f -name "*.sam" -print0)

# Check if we found any files
if [ ${#sam_files[@]} -eq 0 ]; then
    log_error "No SAM files found in directory: $sam_dir"
    exit 1
fi

log_message "Found ${#sam_files[@]} SAM files to process"

# --- Step 4: Define SAM to BAM Conversion Function ---
sam_to_bam() {
    local sam_file="$1"
    local sample_name=$(basename "$sam_file" .sam)
    local temp_bam="${bam_dir}/${sample_name}.temp.bam"
    local sorted_bam="${bam_dir}/${sample_name}.bam"
    local indexed_bam="${indexed_bam_dir}/${sample_name}.bam"
    local sample_log="${logs_dir}/sample_logs/${sample_name}.log"
    
    echo "Processing: $sample_name" > "$sample_log"
    
    # Record start time for this sample
    local sample_start=$(date +%s)
    
    # Step 1: Convert SAM to BAM (without sorting)
    echo "Converting SAM to BAM..." >> "$sample_log"
    samtools view -@ 2 -S -b "$sam_file" -o "$temp_bam" &>> "$sample_log"
    
    if [ $? -ne 0 ]; then
        echo "SAM to BAM conversion failed for $sample_name" >> "$error_log"
        return 1
    fi
    
    # Step 2: Sort BAM
    echo "Sorting BAM file..." >> "$sample_log"
    samtools sort -@ "$threads_per_job" -m "$memory_per_job" -o "$sorted_bam" "$temp_bam" &>> "$sample_log"
    
    if [ $? -ne 0 ]; then
        echo "BAM sorting failed for $sample_name" >> "$error_log"
        return 1
    fi
    
    # Step 3: Remove temporary BAM file
    rm -f "$temp_bam"
    
    # Step 4: Create index for sorted BAM
    echo "Indexing sorted BAM file..." >> "$sample_log"
    samtools index -@ "$threads_per_job" "$sorted_bam" &>> "$sample_log"
    
    if [ $? -ne 0 ]; then
        echo "BAM indexing failed for $sample_name" >> "$error_log"
        return 1
    fi
    
    # Step 5: Copy indexed BAM to the indexed directory (optional)
    cp "$sorted_bam" "$indexed_bam"
    cp "$sorted_bam.bai" "$indexed_bam.bai"
    
    # Get file sizes for logging
    sam_size=$(stat -c%s "$sam_file")
    bam_size=$(stat -c%s "$sorted_bam")
    bai_size=$(stat -c%s "$sorted_bam.bai")
    
    # Calculate compression ratio
    compression_ratio=$(echo "scale=2; $sam_size / $bam_size" | bc)
    
    # Record end time and calculate duration for this sample
    local sample_end=$(date +%s)
    local sample_duration=$((sample_end - sample_start))
    local sample_duration_formatted=$(printf '%dm:%ds' $((sample_duration/60)) $((sample_duration%60)))
    
    echo "Completed processing $sample_name in $sample_duration_formatted" >> "$progress_log"
    echo "  SAM size: $(numfmt --to=iec-i --suffix=B $sam_size)" >> "$progress_log"
    echo "  BAM size: $(numfmt --to=iec-i --suffix=B $bam_size)" >> "$progress_log"
    echo "  Index size: $(numfmt --to=iec-i --suffix=B $bai_size)" >> "$progress_log"
    echo "  Compression ratio: ${compression_ratio}:1" >> "$progress_log"
    
    # Add detailed stats to sample log
    echo "Conversion Statistics:" >> "$sample_log"
    echo "  Processing time: $sample_duration_formatted" >> "$sample_log"
    echo "  SAM size: $(numfmt --to=iec-i --suffix=B $sam_size)" >> "$sample_log"
    echo "  BAM size: $(numfmt --to=iec-i --suffix=B $bam_size)" >> "$sample_log"
    echo "  Index size: $(numfmt --to=iec-i --suffix=B $bai_size)" >> "$sample_log"
    echo "  Compression ratio: ${compression_ratio}:1" >> "$sample_log"
    
    return 0
}

export -f sam_to_bam
export bam_dir indexed_bam_dir threads_per_job memory_per_job error_log progress_log logs_dir

# --- Step 5: Run Conversion in Parallel ---
log_message "Starting parallel SAM to BAM conversion with $max_parallel_jobs jobs"

parallel -u --env bam_dir --env indexed_bam_dir --env threads_per_job --env memory_per_job --env error_log --env progress_log --env logs_dir \
    --jobs "$max_parallel_jobs" --no-notice --joblog "$parallel_joblog" sam_to_bam ::: "${sam_files[@]}"

# --- Step 6: Validate Converted Files ---
log_message "Validating BAM and index files"
total_files=${#sam_files[@]}
successful_bams=0
successful_indexes=0
failed_conversions=0

for sam_file in "${sam_files[@]}"; do
    sample_name=$(basename "$sam_file" .sam)
    bam_file="${bam_dir}/${sample_name}.bam"
    bai_file="${bam_dir}/${sample_name}.bam.bai"
    
    if [ -f "$bam_file" ] && [ -s "$bam_file" ]; then
        # Check if BAM file can be read with samtools
        if samtools view -H "$bam_file" &>/dev/null; then
            successful_bams=$((successful_bams + 1))
        else
            log_error "BAM file for $sample_name exists but appears corrupted"
            failed_conversions=$((failed_conversions + 1))
        fi
    else
        log_error "BAM file for $sample_name is missing or empty"
        failed_conversions=$((failed_conversions + 1))
    fi
    
    if [ -f "$bai_file" ] && [ -s "$bai_file" ]; then
        successful_indexes=$((successful_indexes + 1))
    else
        log_error "Index file for $sample_name is missing or empty"
    fi
done

# --- Step 7: Calculate Size and Compression Statistics ---
log_message "Calculating file size statistics"
total_sam_size=0
total_bam_size=0

for sam_file in "${sam_files[@]}"; do
    sample_name=$(basename "$sam_file" .sam)
    bam_file="${bam_dir}/${sample_name}.bam"
    
    if [ -f "$sam_file" ] && [ -f "$bam_file" ]; then
        # Get file sizes
        sam_size=$(stat -c%s "$sam_file")
        bam_size=$(stat -c%s "$bam_file")
        
        # Add to totals
        total_sam_size=$((total_sam_size + sam_size))
        total_bam_size=$((total_bam_size + bam_size))
    fi
done

# Calculate overall compression ratio
if [ $total_bam_size -gt 0 ]; then
    overall_ratio=$(echo "scale=2; $total_sam_size / $total_bam_size" | bc)
else
    overall_ratio="N/A"
fi

# Record end time and calculate duration
end_time=$(date +%s)
duration=$((end_time - start_time))
duration_formatted=$(printf '%dh:%dm:%ds' $((duration/3600)) $((duration%3600/60)) $((duration%60)))

# --- Step 8: Generate Summary Report ---
failed_jobs=$(grep -c 'exit code: [^0]' "$parallel_joblog" 2>/dev/null || echo 0)

summary_file="${logs_dir}/10.sam_to_bam_summary.txt"
{
    echo "===== SAM to BAM Conversion Summary ====="
    echo "Date: $(date)"
    echo "Execution time: $duration_formatted"
    echo
    echo "Processing Statistics:"
    echo "  Total files processed: $total_files"
    echo "  Successfully converted BAMs: $successful_bams"
    echo "  Successfully indexed BAMs: $successful_indexes"
    echo "  Failed conversions: $failed_conversions"
    echo "  Failed jobs (from parallel): $failed_jobs"
    echo
    echo "Storage Statistics:"
    echo "  Total SAM size: $(numfmt --to=iec-i --suffix=B $total_sam_size)"
    echo "  Total BAM size: $(numfmt --to=iec-i --suffix=B $total_bam_size)"
    echo "  Overall compression ratio: ${overall_ratio}:1"
    echo
    if [ $failed_conversions -gt 0 ]; then
        echo "Failed conversions:"
        grep -E 'conversion failed|sorting failed|indexing failed|exists but appears corrupted|is missing or empty' "$error_log" | \
        sed -E 's/.*failed for ([^[:space:]]+).*/  - \1/' | \
        sed -E 's/.*file for ([^[:space:]]+) is.*/  - \1/' | \
        sed -E 's/.*file for ([^[:space:]]+) exists.*/  - \1/' | \
        sort | uniq
    fi
    echo
    echo "Output:"
    echo "  BAM files location: $bam_dir"
    echo "  Indexed BAM files: $indexed_bam_dir"
    echo "  Sample logs location: ${logs_dir}/sample_logs"
    echo "==============================================="
} > "$summary_file"

# --- Step 9: Log Completion ---
if [ $failed_conversions -gt 0 ]; then
    log_error "SAM to BAM conversion completed with $failed_conversions errors. See $error_log for details."
    cat "$summary_file"
    exit 1
else
    log_message "SAM to BAM conversion completed successfully for all $total_files files in $duration_formatted"
    cat "$summary_file"
    exit 0
fi
