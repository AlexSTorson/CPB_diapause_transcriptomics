#!/bin/bash

#SBATCH --job-name="read_mapping"          # Job name
#SBATCH --output="09.read_mapping_logs/09.read_mapping_%j.out"
#SBATCH --error="09.read_mapping_logs/09.read_mapping_%j.err"
#SBATCH --time=48:00:00                   # Walltime limit
#SBATCH --nodes=2                         # Two nodes for distributed processing
#SBATCH --ntasks-per-node=40              # 40 tasks per node
#SBATCH --partition=short                 # Partition to use
#SBATCH --mem=0                           # Request all available memory on the node

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
scripts_dir="${base_dir}/scripts"
genome_index="$GENOME_INDEX"
trimmed_reads_dir="${base_dir}/trimmed_reads"
alignments_dir="${base_dir}/alignments_sam"

# Standard log setup
current_dir="$(pwd)"
logs_dir="${current_dir}/09.read_mapping_logs"
sample_logs_dir="${logs_dir}/sample_logs"

# Define log files with consistent naming convention
error_log="${logs_dir}/09.read_mapping_errors.log"
progress_log="${logs_dir}/09.read_mapping_progress.log"
parallel_joblog="${logs_dir}/09.read_mapping_parallel_joblog.txt"

# Performance tuning
threads_per_job=4       # Number of threads per HISAT2 job
max_parallel_jobs=15    # Number of parallel jobs across both nodes (80/4=20 but limit to avoid memory issues)

# Create log directories
mkdir -p "$logs_dir" "$sample_logs_dir"

# Clear logs at the start of the script
> "$error_log"
> "$progress_log"

# --- Logging Functions ---
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$progress_log"
}

log_error() {
    echo "[ERROR] [$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$error_log" "$progress_log"
}

# Record start time for performance tracking
start_time=$(date +%s)

# --- Step 1: Load Necessary Modules ---
log_message "Loading required modules"
$MODULE_HISAT2
$MODULE_PARALLEL

# --- Step 2: Create Required Directories ---
log_message "Creating output directory: $alignments_dir"
mkdir -p "$alignments_dir"
if [ $? -ne 0 ]; then
    log_error "Failed to create output directory: $alignments_dir"
    exit 1
fi

# --- Step 3: Validate Genome Index Files ---
log_message "Validating genome index files"
if [ ! -f "${genome_index}.1.ht2" ]; then
    log_error "Genome index files not found: ${genome_index}"
    exit 1
else
    # Count index files to verify completeness
    index_count=$(ls "${genome_index}".*.ht2 2>/dev/null | wc -l)
    log_message "Found $index_count genome index files"
    
    if [ $index_count -lt 8 ]; then
        log_error "Expected 8 index files but found only $index_count"
        exit 1
    fi
fi

# --- Step 4: Get List of Paired-End FASTQ Files ---
log_message "Locating paired-end FASTQ files"

# Find R1 files first
r1_files=()
while IFS= read -r -d '' file; do
    r1_files+=("$file")
done < <(find "$trimmed_reads_dir" -type f -name "*_R1_T.fastq.gz" -print0)

# Check if we found any files
if [ ${#r1_files[@]} -eq 0 ]; then
    log_error "No R1 FASTQ files found in directory: $trimmed_reads_dir"
    exit 1
fi

log_message "Found ${#r1_files[@]} paired-end samples to process"

# --- Step 5: Validate Paired-End Files Exist ---
log_message "Validating paired-end files"
valid_samples=()
missing_r2=0

for r1_file in "${r1_files[@]}"; do
    sample_name=$(basename "$r1_file" | sed 's/_R1_T.fastq.gz//')
    r2_file="${r1_file/_R1_T.fastq.gz/_R2_T.fastq.gz}"
    
    if [ ! -f "$r2_file" ]; then
        log_error "Missing R2 file for sample $sample_name"
        missing_r2=$((missing_r2 + 1))
    else
        valid_samples+=("$sample_name")
    fi
done

# Check if we have valid paired samples
if [ ${#valid_samples[@]} -eq 0 ]; then
    log_error "No valid paired-end samples found"
    exit 1
elif [ $missing_r2 -gt 0 ]; then
    log_error "Found $missing_r2 samples with missing R2 files"
    log_message "Continuing with ${#valid_samples[@]} valid paired-end samples"
else
    log_message "Validated ${#valid_samples[@]} paired-end samples"
fi

# --- Step 6: Define HISAT2 Function for Parallel Execution ---
run_hisat2() {
    local sample_id="$1"
    local r1_file="${trimmed_reads_dir}/${sample_id}_R1_T.fastq.gz"
    local r2_file="${trimmed_reads_dir}/${sample_id}_R2_T.fastq.gz"
    local output_file="${alignments_dir}/${sample_id}.sam"
    local sample_log="${sample_logs_dir}/${sample_id}_hisat2.log"
    local metrics_file="${sample_logs_dir}/${sample_id}_metrics.txt"
    local full_output_file="${sample_logs_dir}/${sample_id}_hisat2_full_output.txt"
    
    # Initialize the log files
    echo "Processing: $sample_id" > "$sample_log"
    echo "R1 file: $r1_file" >> "$sample_log"
    echo "R2 file: $r2_file" >> "$sample_log"
    echo "Output file: $output_file" >> "$sample_log"
    
    # Record start time for this sample
    local sample_start=$(date +%s)
    
    # Prepare the HISAT2 command
    local hisat_cmd="hisat2 -p $threads_per_job --mm --new-summary --dta --time --no-unal --no-mixed --no-discordant -x $genome_index -1 $r1_file -2 $r2_file -S $output_file"
    
    # Run HISAT2 and capture both stderr and stdout
    # stderr contains the alignment statistics, which is what we want to capture in metrics_file
    # We'll use a temporary file to capture the full output
    local temp_output_file=$(mktemp)
    local temp_metrics_file=$(mktemp)
    
    # Record the command in the log
    echo "Alignment Command:" >> "$sample_log"
    echo "  $hisat_cmd" >> "$sample_log"
    
    # Execute the command and capture output
    # Note: we're redirecting stderr to temp_metrics_file (since that's where HISAT2 puts stats)
    # and stdout to temp_output_file
    eval "$hisat_cmd" 2> "$temp_metrics_file" > "$temp_output_file"
    local exit_code=$?
    
    # Copy the metrics to the metrics file (this captures the alignment stats)
    cp "$temp_metrics_file" "$metrics_file"
    
    # Append the metrics to the full output file for completeness
    cat "$temp_metrics_file" > "$full_output_file"
    cat "$temp_output_file" >> "$full_output_file"
    
    # Clean up temp files
    rm -f "$temp_metrics_file" "$temp_output_file"
    
    # Record end time and calculate duration for this sample
    local sample_end=$(date +%s)
    local sample_duration=$((sample_end - sample_start))
    local sample_duration_formatted=$(printf '%dh:%dm:%ds' $((sample_duration/3600)) $((sample_duration%3600/60)) $((sample_duration%60)))
    
    if [ $exit_code -ne 0 ]; then
        echo "HISAT2 failed for sample: $sample_id (Exit code: $exit_code)" >> "$error_log"
        echo "Error details:" >> "$sample_log"
        echo "Full output file contents:" >> "$sample_log"
        cat "$full_output_file" >> "$sample_log"
        return 1
    else
        # Extract alignment statistics with improved parsing
        
        # Try to extract from metrics file
        total_reads=$(grep -m 1 "Total pairs:" "$metrics_file" | awk '{print $3}')
        
        # Different HISAT2 versions might output different formats, try both patterns
        aligned_pairs=$(grep -m 1 "Aligned concordantly 1 time:" "$metrics_file" | awk '{print $5}' | sed 's/[()%]//g')
        if [ -z "$aligned_pairs" ]; then
            aligned_pairs=$(grep -m 1 "Aligned concordantly exactly 1 time:" "$metrics_file" | awk '{print $6}' | sed 's/[()%]//g')
        fi
        
        # Get alignment rate
        alignment_rate=$(grep -m 1 "overall alignment rate:" "$metrics_file" | awk '{print $NF}')
        
        # If parsing failed, try alternative patterns
        if [ -z "$total_reads" ] || [ -z "$aligned_pairs" ] || [ -z "$alignment_rate" ]; then
            # Try parsing from full output
            total_reads=$(grep -m 1 "Total pairs:" "$full_output_file" | awk '{print $3}')
            aligned_pairs=$(grep -m 1 "Aligned concordantly 1 time:" "$full_output_file" | awk '{print $5}' | sed 's/[()%]//g')
            if [ -z "$aligned_pairs" ]; then
                aligned_pairs=$(grep -m 1 "Aligned concordantly exactly 1 time:" "$full_output_file" | awk '{print $6}' | sed 's/[()%]//g')
            fi
            alignment_rate=$(grep -m 1 "overall alignment rate:" "$full_output_file" | awk '{print $NF}')
        fi
        
        # If still missing data, set default values
        [ -z "$total_reads" ] && total_reads="N/A"
        [ -z "$aligned_pairs" ] && aligned_pairs="N/A"
        [ -z "$alignment_rate" ] && alignment_rate="N/A"
        
        # Log success to progress log
        {
            echo "HISAT2 successfully processed sample: $sample_id in $sample_duration_formatted"
            echo "  Total read pairs: $total_reads"
            echo "  Uniquely aligned pairs: $aligned_pairs"
            echo "  Overall alignment rate: $alignment_rate"
        } >> "$progress_log"
        
        # Check SAM file
        if [ -f "$output_file" ]; then
            sam_size=$(stat -c%s "$output_file" 2>/dev/null || echo "N/A")
            if [ "$sam_size" != "N/A" ] && [ "$sam_size" -lt 1000 ]; then  # Arbitrary small size to detect potential issues
                echo "WARNING: SAM file for $sample_id is suspiciously small: $sam_size bytes" >> "$error_log"
                sam_size_formatted="${sam_size}B"
            else
                if [ "$sam_size" != "N/A" ]; then
                    # Convert to human-readable format (GB/MB)
                    sam_size_formatted=$(numfmt --to=iec-i --suffix=B "$sam_size" 2>/dev/null || echo "${sam_size}B")
                else
                    sam_size_formatted="N/A"
                fi
                echo "  SAM file size: $sam_size_formatted" >> "$progress_log"
            fi
        else
            sam_size_formatted="File not found"
            echo "WARNING: SAM file for $sample_id was not created" >> "$error_log"
        fi
        
        # Save detailed stats to sample log
        {
            echo "Alignment Statistics:"
            echo "  Processing time: $sample_duration_formatted"
            echo "  Total read pairs: $total_reads"
            echo "  Uniquely aligned pairs: $aligned_pairs"
            echo "  Overall alignment rate: $alignment_rate"
            echo "  SAM file size: $sam_size_formatted"
            echo ""
            echo "--- Full HISAT2 Output ---"
        } >> "$sample_log"
        
        # Append full metrics to sample log
        cat "$metrics_file" >> "$sample_log"
        
        return 0
    fi
}

# Export functions and variables needed by parallel
export -f run_hisat2 log_error
export trimmed_reads_dir alignments_dir sample_logs_dir error_log progress_log threads_per_job genome_index

# --- Step 7: Run HISAT2 in Parallel ---
log_message "Starting parallel HISAT2 mapping with $max_parallel_jobs jobs"

# Use GNU Parallel to run the HISAT2 jobs
parallel --env run_hisat2 --env log_error --env threads_per_job --env genome_index --env trimmed_reads_dir --env alignments_dir --env error_log --env progress_log --env sample_logs_dir \
    --jobs "$max_parallel_jobs" --joblog "$parallel_joblog" \
    run_hisat2 ::: "${valid_samples[@]}"

# --- Step 8: Validate Output SAM Files ---
log_message "Validating output SAM files"
total_samples=${#valid_samples[@]}
successful_samples=0
failed_samples=0

for sample in "${valid_samples[@]}"; do
    sam_file="${alignments_dir}/${sample}.sam"
    
    if [ -f "$sam_file" ] && [ -s "$sam_file" ]; then
        # Simple validation - check if file starts with SAM header
        if head -n 1 "$sam_file" | grep -q "^@HD"; then
            successful_samples=$((successful_samples + 1))
        else
            log_error "SAM file for $sample appears to be corrupted (no valid header)"
            failed_samples=$((failed_samples + 1))
        fi
    else
        log_error "SAM file for $sample is missing or empty"
        failed_samples=$((failed_samples + 1))
    fi
done

# --- Step 9: Extract Overall Alignment Statistics ---
log_message "Calculating overall alignment statistics"
overall_reads=0
overall_aligned=0
overall_aligned_rate="N/A"

# Create a detailed alignment statistics table
alignment_table="${logs_dir}/09.read_mapping_alignment_table.tsv"
echo -e "Sample\tTotal Reads\tAligned Reads\tAlignment Rate\tSAM Size" > "$alignment_table"

for sample in "${valid_samples[@]}"; do
    metrics_file="${sample_logs_dir}/${sample}_metrics.txt"
    sample_total="N/A"
    sample_aligned="N/A"
    sample_rate="N/A"
    sam_size="N/A"
    
    if [ -f "$metrics_file" ]; then
        # Extract metrics with improved parsing
        sample_total=$(grep -m 1 "Total pairs:" "$metrics_file" | awk '{print $3}')
        
        # Get aligned reads count
        sample_aligned=$(grep -m 1 "Aligned concordantly 1 time:" "$metrics_file" | awk '{print $5}' | sed 's/[()%]//g')
        if [ -z "$sample_aligned" ]; then
            sample_aligned=$(grep -m 1 "Aligned concordantly exactly 1 time:" "$metrics_file" | awk '{print $6}' | sed 's/[()%]//g')
        fi
        
        # Explicitly get alignment rate from overall alignment rate line
        sample_rate=$(grep -m 1 "overall alignment rate" "$metrics_file" | awk '{print $NF}')
        
        # If that doesn't work, calculate it ourselves
        if [ -z "$sample_rate" ] && [ "$sample_total" != "N/A" ] && [ "$sample_aligned" != "N/A" ]; then
            # Remove commas from numbers
            sample_total_clean=$(echo "$sample_total" | tr -d ',')
            sample_aligned_clean=$(echo "$sample_aligned" | tr -d ',')
            
            # Only calculate if we have valid numbers
            if [[ "$sample_total_clean" =~ ^[0-9]+$ ]] && [[ "$sample_aligned_clean" =~ ^[0-9]+$ ]] && [ "$sample_total_clean" -gt 0 ]; then
                sample_rate_val=$(echo "scale=2; ($sample_aligned_clean * 100) / $sample_total_clean" | bc)
                sample_rate="${sample_rate_val}%"
            fi
        fi
    fi
    
    # Get SAM file size
    sam_file="${alignments_dir}/${sample}.sam"
    if [ -f "$sam_file" ]; then
        sam_bytes=$(stat -c%s "$sam_file" 2>/dev/null || echo "N/A")
        if [ "$sam_bytes" != "N/A" ]; then
            sam_size=$(numfmt --to=iec-i --suffix=B "$sam_bytes" 2>/dev/null || echo "${sam_bytes}B")
        fi
    fi
    
    # Make sure we have tab-separated values
    echo -e "${sample}\t${sample_total}\t${sample_aligned}\t${sample_rate}\t${sam_size}" >> "$alignment_table"
    
    # Convert to numeric values for totals (remove commas and non-numeric characters)
    sample_total_num=$(echo "$sample_total" | tr -d ',' | grep -o '[0-9]*' || echo 0)
    sample_aligned_num=$(echo "$sample_aligned" | tr -d ',' | grep -o '[0-9]*' || echo 0)
    
    # Add to totals if numeric and non-zero
    if [ -n "$sample_total_num" ] && [ "$sample_total_num" -gt 0 ]; then
        overall_reads=$((overall_reads + sample_total_num))
    fi
    
    if [ -n "$sample_aligned_num" ] && [ "$sample_aligned_num" -gt 0 ]; then
        overall_aligned=$((overall_aligned + sample_aligned_num))
    fi
done

# Calculate overall alignment rate
if [ "$overall_reads" -gt 0 ] && [ "$overall_aligned" -gt 0 ]; then
    overall_aligned_rate=$(echo "scale=2; ($overall_aligned * 100) / $overall_reads" | bc)"%"
fi

# Add diagnostic info to check file format
log_message "Verifying table format"
head -n 5 "$alignment_table" | cat -A >> "$progress_log"  # Shows tab characters as ^I

# Record end time and calculate duration
end_time=$(date +%s)
duration=$((end_time - start_time))
duration_formatted=$(printf '%dh:%dm:%ds' $((duration/3600)) $((duration%3600/60)) $((duration%60)))

# --- Step 10: Generate Summary Report ---
log_message "Generating summary report"
failed_jobs=$(grep -c 'exit code: [^0]' "$parallel_joblog" 2>/dev/null || echo 0)

summary_file="${logs_dir}/09.read_mapping_summary.txt"
{
    echo "===== HISAT2 Read Mapping Summary ====="
    echo "Date: $(date)"
    echo "Execution time: $duration_formatted"
    echo ""
    echo "Processing Statistics:"
    echo "  Total samples processed: $total_samples"
    echo "  Successfully processed: $successful_samples"
    echo "  Failed samples: $failed_samples"
    echo "  Failed jobs (from parallel): $failed_jobs"
    echo ""
    echo "Alignment Statistics:"
    echo "  Total read pairs: $overall_reads"
    echo "  Uniquely aligned pairs: $overall_aligned"
    echo "  Overall alignment rate: $overall_aligned_rate"
    echo ""
    echo "Detailed alignment statistics per sample are available in: $alignment_table"
    echo ""
    echo "===== End of Summary ====="
} > "$summary_file"

# Also output to progress log

log_message "Read mapping completed. See summary at: $summary_file" false "Read mapping completed. Summary available at: $summary_file"

exit 0
