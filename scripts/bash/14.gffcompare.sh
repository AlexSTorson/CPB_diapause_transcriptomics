#!/bin/bash

#SBATCH --time=04:00:00                 # Walltime for gffcompare
#SBATCH --nodes=1                       # Single node
#SBATCH --ntasks=1                      # Single task
#SBATCH --cpus-per-task=8               # Use 8 threads
#SBATCH --partition=short               # Partition type
#SBATCH --job-name="gffcompare"         # Descriptive job name
#SBATCH --output="14.gffcompare_logs/14.gffcompare_%j.out"
#SBATCH --error="14.gffcompare_logs/14.gffcompare_%j.err"
#SBATCH --mem=32G                       # Request sufficient memory

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
gffcompare_exec="$GFFCOMPARE_EXEC"   # Path to gffcompare executable
reference_gtf="$REFERENCE_GTF"
reference_fna="$REFERENCE_GENOME"
merged_gtf="$MERGED_GTF"              # Use the merged GTF path from config
output_dir="${base_dir}/gffcompare_output"  # Output directory
output_prefix="${output_dir}/gffcompare_result"                # Output prefix

# Standard log setup
current_dir="$(pwd)"
logs_dir="${current_dir}/14.gffcompare_logs"
debug_dir="${logs_dir}/debug_files"

# Define log files with consistent naming convention
error_log="${logs_dir}/14.gffcompare_errors.log"
progress_log="${logs_dir}/14.gffcompare_progress.log"
gffcompare_log="${logs_dir}/14.gffcompare_output.log"
debug_log="${logs_dir}/14.gffcompare_debug.log"

# Create required directories
mkdir -p "$logs_dir" "$debug_dir"

# Clear logs at the start of the script
> "$error_log"
> "$progress_log"
> "$gffcompare_log"
> "$debug_log"

# --- Logging Functions ---
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$progress_log"
}

log_error() {
    echo "[ERROR] [$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$error_log"
}

log_debug() {
    echo "[DEBUG] [$(date '+%Y-%m-%d %H:%M:%S')] $1" >> "$debug_log"
}

# Function to extract more accurate counts from gffcompare outputs
extract_accurate_counts() {
    local stats_file="$1"
    local gtf_file="$2"
    local name="$3"
    
    # Log the stats file content for debugging
    log_debug "$name stats file content:"
    log_debug "$(cat "$stats_file")"
    
    # Try to extract transcript and gene counts using various patterns
    # Extract transcript count
    local transcript_count=""
    
    # Pattern 1: Query mRNAs
    if [ -z "$transcript_count" ]; then
        transcript_count=$(grep "Query mRNAs:" "$stats_file" | awk '{print $4}')
        if [ -n "$transcript_count" ]; then
            log_debug "Found transcript count using 'Query mRNAs:' pattern: $transcript_count"
        fi
    fi
    
    # Pattern 2: Matching intron chains
    if [ -z "$transcript_count" ]; then
        transcript_count=$(grep "Matching intron chains:" "$stats_file" | awk '{print $4}')
        if [ -n "$transcript_count" ]; then
            log_debug "Found transcript count using 'Matching intron chains:' pattern: $transcript_count"
        fi
    fi
    
    # Pattern 3: Total query mRNAs (some versions)
    if [ -z "$transcript_count" ]; then
        transcript_count=$(grep "Total query mRNAs:" "$stats_file" | awk '{print $4}')
        if [ -n "$transcript_count" ]; then
            log_debug "Found transcript count using 'Total query mRNAs:' pattern: $transcript_count"
        fi
    fi
    
    # Extract gene count
    local gene_count=""
    
    # Pattern 1: Query loci
    if [ -z "$gene_count" ]; then
        gene_count=$(grep "Query loci:" "$stats_file" | awk '{print $3}')
        if [ -n "$gene_count" ]; then
            log_debug "Found gene count using 'Query loci:' pattern: $gene_count"
        fi
    fi
    
    # Pattern 2: Matching loci
    if [ -z "$gene_count" ]; then
        gene_count=$(grep "Matching loci:" "$stats_file" | awk '{print $3}')
        if [ -n "$gene_count" ]; then
            log_debug "Found gene count using 'Matching loci:' pattern: $gene_count"
        fi
    fi
    
    # Pattern 3: Total loci (some versions)
    if [ -z "$gene_count" ]; then
        gene_count=$(grep "Total loci:" "$stats_file" | awk '{print $3}')
        if [ -n "$gene_count" ]; then
            log_debug "Found gene count using 'Total loci:' pattern: $gene_count"
        fi
    fi
    
    # Fallback to alternative methods if needed
    if [ -z "$gene_count" ]; then
        log_debug "Using alternative method to count genes for $name"
        
        # Try to extract from the tmap file
        local tmap_file="${output_prefix}.tmap"
        if [ -f "$tmap_file" ]; then
            # Count unique ref_gene_id entries
            gene_count=$(grep -v "^#" "$tmap_file" | awk '{print $2}' | sort -u | wc -l)
            log_debug "Gene count from tmap method: $gene_count"
        else
            log_debug "No tmap file found at $tmap_file"
        fi
    fi
    
    # If still empty, try to count loci from the GTF file directly
    if [ -z "$gene_count" ]; then
        log_debug "Counting genes directly from GTF file for $name"
        
        # Count unique gene_id entries
        gene_count=$(grep -v "^#" "$gtf_file" | grep -o 'gene_id "[^"]*"' | sort -u | wc -l)
        log_debug "Gene count from GTF gene_id attributes: $gene_count"
    fi
    
    # Set defaults if still empty
    if [ -z "$transcript_count" ]; then
        log_debug "Setting default transcript count for $name"
        # Count transcripts directly from GTF
        transcript_count=$(grep -v "^#" "$gtf_file" | grep -o 'transcript_id "[^"]*"' | sort -u | wc -l)
        if [ -z "$transcript_count" ]; then
            transcript_count=0
        fi
    fi
    
    if [ -z "$gene_count" ]; then
        log_debug "Setting default gene count for $name"
        gene_count=0
    fi
    
    # Return the counts
    echo "${transcript_count}:${gene_count}"
}

# Record start time for performance tracking
start_time=$(date +%s)

# --- Step 1: Validate Input Files ---
log_message "Validating input files"
log_message "Using merge method: $MERGE_METHOD (Merged GTF: $merged_gtf)"

# Check gffcompare executable
if [ ! -f "$gffcompare_exec" ]; then
    log_error "gffcompare executable not found: $gffcompare_exec"
    exit 1
fi

if [ ! -x "$gffcompare_exec" ]; then
    log_error "gffcompare executable is not executable: $gffcompare_exec"
    exit 1
fi

# Check reference GTF
if [ ! -f "$reference_gtf" ]; then
    log_error "Reference GTF file not found: $reference_gtf"
    exit 1
fi

if [ ! -s "$reference_gtf" ]; then
    log_error "Reference GTF file is empty: $reference_gtf"
    exit 1
fi

# Get reference GTF stats for logging
ref_gtf_size=$(stat -c%s "$reference_gtf")
ref_gtf_size_mb=$(echo "scale=2; $ref_gtf_size / 1048576" | bc)

# Check reference genome
if [ ! -f "$reference_fna" ]; then
    log_error "Reference genome file not found: $reference_fna"
    exit 1
fi

if [ ! -s "$reference_fna" ]; then
    log_error "Reference genome file is empty: $reference_fna"
    exit 1
fi

ref_fna_size=$(stat -c%s "$reference_fna")
ref_fna_size_mb=$(echo "scale=2; $ref_fna_size / 1048576" | bc)

log_message "Reference genome: $reference_fna"
log_message "  Size: ${ref_fna_size_mb} MB"

# Check merged GTF
if [ ! -f "$merged_gtf" ]; then
    log_error "Merged GTF file not found: $merged_gtf"
    exit 1
fi

if [ ! -s "$merged_gtf" ]; then
    log_error "Merged GTF file is empty: $merged_gtf"
    exit 1
fi

# Get merged GTF stats for logging
merged_gtf_size=$(stat -c%s "$merged_gtf")
merged_gtf_size_mb=$(echo "scale=2; $merged_gtf_size / 1048576" | bc)

log_message "All input files validated successfully"

# --- Step 2: Prepare Output Directory ---
log_message "Creating output directory: $output_dir"
mkdir -p "$output_dir"
if [ $? -ne 0 ]; then
    log_error "Failed to create output directory: $output_dir"
    exit 1
fi

# --- Step 3: Run gffcompare ---
log_message "Starting gffcompare with command:"
log_message "$gffcompare_exec -R -r $reference_gtf -s $reference_fna -o $output_prefix $merged_gtf"

# Run gffcompare and capture output
"$gffcompare_exec" -R -r "$reference_gtf" -s "$reference_fna" -o "$output_prefix" "$merged_gtf" &> "$gffcompare_log"

# Check if gffcompare ran successfully
if [ $? -ne 0 ]; then
    log_error "gffcompare failed. Check $gffcompare_log for details."
    exit 1
else
    log_message "gffcompare command completed successfully"
fi

# --- Step 4: Validate Output Files ---
log_message "Validating gffcompare output files"

# Define expected output files
expected_files=(
    "${output_prefix}.tracking"
    "${output_prefix}.stats"
    "${output_prefix}.loci"
    "${output_prefix}.annotated.gtf"
    "${output_prefix}.tmap"
    "${output_prefix}.refmap"
)

missing_outputs=0
output_sizes=()

for file in "${expected_files[@]}"; do
    if [ ! -f "$file" ]; then
        log_error "Expected output file not found: $file"
        missing_outputs=$((missing_outputs + 1))
    elif [ ! -s "$file" ]; then
        log_error "Output file is empty: $file"
        missing_outputs=$((missing_outputs + 1))
    else
        # Get file size for reporting
        file_size=$(stat -c%s "$file")
        file_size_kb=$(echo "scale=2; $file_size / 1024" | bc)
        output_sizes+=("$(basename "$file"): ${file_size_kb} KB")
        log_message "Output file: $(basename "$file") (${file_size_kb} KB)"
    fi
done

if [ $missing_outputs -gt 0 ]; then
    log_error "Found $missing_outputs missing or empty output files"
else
    log_message "All expected output files are present and non-empty"
fi

# --- Step 5: Extract Key Statistics using the same method as other scripts ---
log_message "Extracting key statistics from output files"

# Run the extraction function for the reference GTF
stats_file="${output_prefix}.stats"
if [ -f "$stats_file" ]; then
    # Get reference counts 
    log_message "Getting accurate reference statistics"
    ref_counts=$(extract_accurate_counts "$stats_file" "$reference_gtf" "Reference")
    ref_transcript_count=$(echo "$ref_counts" | cut -d':' -f1)
    ref_gene_count=$(echo "$ref_counts" | cut -d':' -f2)
    
    # Get merged GTF counts
    log_message "Getting accurate merged assembly statistics"
    merged_counts=$(extract_accurate_counts "$stats_file" "$merged_gtf" "Merged")
    merged_transcript_count=$(echo "$merged_counts" | cut -d':' -f1)
    merged_gene_count=$(echo "$merged_counts" | cut -d':' -f2)

    # Count matching statistics
    matching_transcripts=$(grep "Matching intron chains" "$stats_file" | awk '{print $4}')
    matching_loci=$(grep "Matching loci" "$stats_file" | awk '{print $4}')
    
    log_message "Key statistics from gffcompare:"
    log_message "  Reference transcripts: $ref_transcript_count"
    log_message "  Reference genes: $ref_gene_count"
    log_message "  Merged transcripts: $merged_transcript_count"
    log_message "  Merged genes: $merged_gene_count"
    log_message "  Matching transcripts: $matching_transcripts"
    log_message "  Matching loci: $matching_loci"
else
    log_error "Stats file not found: $stats_file"
    ref_transcript_count="N/A"
    ref_gene_count="N/A"
    merged_transcript_count="N/A"
    merged_gene_count="N/A"
    matching_transcripts="N/A"
    matching_loci="N/A"
fi

# Record end time and calculate duration
end_time=$(date +%s)
duration=$((end_time - start_time))
duration_formatted=$(printf '%dm:%ds' $((duration/60)) $((duration%60)))

# --- Step 6: Generate Summary Report ---
summary_file="${logs_dir}/14.gffcompare_summary.txt"
{
    echo "===== gffcompare Summary ====="
    echo "Date: $(date)"
    echo "Execution time: $duration_formatted"
    echo
    echo "Pipeline Configuration:"
    echo "  Merge method: $MERGE_METHOD"
    echo "  Merged GTF file: $merged_gtf"
    echo
    echo "Input Files:"
    echo "  Reference GTF: $reference_gtf"
    echo "  Reference GTF size: ${ref_gtf_size_mb} MB"
    echo "  Reference transcripts: $ref_transcript_count"
    echo "  Reference genes: $ref_gene_count"
    echo
    echo "  Reference genome: $reference_fna"
    echo "  Reference genome size: ${ref_fna_size_mb} MB"
    echo
    echo "  Merged GTF: $merged_gtf"
    echo "  Merged GTF size: ${merged_gtf_size_mb} MB"
    echo "  Merged transcripts: $merged_transcript_count"
    echo "  Merged genes: $merged_gene_count"
    echo
    echo "Comparison Results:"
    echo "  Matching transcripts: $matching_transcripts"
    echo "  Matching loci: $matching_loci"  
    echo
    echo "Output:"
    echo "  Output directory: $output_dir"
    echo "  Output prefix: $output_prefix"
    echo
    echo "Output Files:"
    for size in "${output_sizes[@]}"; do
        echo "  $size"
    done
    echo
    if [ $missing_outputs -gt 0 ]; then
        echo "WARNING: $missing_outputs expected output files are missing or empty"
    fi
    echo "======================================="
} > "$summary_file"

# --- Step 7: Extract Detailed Statistics to Separate File ---
if [ -f "${output_prefix}.stats" ]; then
    detailed_stats_file="${logs_dir}/14.gffcompare_detailed_stats.txt"
    {
        echo "===== Detailed gffcompare Statistics ====="
        echo "Date: $(date)"
        echo
        echo "Pipeline Configuration:"
        echo "  Merge method: $MERGE_METHOD"
        echo "  Merged GTF file: $merged_gtf"
        echo
        echo "Statistics from ${output_prefix}.stats:"
        echo "----------------------------------------"
        cat "${output_prefix}.stats"
        echo
        echo "======================================="
    } > "$detailed_stats_file"
    
    log_message "Detailed statistics saved to: $detailed_stats_file"
fi

# --- Step 8: Log Completion ---
if [ $missing_outputs -gt 0 ]; then
    log_error "gffcompare completed with warnings. Some expected output files are missing or empty."
    cat "$summary_file"
    exit 0  # Don't fail as at least some output was produced
elif [ -s "$error_log" ]; then
    log_error "gffcompare completed with errors. Check $error_log for details."
    cat "$summary_file"
    exit 1
else
    log_message "gffcompare completed successfully in $duration_formatted"
    cat "$summary_file"
    exit 0
fi
