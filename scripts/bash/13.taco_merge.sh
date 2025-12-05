#!/bin/bash

#SBATCH --time=04:00:00                 # Walltime limit (HH:MM:SS)
#SBATCH --nodes=1                       # Single node
#SBATCH --ntasks=1                      # Single task (merging is a single process)
#SBATCH --cpus-per-task=20              # Increased threads for TACO
#SBATCH --partition=short               # Partition type
#SBATCH --job-name="taco_merge"         # Descriptive job name
#SBATCH --output="13.taco_merge_logs/13.taco_merge_%j.out"
#SBATCH --error="13.taco_merge_logs/13.taco_merge_%j.err"
#SBATCH --mem=120G                      # Request sufficient memory for merge operation

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
assemblies_dir="${base_dir}/assemblies"          # Assemblies directory
mergelist_file="${assemblies_dir}/mergelist.txt" # Mergelist file
reference_gtf="$REFERENCE_GTF"
gffcompare_exec="$GFFCOMPARE_EXEC"  # Path to gffcompare executable from config

# Use a fixed output directory
merged_output="${base_dir}/assemblies/taco_assembly_merge"  # Fixed TACO output directory
taco_executable="$TACO_EXEC"      # Path to TACO executable

# Standard log setup
current_dir="$(pwd)"
logs_dir="${current_dir}/13.taco_merge_logs"
gffcompare_dir="${logs_dir}/gffcompare_output"

# Define log files with consistent naming convention
error_log="${logs_dir}/13.taco_merge_errors.log"
progress_log="${logs_dir}/13.taco_merge_progress.log"
taco_log="${logs_dir}/13.taco_merge_output.log"
gffcompare_log="${logs_dir}/13.taco_merge_gffcompare.log"
debug_log="${logs_dir}/13.taco_merge_debug.log"

# Performance tuning
threads=20  # Number of threads for TACO
max_taco_runtime=10800  # Maximum TACO runtime in seconds (3 hours)

# Create necessary directories
mkdir -p "$logs_dir" "$gffcompare_dir"

# Clear logs at the start of the script
> "$error_log"
> "$progress_log"
> "$taco_log"
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

# Function to run GFFCompare for accurate gene and transcript counting
analyze_with_gffcompare() {
    local gtf_file="$1"
    local output_prefix="$2"
    
    # Log the command we're about to run
    log_debug "Running GFFCompare with command: $gffcompare_exec -r $reference_gtf -o $output_prefix $gtf_file"
    
    # Run GFFCompare to analyze the merged GTF
    "$gffcompare_exec" -r "$reference_gtf" -o "$output_prefix" "$gtf_file" &> "$gffcompare_log"
    
    if [ $? -ne 0 ]; then
        log_error "GFFCompare analysis failed. Check $gffcompare_log for details."
        return 1
    fi
    
    # Extract stats from the GFFCompare output
    local stats_file="${output_prefix}.stats"
    
    if [ ! -f "$stats_file" ]; then
        log_error "GFFCompare stats file not found: $stats_file"
        return 1
    fi
    
    # Copy stats file to log directory for reference
    cp "$stats_file" "${logs_dir}/gffcompare_stats.txt"
    
    # Log the stats file content for debugging
    log_debug "GFFCompare stats file content:"
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
        log_debug "Using alternative method to count genes"
        
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
        log_debug "Counting genes directly from GTF file"
        
        # Count unique gene_id entries
        gene_count=$(grep -v "^#" "$gtf_file" | grep -o 'gene_id "[^"]*"' | sort -u | wc -l)
        log_debug "Gene count from GTF gene_id attributes: $gene_count"
    fi
    
    # Set defaults if still empty
    if [ -z "$transcript_count" ]; then
        log_debug "Setting default transcript count"
        # Count transcripts directly from GTF
        transcript_count=$(grep -v "^#" "$gtf_file" | grep -o 'transcript_id "[^"]*"' | sort -u | wc -l)
        if [ -z "$transcript_count" ]; then
            transcript_count=0
        fi
    fi
    
    if [ -z "$gene_count" ]; then
        log_debug "Setting default gene count"
        gene_count=0
    fi
    
    # Return the counts
    echo "${transcript_count}:${gene_count}"
}

# Record start time for performance tracking
start_time=$(date +%s)

# --- Step 1: Validate TACO Executable ---
log_message "Validating TACO executable"

if [ ! -f "$taco_executable" ]; then
    log_error "TACO executable not found: $taco_executable"
    exit 1
fi

if [ ! -x "$taco_executable" ]; then
    log_error "TACO executable is not executable: $taco_executable"
    exit 1
fi

# --- Step 2: Validate Input Files ---
log_message "Validating input files"

# Check reference GTF
if [ ! -f "$reference_gtf" ]; then
    log_error "Reference GTF file not found: $reference_gtf"
    exit 1
fi

# Check GFFCompare executable
if [ ! -f "$gffcompare_exec" ]; then
    log_error "GFFCompare executable not found: $gffcompare_exec"
    exit 1
fi

if [ ! -x "$gffcompare_exec" ]; then
    log_error "GFFCompare executable is not executable: $gffcompare_exec"
    exit 1
fi

# Get reference GTF stats for logging
ref_size=$(stat -c%s "$reference_gtf")
ref_size_mb=$(echo "scale=2; $ref_size / 1048576" | bc)

# Check mergelist file
if [ ! -f "$mergelist_file" ]; then
    log_error "Mergelist file not found: $mergelist_file"
    exit 1
fi

# Check content of mergelist
mergelist_count=$(wc -l < "$mergelist_file")
if [ $mergelist_count -eq 0 ]; then
    log_error "Mergelist file is empty: $mergelist_file"
    exit 1
fi

log_message "Mergelist contains $mergelist_count GTF files"

# Create a local copy of the mergelist with absolute paths
local_mergelist="${logs_dir}/local_mergelist.txt"
> "$local_mergelist"
while IFS= read -r gtf_file; do
    # Convert to absolute path if not already
    if [[ "$gtf_file" != /* ]]; then
        gtf_file="$(readlink -f "$gtf_file")"
    fi
    echo "$gtf_file" >> "$local_mergelist"
done < "$mergelist_file"

# Sample check: Log the first few GTF files in the mergelist
log_message "First 5 entries in mergelist (or all if fewer):"
head -n 5 "$local_mergelist" | tee -a "$progress_log"

# Verify all files in mergelist exist
missing_files=0
total_gtf_size=0
list_of_missing_files=""

while IFS= read -r gtf_file; do
    if [ ! -f "$gtf_file" ]; then
        log_error "GTF file listed in mergelist does not exist: $gtf_file"
        missing_files=$((missing_files + 1))
        list_of_missing_files="${list_of_missing_files}${gtf_file}\n"
    else
        # Add to total size
        file_size=$(stat -c%s "$gtf_file")
        total_gtf_size=$((total_gtf_size + file_size))
    fi
done < "$local_mergelist"

if [ $missing_files -gt 0 ]; then
    log_error "Found $missing_files missing GTF files in mergelist"
    exit 1
fi

total_gtf_size_mb=$(echo "scale=2; $total_gtf_size / 1048576" | bc)
log_message "Total size of all GTF files to merge: ${total_gtf_size_mb} MB"

# Log the output directory we'll be using
log_message "Using fixed output directory: $merged_output"

# --- Step 3: Remove any existing output directory ---
log_message "Removing existing output directory if it exists: $merged_output"
if [ -d "$merged_output" ] || [ -e "$merged_output" ]; then
    log_message "Found existing directory, removing it"
    rm -rf "$merged_output"
    if [ $? -ne 0 ]; then
        log_error "Failed to remove existing output directory: $merged_output"
        exit 1
    fi
    # Verify removal was successful
    if [ -d "$merged_output" ] || [ -e "$merged_output" ]; then
        log_error "Directory removal failed: $merged_output still exists"
        exit 1
    fi
    log_message "Successfully removed existing directory"
else
    log_message "No existing directory found"
fi

# --- Step 4: Run TACO ---
log_message "Starting TACO merge process with $threads threads"
taco_cmd="$taco_executable -p $threads -o $merged_output $local_mergelist"
log_message "TACO command: $taco_cmd"

# Run TACO with monitoring
timeout $max_taco_runtime $taco_cmd &> "$taco_log" &
taco_pid=$!

# Log the PID for reference
log_message "TACO process started with PID: $taco_pid"

# Monitor progress
monitor_interval=300  # Check every 5 minutes
last_size=0

while kill -0 $taco_pid 2>/dev/null; do
    # Process is still running
    elapsed=$(($(date +%s) - start_time))
    elapsed_formatted=$(printf '%dh:%dm:%ds' $((elapsed/3600)) $((elapsed%3600/60)) $((elapsed%60)))
    
    # Check if output directory size is growing
    if [ -d "$merged_output" ]; then
        dir_size_bytes=$(du -sb "$merged_output" 2>/dev/null | awk '{print $1}')
        dir_size=$(du -sh "$merged_output" 2>/dev/null | awk '{print $1}')
        
        if [ -n "$dir_size_bytes" ] && [ -n "$last_size" ]; then
            size_diff=$((dir_size_bytes - last_size))
            size_diff_mb=$(echo "scale=2; $size_diff / 1048576" | bc)
            log_message "TACO still running after $elapsed_formatted. Current output size: $dir_size (changed by ${size_diff_mb} MB)"
            last_size=$dir_size_bytes
        else
            log_message "TACO still running after $elapsed_formatted. Current output size: $dir_size"
            last_size=${dir_size_bytes:-0}
        fi
    else
        log_message "TACO still running after $elapsed_formatted. Output directory not found yet."
    fi
    
    # Check if there's recent activity in the log file
    if [ -f "$taco_log" ]; then
        recent_log=$(tail -n 5 "$taco_log" 2>/dev/null)
        log_message "Recent log output: ${recent_log:-No recent log output}"
    else
        log_message "No log file found yet."
    fi
    
    # Check if maximum runtime exceeded
    if [ $elapsed -gt $max_taco_runtime ]; then
        log_error "TACO exceeded maximum runtime of $max_taco_runtime seconds. Terminating process."
        kill -9 $taco_pid 2>/dev/null
        break
    fi
    
    sleep $monitor_interval
done

# Check if TACO completed successfully
wait $taco_pid
taco_exit_code=$?

if [ $taco_exit_code -ne 0 ]; then
    log_error "TACO merge failed with exit code $taco_exit_code. Check $taco_log for details."
    log_error "Last 20 lines of TACO log:"
    tail -n 20 "$taco_log" | tee -a "$error_log"
    exit 1
else
    log_message "TACO merge command completed successfully"
fi

# --- Step 5: Validate Output Files ---
log_message "Validating TACO output files"

assembly_gtf="${merged_output}/assembly.gtf"
if [ ! -f "$assembly_gtf" ]; then
    log_error "TACO output assembly.gtf file not found: $assembly_gtf"
    exit 1
fi

if [ ! -s "$assembly_gtf" ]; then
    log_error "TACO output assembly.gtf file is empty: $assembly_gtf"
    exit 1
fi

# Only check for the essential file
expected_files=(
    "${merged_output}/assembly.gtf"
)

missing_outputs=0
for file in "${expected_files[@]}"; do
    if [ ! -f "$file" ]; then
        log_error "Expected TACO output file not found: $file"
        missing_outputs=$((missing_outputs + 1))
    fi
done

if [ $missing_outputs -gt 0 ]; then
    log_error "Assembly GTF file missing - TACO merge failed"
    exit 1
else
    log_message "Assembly GTF file found and validated"
fi

# --- Step 6: Run GFFCompare on reference GTF ---
log_message "Analyzing reference GTF with GFFCompare"
ref_gffcompare_prefix="${gffcompare_dir}/reference"
ref_counts=$(analyze_with_gffcompare "$reference_gtf" "$ref_gffcompare_prefix")
ref_transcript_count=$(echo "$ref_counts" | cut -d':' -f1)
ref_gene_count=$(echo "$ref_counts" | cut -d':' -f2)

log_message "Reference GTF file: $reference_gtf"
log_message "  Size: ${ref_size_mb} MB"
log_message "  Transcripts (from GFFCompare): $ref_transcript_count"
log_message "  Genes (from GFFCompare): $ref_gene_count"

# --- Step 7: Analyze Assembly GTF with GFFCompare ---
log_message "Analyzing TACO assembly GTF with GFFCompare for accurate counting"
log_message "Note: For more accurate transcript and gene counts, refer to the gffcompare output after running script 14."

taco_gffcompare_prefix="${gffcompare_dir}/taco_assembly"
taco_counts=$(analyze_with_gffcompare "$assembly_gtf" "$taco_gffcompare_prefix")
taco_transcript_count=$(echo "$taco_counts" | cut -d':' -f1)
taco_gene_count=$(echo "$taco_counts" | cut -d':' -f2)

# Count exons directly from the GTF
taco_exons=$(grep -v "^#" "$assembly_gtf" | awk '$3=="exon"' | wc -l)

# Get statistics on merged GTF
merged_size=$(stat -c%s "$assembly_gtf")
merged_size_mb=$(echo "scale=2; $merged_size / 1048576" | bc)

log_message "TACO merged assembly statistics (from GFFCompare):"
log_message "  Size: ${merged_size_mb} MB"
log_message "  Transcripts: $taco_transcript_count"
log_message "  Genes: $taco_gene_count"
log_message "  Exons: $taco_exons"

# --- Step 8: Analyze TACO Log for Statistics ---
log_message "Analyzing TACO log for performance statistics"
taco_report=""

if [ -f "${merged_output}/taco_run.log" ]; then
    # Extract timing information
    taco_time=$(grep "Total time:" "${merged_output}/taco_run.log" | tail -1 | awk '{print $NF}')
    taco_report="  TACO reported total time: $taco_time seconds"
    
    # Extract memory information if available
    taco_memory=$(grep "Maximum resident set size" "${merged_output}/taco_run.log" | awk '{print $NF}')
    if [ -n "$taco_memory" ]; then
        taco_memory_gb=$(echo "scale=2; $taco_memory / 1048576" | bc)
        taco_report="${taco_report}\n  TACO maximum memory usage: ${taco_memory_gb} GB"
    fi
    
    # Extract other interesting metrics if available
    taco_loci=$(grep "Total loci:" "${merged_output}/taco_run.log" | awk '{print $NF}')
    if [ -n "$taco_loci" ]; then
        taco_report="${taco_report}\n  TACO total loci: $taco_loci"
    fi
fi

# Record end time and calculate duration
end_time=$(date +%s)
duration=$((end_time - start_time))
duration_formatted=$(printf '%dh:%dm:%ds' $((duration/3600)) $((duration%3600/60)) $((duration%60)))

# --- Step 9: Generate Summary Report ---
summary_file="${logs_dir}/13.taco_merge_summary.txt"
{
    echo "===== TACO Merge Summary ====="
    echo "Date: $(date)"
    echo "Execution time: $duration_formatted"
    echo
    echo "Input:"
    echo "  Mergelist file: $mergelist_file"
    echo "  GTF files merged: $mergelist_count"
    echo "  Total size of input GTFs: ${total_gtf_size_mb} MB"
    echo "  Reference GTF file: $reference_gtf"
    echo "  Reference GTF size: ${ref_size_mb} MB"
    echo "  Reference transcripts (GFFCompare): $ref_transcript_count"
    echo "  Reference genes (GFFCompare): $ref_gene_count"
    echo
    echo "Output:"
    echo "  TACO output directory: $merged_output"
    echo "  Merged assembly GTF: $assembly_gtf"
    echo "  Merged GTF size: ${merged_size_mb} MB"
    echo "  Merged transcripts (GFFCompare): $taco_transcript_count"
    echo "  Merged genes (GFFCompare): $taco_gene_count"
    echo "  Merged exons: $taco_exons"
    echo
    if [ -n "$taco_report" ]; then
        echo "TACO Performance Metrics:"
        echo -e "$taco_report"
        echo
    fi
    if [ $missing_files -gt 0 ]; then
        echo "Missing files from mergelist:"
        echo -e "$list_of_missing_files"
    fi
    echo "======================================="
} > "$summary_file"

# --- Step 10: Log Completion ---
if [ -s "$error_log" ]; then
    log_error "TACO merge process completed with errors. Check $error_log for details."
    cat "$summary_file"
    exit 1
else
    log_message "TACO merge process completed successfully in $duration_formatted"
    cat "$summary_file"
    exit 0
fi
