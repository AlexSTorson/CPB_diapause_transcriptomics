#!/bin/bash

#SBATCH --time=04:00:00                 # Walltime limit (HH:MM:SS)
#SBATCH --nodes=1                       # Single node
#SBATCH --ntasks=1                      # Single task (merging is a single process)
#SBATCH --cpus-per-task=8               # Use 8 threads for StringTie
#SBATCH --partition=short               # Partition type
#SBATCH --job-name="stringtie_merge"    # Descriptive job name
#SBATCH --output="13.stringtie_merge_logs/13.stringtie_merge_%j.out"
#SBATCH --error="13.stringtie_merge_logs/13.stringtie_merge_%j.err"

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
mergelist_file="${assemblies_dir}/mergelist.txt"  # Mergelist file
merged_gtf_file="${assemblies_dir}/stringtie_merged.gtf"  # Output GTF file
reference_gtf="$REFERENCE_GTF"
gffcompare_exec="$GFFCOMPARE_EXEC"  # Path to gffcompare executable from config

# Standard log setup
current_dir="$(pwd)"
logs_dir="${current_dir}/13.stringtie_merge_logs"
gffcompare_dir="${logs_dir}/gffcompare_output"

# Define log files with consistent naming convention
error_log="${logs_dir}/13.stringtie_merge_errors.log"
progress_log="${logs_dir}/13.stringtie_merge_progress.log"
stringtie_log="${logs_dir}/13.stringtie_merge_output.log"
gffcompare_log="${logs_dir}/13.stringtie_merge_gffcompare.log"
debug_log="${logs_dir}/13.stringtie_merge_debug.log"

# Performance tuning
threads=8  # Number of threads for StringTie -merge

# Create necessary directories
mkdir -p "$logs_dir" "$gffcompare_dir"

# Clear logs at the start of the script
> "$error_log"
> "$progress_log"
> "$stringtie_log"
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
    # This is similar to the approach in script 11
    
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

# --- Step 1: Load Necessary Modules ---
log_message "Loading required modules"
$MODULE_STRINGTIE

# --- Step 2: Validate Input Files ---
log_message "Validating input files"

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

# Run GFFCompare on reference GTF to get accurate counts
log_message "Analyzing reference GTF with GFFCompare"
ref_gffcompare_prefix="${gffcompare_dir}/reference"
ref_counts=$(analyze_with_gffcompare "$reference_gtf" "$ref_gffcompare_prefix")
ref_transcript_count=$(echo "$ref_counts" | cut -d':' -f1)
ref_gene_count=$(echo "$ref_counts" | cut -d':' -f2)

log_message "Reference GTF file: $reference_gtf"
log_message "  Size: ${ref_size_mb} MB"
log_message "  Transcripts (from GFFCompare): $ref_transcript_count"
log_message "  Genes (from GFFCompare): $ref_gene_count"

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
done < "$mergelist_file"

if [ $missing_files -gt 0 ]; then
    log_error "Found $missing_files missing GTF files in mergelist"
    exit 1
fi

total_gtf_size_mb=$(echo "scale=2; $total_gtf_size / 1048576" | bc)
log_message "Total size of all GTF files to merge: ${total_gtf_size_mb} MB"

# --- Step 3: Create Backup of Existing Merged GTF (if present) ---
if [ -f "$merged_gtf_file" ]; then
    backup_file="${merged_gtf_file}.backup_$(date +%Y%m%d_%H%M%S)"
    log_message "Creating backup of existing merged GTF file: $backup_file"
    cp "$merged_gtf_file" "$backup_file"
    if [ $? -ne 0 ]; then
        log_error "Failed to create backup of existing merged GTF file"
    fi
fi

# --- Step 4: Perform the Merge ---
log_message "Starting StringTie merge process with $threads threads"
log_message "Command: stringtie --merge -p $threads -G $reference_gtf -o $merged_gtf_file $mergelist_file"

# Run StringTie merge
stringtie --merge -p $threads -G "$reference_gtf" -o "$merged_gtf_file" "$mergelist_file" &> "$stringtie_log"

if [ $? -ne 0 ]; then
    log_error "StringTie merge failed. Check $stringtie_log for details."
    exit 1
else
    log_message "StringTie merge command completed successfully"
fi

# --- Step 5: Validate Merged GTF ---
log_message "Validating merged GTF file"

if [ ! -f "$merged_gtf_file" ]; then
    log_error "Merged GTF file was not created: $merged_gtf_file"
    exit 1
fi

if [ ! -s "$merged_gtf_file" ]; then
    log_error "Merged GTF file is empty: $merged_gtf_file"
    exit 1
fi

# --- Step 6: Analyze Merged GTF with GFFCompare ---
log_message "Analyzing merged GTF with GFFCompare for accurate counting"
merge_gffcompare_prefix="${gffcompare_dir}/merged"
merge_counts=$(analyze_with_gffcompare "$merged_gtf_file" "$merge_gffcompare_prefix")
merge_transcript_count=$(echo "$merge_counts" | cut -d':' -f1)
merge_gene_count=$(echo "$merge_counts" | cut -d':' -f2)

# Count exons directly from the GTF
merged_exons=$(grep -v "^#" "$merged_gtf_file" | awk '$3=="exon"' | wc -l)

# Get statistics on merged GTF
merged_size=$(stat -c%s "$merged_gtf_file")
merged_size_mb=$(echo "scale=2; $merged_size / 1048576" | bc)

log_message "Merged GTF file statistics from GFFCompare:"
log_message "  Size: ${merged_size_mb} MB"
log_message "  Transcripts: $merge_transcript_count"
log_message "  Genes: $merge_gene_count"
log_message "  Exons: $merged_exons"

# Record end time and calculate duration
end_time=$(date +%s)
duration=$((end_time - start_time))
duration_formatted=$(printf '%dh:%dm:%ds' $((duration/3600)) $((duration%3600/60)) $((duration%60)))

# --- Step 7: Generate Summary Report ---
summary_file="${logs_dir}/13.stringtie_merge_summary.txt"
{
    echo "===== StringTie Merge Summary ====="
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
    echo "  Merged GTF file: $merged_gtf_file"
    echo "  Merged GTF size: ${merged_size_mb} MB"
    echo "  Merged transcripts (GFFCompare): $merge_transcript_count"
    echo "  Merged genes (GFFCompare): $merge_gene_count"
    echo "  Merged exons: $merged_exons"
    echo
    if [ $missing_files -gt 0 ]; then
        echo "Missing files from mergelist:"
        echo -e "$list_of_missing_files"
    fi
    echo "======================================="
} > "$summary_file"

# --- Step 8: Log Completion ---
if [ -s "$error_log" ]; then
    log_error "StringTie merge process completed with errors. Check $error_log for details."
    cat "$summary_file"
    exit 1
else
    log_message "StringTie merge process completed successfully in $duration_formatted"
    cat "$summary_file"
    exit 0
fi
