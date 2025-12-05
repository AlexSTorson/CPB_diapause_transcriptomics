#!/bin/bash

#SBATCH --time=12:00:00                 # Walltime limit (HH:MM:SS)
#SBATCH --nodes=2                       # Two nodes for parallel processing
#SBATCH --ntasks-per-node=40            # 40 tasks per node
#SBATCH --partition=short               # Partition type
#SBATCH --job-name="count_estimation"   # Descriptive job name
#SBATCH --output="15.count_estimation_logs/15.count_estimation_%j.out"
#SBATCH --error="15.count_estimation_logs/15.count_estimation_%j.err"
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
base_dir="$PROJECT_DIR"
alignments_dir="${base_dir}/alignments_bam"           # BAM file location
ballgown_dir="${base_dir}/ballgown"                   # Output directory for Ballgown
gtf_file="$MERGED_GTF"                                # Use merged GTF from config
gffcompare_exec="$GFFCOMPARE_EXEC"                    # Path to gffcompare executable from config

# Standard log setup
current_dir="$(pwd)"
logs_dir="${current_dir}/15.count_estimation_logs"
sample_logs_dir="${logs_dir}/sample_logs"
gffcompare_dir="${logs_dir}/gffcompare_output"
debug_dir="${logs_dir}/debug_files"

# Define log files with consistent naming convention
error_log="${logs_dir}/15.count_estimation_errors.log"
progress_log="${logs_dir}/15.count_estimation_progress.log"
parallel_joblog="${logs_dir}/15.count_estimation_parallel_joblog.txt"
debug_log="${logs_dir}/15.count_estimation_debug.log"

# Performance tuning
threads_per_job=4       # Threads for each StringTie job
max_parallel_jobs=15    # Number of parallel jobs (adjust based on available resources)

# Create required directories
mkdir -p "$logs_dir" "$sample_logs_dir" "$gffcompare_dir" "$debug_dir"

# Clear logs at the start of the script
> "$error_log"
> "$progress_log"
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
count_with_gffcompare() {
    local gtf_file="$1"
    local sample_name="$2"
    local output_dir="${gffcompare_dir}/${sample_name}"
    local output_prefix="${output_dir}/${sample_name}"
    local gffcompare_log_file="${output_dir}/gffcompare.log"
    
    # Create output directory
    mkdir -p "$output_dir"
    
    # Save detailed debug info
    log_debug "Running gffcompare for $sample_name"
    log_debug "Output directory: $output_dir"
    log_debug "Output prefix: $output_prefix"
    
    # Run gffcompare with the reference GTF
    "$gffcompare_exec" -r "$gtf_file" -o "$output_prefix" "$gtf_file" &> "$gffcompare_log_file"
    gff_exit_code=$?
    
    if [ $gff_exit_code -ne 0 ]; then
        log_error "GFFCompare failed for sample: $sample_name (Exit code: $gff_exit_code)"
        # Capture the error log for debugging
        cat "$gffcompare_log_file" >> "$debug_log"
        return 1
    fi
    
    # Extract counts from the stats file
    local stats_file="${output_prefix}.stats"
    
    # Save the stats file for debugging
    if [ -f "$stats_file" ]; then
        cp "$stats_file" "${debug_dir}/${sample_name}_stats.txt"
        log_debug "Stats file saved to ${debug_dir}/${sample_name}_stats.txt"
        
        # Try multiple patterns for transcript counts
        transcript_count=""
        
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
        
        # Try multiple patterns for gene counts
        gene_count=""
        
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
        
        # Fallback method for older gffcompare versions
        if [ -z "$gene_count" ]; then
            log_debug "Using alternative method to count genes from tmap file"
            
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
            log_debug "Attempting to count genes directly from GTF loci attribute"
            
            # Try counting from the annotated GTF
            local annotated_gtf="${output_prefix}.annotated.gtf"
            if [ -f "$annotated_gtf" ]; then
                # Extract gene_id or locus ID and count unique values
                gene_count=$(grep -v "^#" "$annotated_gtf" | grep -o 'locus "[^"]*"' | sort -u | wc -l)
                log_debug "Gene count from annotated GTF locus attributes: $gene_count"
                
                # If that fails, try gene_id
                if [ -z "$gene_count" ] || [ "$gene_count" -eq 0 ]; then
                    gene_count=$(grep -v "^#" "$annotated_gtf" | grep -o 'gene_id "[^"]*"' | sort -u | wc -l)
                    log_debug "Gene count from annotated GTF gene_id attributes: $gene_count"
                fi
            else
                log_debug "No annotated GTF file found at $annotated_gtf"
            fi
        fi
        
        # Default to transcript count / 1.5 as a rough estimate if all else fails
        if [ -z "$gene_count" ] || [ "$gene_count" -eq 0 ]; then
            if [ -n "$transcript_count" ] && [ "$transcript_count" -gt 0 ]; then
                gene_count=$(echo "$transcript_count * 0.66" | bc | awk '{printf("%d\n",$1 + 0.5)}')
                log_debug "Using estimated gene count (transcript/1.5): $gene_count"
            else
                gene_count=0
                log_debug "Setting gene count to 0 as fallback"
            fi
        fi
        
        # Ensure default values if still empty
        if [ -z "$transcript_count" ]; then
            transcript_count=0
            log_debug "Setting transcript count to 0 as fallback"
        fi
        
        if [ -z "$gene_count" ]; then
            gene_count=0
            log_debug "Setting gene count to 0 as fallback"
        fi
        
        echo "${transcript_count}:${gene_count}"
    else
        log_error "Stats file not found for $sample_name: $stats_file"
        log_debug "Stats file not found: $stats_file"
        echo "0:0"
    fi
}

# Record start time for performance tracking
start_time=$(date +%s)

# --- Step 1: Load Necessary Modules ---
log_message "Loading required modules"
$MODULE_STRINGTIE  # Load the StringTie module
$MODULE_PARALLEL   # Load GNU Parallel for job management
$MODULE_SAMTOOLS   # For BAM validation

# --- Step 2: Validate Required Executables ---
log_message "Validating required executables"

# Check gffcompare executable
if [ ! -f "$gffcompare_exec" ]; then
    log_error "GFFCompare executable not found: $gffcompare_exec"
    exit 1
fi

if [ ! -x "$gffcompare_exec" ]; then
    log_error "GFFCompare executable is not executable: $gffcompare_exec"
    exit 1
fi

log_message "GFFCompare executable found: $gffcompare_exec"

# --- Step 3: Create Output Directory ---
log_message "Creating output directory: $ballgown_dir"
mkdir -p "$ballgown_dir"
if [ $? -ne 0 ]; then
    log_error "Failed to create Ballgown directory: $ballgown_dir"
    exit 1
fi

# --- Step 4: Validate Input Files ---
log_message "Validating input files"
log_message "Using merge method: $MERGE_METHOD"
log_message "Merged GTF file: $gtf_file"

# Check merged GTF file
if [ ! -f "$gtf_file" ]; then
    log_error "Merged GTF file not found: $gtf_file"
    exit 1
fi

if [ ! -s "$gtf_file" ]; then
    log_error "Merged GTF file is empty: $gtf_file"
    exit 1
fi

# Get reference GTF stats for logging using gffcompare for accurate counts
log_message "Analyzing GTF file with gffcompare for accurate counts"
ref_counts=$(count_with_gffcompare "$gtf_file" "reference")
ref_transcript_count=$(echo "$ref_counts" | cut -d':' -f1)
ref_gene_count=$(echo "$ref_counts" | cut -d':' -f2)

# Get file sizes for logging
gtf_size=$(stat -c%s "$gtf_file")
gtf_size_mb=$(echo "scale=2; $gtf_size / 1048576" | bc)

log_message "Merged GTF file: $gtf_file"
log_message "  Size: ${gtf_size_mb} MB"
log_message "  Transcripts (from GFFCompare): $ref_transcript_count"
log_message "  Genes (from GFFCompare): $ref_gene_count"

# Check alignments directory
if [ ! -d "$alignments_dir" ]; then
    log_error "Alignments directory not found: $alignments_dir"
    exit 1
fi

# Find BAM files
log_message "Locating BAM files for count estimation"
bam_files=()
while IFS= read -r -d '' file; do
    bam_files+=("$file")
done < <(find "$alignments_dir" -name "*.bam" -not -name "*.bai" -print0)

# Check if we found any BAM files
if [ ${#bam_files[@]} -eq 0 ]; then
    log_error "No BAM files found in directory: $alignments_dir"
    exit 1
fi

log_message "Found ${#bam_files[@]} BAM files to process"

# --- Step 5: Define Count Estimation Function ---
run_stringtie_counts() {
    local bam_file="$1"
    local sample_name=$(basename "$bam_file" .bam)
    local sample_dir="${ballgown_dir}/${sample_name}"
    local sample_log="${sample_logs_dir}/${sample_name}.log"
    
    echo "Processing: $sample_name" > "$sample_log"
    echo "BAM file: $bam_file" >> "$sample_log"
    
    # Record start time for this sample
    local sample_start=$(date +%s)
    
    # Validate BAM file
    echo "Validating BAM file..." >> "$sample_log"
    if ! samtools quickcheck -v "$bam_file" &>> "$sample_log"; then
        echo "BAM file validation failed for sample: $sample_name" >> "$error_log"
        return 1
    fi
    
    # Create sample directory
    mkdir -p "$sample_dir"
    if [ $? -ne 0 ]; then
        echo "Failed to create sample directory: $sample_dir" >> "$error_log"
        return 1
    fi
    
    # Run StringTie with the specified parameters for Ballgown
    echo "Running StringTie for count estimation..." >> "$sample_log"
    stringtie "$bam_file" \
        -e \
        -B \
        -p "$threads_per_job" \
        -G "$gtf_file" \
        -o "${sample_dir}/${sample_name}.gtf" \
        -A "${sample_dir}/${sample_name}.gene_abund.tab" \
        -C "${sample_dir}/${sample_name}.cov_refs.gtf" \
        -v &>> "$sample_log"
    
    local exit_code=$?
    
    # Record end time and calculate duration for this sample
    local sample_end=$(date +%s)
    local sample_duration=$((sample_end - sample_start))
    local sample_duration_formatted=$(printf '%dm:%ds' $((sample_duration/60)) $((sample_duration%60)))
    
    if [ $exit_code -ne 0 ]; then
        echo "StringTie count estimation failed for sample: $sample_name (Exit code: $exit_code)" >> "$error_log"
        return 1
    else
        # Validate output files
        echo "Validating output files..." >> "$sample_log"
        if [ ! -f "${sample_dir}/${sample_name}.gtf" ]; then
            echo "StringTie did not produce a GTF file for sample: $sample_name" >> "$error_log"
            return 1
        fi
        
        # Get transcript and gene counts using GFFCompare
        echo "Running GFFCompare for accurate counting..." >> "$sample_log"
        local sample_gtf="${sample_dir}/${sample_name}.gtf"
        local counts=$(count_with_gffcompare "$sample_gtf" "$sample_name")
        local transcript_count=$(echo "$counts" | cut -d':' -f1)
        local gene_count=$(echo "$counts" | cut -d':' -f2)
        
        # Check for Ballgown files
        local ballgown_files=(
            "${sample_dir}/e2t.ctab"
            "${sample_dir}/e_data.ctab"
            "${sample_dir}/i2t.ctab"
            "${sample_dir}/i_data.ctab"
            "${sample_dir}/t_data.ctab"
        )
        
        local missing_bg_files=0
        for bg_file in "${ballgown_files[@]}"; do
            if [ ! -f "$bg_file" ]; then
                echo "Missing Ballgown file: $bg_file" >> "$sample_log"
                missing_bg_files=$((missing_bg_files + 1))
            fi
        done
        
        if [ $missing_bg_files -gt 0 ]; then
            echo "StringTie did not produce all required Ballgown files for sample: $sample_name (Missing: $missing_bg_files)" >> "$error_log"
            # Don't fail, this is expected behavior in some StringTie versions
        fi
        
        # Get output file sizes for reporting
        local gtf_size=$(stat -c%s "${sample_dir}/${sample_name}.gtf")
        local gtf_size_kb=$(echo "scale=2; $gtf_size / 1024" | bc)
        
        echo "StringTie count estimation completed for sample: $sample_name in $sample_duration_formatted" >> "$progress_log"
        echo "  Transcripts (GFFCompare): $transcript_count" >> "$progress_log"
        echo "  Genes (GFFCompare): $gene_count" >> "$progress_log"
        echo "  Output GTF size: ${gtf_size_kb} KB" >> "$progress_log"
        
        # Add detailed stats to sample log
        echo "Count Estimation Statistics:" >> "$sample_log"
        echo "  Processing time: $sample_duration_formatted" >> "$sample_log"
        echo "  Transcripts (GFFCompare): $transcript_count" >> "$sample_log"
        echo "  Genes (GFFCompare): $gene_count" >> "$sample_log"
        echo "  Output GTF size: ${gtf_size_kb} KB" >> "$sample_log"
        
        # Store counts in files for later collection
        echo "$transcript_count" > "${gffcompare_dir}/${sample_name}_transcripts.txt"
        echo "$gene_count" > "${gffcompare_dir}/${sample_name}_genes.txt"
        
        return 0
    fi
}

export -f run_stringtie_counts count_with_gffcompare log_debug
export ballgown_dir gtf_file threads_per_job error_log progress_log logs_dir sample_logs_dir gffcompare_dir gffcompare_exec debug_log debug_dir

# --- Step 6: Perform Count Estimation for Ballgown in Parallel ---
log_message "Starting parallel count estimation with $max_parallel_jobs jobs"

parallel -u --env ballgown_dir --env gtf_file --env threads_per_job --env error_log --env progress_log --env logs_dir --env sample_logs_dir --env count_with_gffcompare --env gffcompare_dir --env gffcompare_exec --env log_debug --env debug_log --env debug_dir \
    --jobs "$max_parallel_jobs" --no-notice --joblog "$parallel_joblog" run_stringtie_counts ::: "${bam_files[@]}"

# --- Step 7: Validate the Ballgown Directory Structure ---
log_message "Validating Ballgown directory structure"
total_samples=${#bam_files[@]}
successful_samples=0
failed_samples=0
missing_bg_files=0

for bam_file in "${bam_files[@]}"; do
    sample_name=$(basename "$bam_file" .bam)
    sample_dir="${ballgown_dir}/${sample_name}"
    
    if [ -d "$sample_dir" ] && [ -f "${sample_dir}/${sample_name}.gtf" ]; then
        # Check if the sample has the required Ballgown files
        ballgown_files=(
            "${sample_dir}/e2t.ctab"
            "${sample_dir}/e_data.ctab"
            "${sample_dir}/i2t.ctab"
            "${sample_dir}/i_data.ctab"
            "${sample_dir}/t_data.ctab"
        )
        
        sample_missing=0
        for bg_file in "${ballgown_files[@]}"; do
            if [ ! -f "$bg_file" ]; then
                sample_missing=$((sample_missing + 1))
            fi
        done
        
        successful_samples=$((successful_samples + 1))
        
        if [ $sample_missing -gt 0 ]; then
            log_error "Sample $sample_name is missing $sample_missing Ballgown files"
            missing_bg_files=$((missing_bg_files + sample_missing))
        fi
    else
        failed_samples=$((failed_samples + 1))
        log_error "Sample directory or GTF file missing for sample: $sample_name"
    fi
done

log_message "Ballgown directory validation complete:"
log_message "  Total samples: $total_samples"
log_message "  Successful samples: $successful_samples"
log_message "  Failed samples: $failed_samples"
log_message "  Missing Ballgown files: $missing_bg_files"

# --- Step 8: Calculate Assembly Statistics from GFFCompare Results ---
log_message "Calculating transcript and gene statistics from GFFCompare results"
# We'll track min, max, and individual sample stats
min_transcripts=999999
max_transcripts=0
min_genes=999999
max_genes=0
sample_stats=()

for bam_file in "${bam_files[@]}"; do
    sample_name=$(basename "$bam_file" .bam)
    transcript_file="${gffcompare_dir}/${sample_name}_transcripts.txt"
    gene_file="${gffcompare_dir}/${sample_name}_genes.txt"
    
    # Print debug info about the files
    log_debug "Checking count files for sample $sample_name"
    log_debug "Transcript file: $transcript_file (exists: $(test -f "$transcript_file" && echo "yes" || echo "no"))"
    log_debug "Gene file: $gene_file (exists: $(test -f "$gene_file" && echo "yes" || echo "no"))"
    
    transcript_count=0
    gene_count=0
    
    if [ -f "$transcript_file" ]; then
        transcript_count=$(cat "$transcript_file")
        log_debug "Read transcript count for $sample_name: $transcript_count"
    else
        log_debug "Transcript file missing for $sample_name"
    fi
    
    if [ -f "$gene_file" ]; then
        gene_count=$(cat "$gene_file")
        log_debug "Read gene count for $sample_name: $gene_count"
    else
        log_debug "Gene file missing for $sample_name"
    fi
    
    log_message "Sample $sample_name counts: $transcript_count transcripts, $gene_count genes"
    
    # Ensure we have valid numbers
    if [[ "$transcript_count" =~ ^[0-9]+$ ]]; then
        # Update transcript min/max
        if [ "$transcript_count" -lt "$min_transcripts" ]; then min_transcripts=$transcript_count; fi
        if [ "$transcript_count" -gt "$max_transcripts" ]; then max_transcripts=$transcript_count; fi
    else
        log_error "Invalid transcript count for sample $sample_name: $transcript_count"
    fi
    
    if [[ "$gene_count" =~ ^[0-9]+$ ]]; then
        # Update gene min/max
        if [ "$gene_count" -lt "$min_genes" ]; then min_genes=$gene_count; fi
        if [ "$gene_count" -gt "$max_genes" ]; then max_genes=$gene_count; fi
    else
        log_error "Invalid gene count for sample $sample_name: $gene_count"
    fi
    
    # Add to sample stats for reporting
    sample_stats+=("  $sample_name: $transcript_count transcripts, $gene_count genes")
done

# If min values were not set, adjust them to 0
if [ "$min_transcripts" -eq 999999 ]; then min_transcripts=0; fi
if [ "$min_genes" -eq 999999 ]; then min_genes=0; fi

# --- Step 9: Create Sample List File ---
log_message "Creating sample list file for downstream analysis"

# Create a file with all sample names for use in R analysis
sample_list_file="${ballgown_dir}/sample_list.txt"
> "$sample_list_file"  # Clear the file if it exists

for bam_file in "${bam_files[@]}"; do
    sample_name=$(basename "$bam_file" .bam)
    sample_dir="${ballgown_dir}/${sample_name}"
    
    if [ -d "$sample_dir" ] && [ -f "${sample_dir}/${sample_name}.gtf" ]; then
        echo "$sample_name" >> "$sample_list_file"
    fi
done

sample_count=$(wc -l < "$sample_list_file")
log_message "Created sample list file: $sample_list_file with $sample_count entries"

# Record end time and calculate duration
end_time=$(date +%s)
duration=$((end_time - start_time))
duration_formatted=$(printf '%dh:%dm:%ds' $((duration/3600)) $((duration%3600/60)) $((duration%60)))

# --- Step 10: Generate Summary Report ---
log_message "Generating summary report"
failed_jobs=$(grep -c 'exit code: [^0]' "$parallel_joblog" 2>/dev/null || echo 0)

summary_file="${logs_dir}/15.count_estimation_summary.txt"
{
    echo "===== Count Estimation for Ballgown Summary ====="
    echo "Date: $(date)"
    echo "Execution time: $duration_formatted"
    echo
    echo "Pipeline Configuration:"
    echo "  Merge method: $MERGE_METHOD"
    echo "  Merged GTF file: $gtf_file"
    echo
    echo "Input:"
    echo "  Merged GTF file: $gtf_file"
    echo "  GTF size: ${gtf_size_mb} MB"
    echo "  Reference transcripts (GFFCompare): $ref_transcript_count"
    echo "  Reference genes (GFFCompare): $ref_gene_count"
    echo "  BAM files processed: $total_samples"
    echo
    echo "Output:"
    echo "  Ballgown directory: $ballgown_dir"
    echo "  Sample list file: $sample_list_file with $sample_count samples"
    echo
    echo "Processing Statistics:"
    echo "  Successfully processed samples: $successful_samples"
    echo "  Failed samples: $failed_samples"
    echo "  Failed jobs (from parallel): $failed_jobs"
    echo "  Samples missing some Ballgown files: $missing_bg_files"
    echo
    echo "Transcript Statistics (from GFFCompare):"
    echo "  Minimum transcripts per sample: $min_transcripts"
    echo "  Maximum transcripts per sample: $max_transcripts"
    echo
    echo "Gene Statistics (from GFFCompare):"
    echo "  Minimum genes per sample: $min_genes"
    echo "  Maximum genes per sample: $max_genes"
    echo
    echo "Sample-specific Statistics:"
    for stat in "${sample_stats[@]}"; do
        echo "$stat"
    done
    echo
    if [ $failed_samples -gt 0 ]; then
        echo "Failed samples:"
        grep 'Sample directory or GTF file missing for sample' "$error_log" | sed 's/.*missing for sample: /  - /'
    fi
    echo "======================================="
} > "$summary_file"

# --- Step 11: Log Completion ---
if [ $failed_samples -gt 0 ]; then
    log_error "Count estimation completed with $failed_samples failed samples. See $error_log for details."
    cat "$summary_file"
    exit 1
else
    log_message "Count estimation completed successfully for all $total_samples samples in $duration_formatted"
    cat "$summary_file"
    exit 0
fi
