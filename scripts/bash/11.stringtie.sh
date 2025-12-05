#!/bin/bash

#SBATCH --time=48:00:00                  # Walltime limit (HH:MM:SS)
#SBATCH --nodes=2                        # Number of nodes
#SBATCH --ntasks-per-node=40             # 40 threads per node
#SBATCH --partition=short                # Partition type
#SBATCH --job-name="stringtie_assembly"  # Descriptive job name
#SBATCH --output="11.stringtie_logs/11.stringtie_%j.out"
#SBATCH --error="11.stringtie_logs/11.stringtie_%j.err"
#SBATCH --mem=0                          # Use all available memory

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
alignments_dir="${base_dir}/alignments_bam"  # BAM file location
assemblies_dir="${base_dir}/assemblies"      # Output directory for GTF files
gtf_file="$REFERENCE_GTF"
gffcompare_exec="$GFFCOMPARE_EXEC"  # Path to gffcompare executable from config

# Standard log setup
current_dir="$(pwd)"
logs_dir="${current_dir}/11.stringtie_logs"
sample_logs_dir="${logs_dir}/sample_logs"
gffcompare_logs_dir="${logs_dir}/gffcompare_logs"
debug_dir="${logs_dir}/debug_files"

# Define log files with consistent naming convention
error_log="${logs_dir}/11.stringtie_errors.log"
progress_log="${logs_dir}/11.stringtie_progress.log"
parallel_joblog="${logs_dir}/11.stringtie_parallel_joblog.txt"
debug_log="${logs_dir}/11.stringtie_debug.log"

# Performance tuning
threads_per_job=4        # Threads per StringTie job
max_parallel_jobs=15     # Number of parallel jobs (adjust based on available resources)

# Create required directories
mkdir -p "$logs_dir" "$sample_logs_dir" "$gffcompare_logs_dir" "$debug_dir"

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

# Check if novel transcript discovery is enabled
# Default to true if not specified through environment variable
ENABLE_NOVEL_DISCOVERY=${ENABLE_NOVEL_DISCOVERY:-true}

log_message "Novel transcript discovery mode: $ENABLE_NOVEL_DISCOVERY"

# Function to process GTF with gffcompare for accurate gene counting
count_with_gffcompare() {
    local gtf_file="$1"
    local sample_name="$2"
    local output_dir="${gffcompare_logs_dir}/${sample_name}"
    local output_prefix="${output_dir}/${sample_name}"
    
    # Create output directory
    mkdir -p "$output_dir"
    
    # Save detailed debug info
    log_debug "Running gffcompare for $sample_name with reference $gtf_file"
    log_debug "Output directory: $output_dir"
    log_debug "Output prefix: $output_prefix"
    
    # Run gffcompare with the reference GTF
    "$gffcompare_exec" -r "$gtf_file" -o "$output_prefix" "$gtf_file" &> "${output_dir}/gffcompare.log"
    gff_exit_code=$?
    
    if [ $gff_exit_code -ne 0 ]; then
        log_error "GFFCompare failed for sample: $sample_name (Exit code: $gff_exit_code)"
        # Capture the error log for debugging
        cat "${output_dir}/gffcompare.log" >> "$debug_log"
        return 1
    fi
    
    # Extract counts from the stats file
    local stats_file="${output_prefix}.stats"
    
    # Save the stats file for debugging
    if [ -f "$stats_file" ]; then
        cp "$stats_file" "${debug_dir}/${sample_name}_stats.txt"
        log_debug "Stats file saved to ${debug_dir}/${sample_name}_stats.txt"
        
        # Log the entire stats file for debugging
        log_debug "Contents of stats file for $sample_name:"
        log_debug "$(cat "$stats_file")"
        
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
log_message "Creating output directory: $assemblies_dir"
mkdir -p "$assemblies_dir"
if [ $? -ne 0 ]; then
    log_error "Failed to create directory: $assemblies_dir"
    exit 1
fi

# --- Step 4: Validate Input Files ---
log_message "Validating reference GTF file"
if [ ! -f "$gtf_file" ]; then
    log_error "Reference GTF file not found: $gtf_file"
    exit 1
else
    # Get the size of the GTF file for reference
    gtf_size=$(stat -c%s "$gtf_file")
    gtf_size_mb=$(echo "scale=2; $gtf_size / 1048576" | bc)
    log_message "Reference GTF file size: ${gtf_size_mb} MB"
    
    # Run gffcompare on reference GTF to get accurate counts
    log_message "Analyzing reference GTF with gffcompare"
    ref_counts=$(count_with_gffcompare "$gtf_file" "reference")
    ref_transcript_count=$(echo "$ref_counts" | cut -d':' -f1)
    ref_gene_count=$(echo "$ref_counts" | cut -d':' -f2)
    
    log_message "Reference GTF features from gffcompare: $ref_gene_count genes, $ref_transcript_count transcripts"
fi

log_message "Validating input BAM files"
bam_files=()
while IFS= read -r -d '' file; do
    bam_files+=("$file")
done < <(find "$alignments_dir" -name "*.bam" -not -name "*.bai" -print0)

# Check if we found any files
if [ ${#bam_files[@]} -eq 0 ]; then
    log_error "No BAM files found in directory: $alignments_dir"
    exit 1
fi

log_message "Found ${#bam_files[@]} BAM files to process"

# --- Step 5: Define StringTie Function for Parallel Execution ---
run_stringtie() {
    local bam_file="$1"
    local sample_name=$(basename "$bam_file" .bam)
    local output_gtf="${assemblies_dir}/${sample_name}.gtf"
    local coverage_file="${assemblies_dir}/${sample_name}.cov"
    local abundance_file="${assemblies_dir}/${sample_name}.abund"
    local sample_log="${sample_logs_dir}/${sample_name}.log"
    
    echo "Processing: $sample_name" > "$sample_log"
    echo "BAM file: $bam_file" >> "$sample_log"
    echo "Output GTF: $output_gtf" >> "$sample_log"
    
    # Record start time for this sample
    local sample_start=$(date +%s)
    
    # Validate BAM file
    echo "Validating BAM file..." >> "$sample_log"
    if ! samtools quickcheck -v "$bam_file" &>> "$sample_log"; then
        echo "BAM file validation failed for $sample_name" >> "$error_log"
        return 1
    else
        echo "BAM file validation successful" >> "$sample_log"
        
        # Get BAM file stats for reference
        samtools flagstat "$bam_file" > "${sample_logs_dir}/${sample_name}_flagstat.txt"
        total_reads=$(grep "in total" "${sample_logs_dir}/${sample_name}_flagstat.txt" | awk '{print $1}')
        echo "Total reads in BAM: $total_reads" >> "$sample_log"
    fi
    
    # Run StringTie with parameters based on novel discovery setting
    echo "Running StringTie with novel transcript discovery: $ENABLE_NOVEL_DISCOVERY..." >> "$sample_log"

    if [ "$ENABLE_NOVEL_DISCOVERY" = "true" ]; then
        # Run with novel transcript discovery enabled
        stringtie "$bam_file" \
            -G "$gtf_file" \
            -o "$output_gtf" \
            -p "$threads_per_job" \
            -A "$abundance_file" \
            -C "$coverage_file" \
            -v &>> "$sample_log"
    else
        # Run in reference-guided mode (disable novel transcript discovery)
        # The --conservative flag and -e option restrict to reference transcripts
        stringtie "$bam_file" \
            -G "$gtf_file" \
            -o "$output_gtf" \
            -p "$threads_per_job" \
            -A "$abundance_file" \
            -C "$coverage_file" \
            -e \
            --conservative \
            -v &>> "$sample_log"
    fi
    
    local exit_code=$?
    
    # Record end time and calculate duration for this sample
    local sample_end=$(date +%s)
    local sample_duration=$((sample_end - sample_start))
    local sample_duration_formatted=$(printf '%dm:%ds' $((sample_duration/60)) $((sample_duration%60)))
    
    if [ $exit_code -ne 0 ]; then
        echo "StringTie failed for sample: $sample_name (Exit code: $exit_code)" >> "$error_log"
        return 1
    else
        # Validate output GTF
        if [ ! -s "$output_gtf" ]; then
            echo "StringTie produced empty GTF for sample: $sample_name" >> "$error_log"
            return 1
        fi
        
        # Get transcript and gene counts using gffcompare
        echo "Running GFFCompare for accurate counting..." >> "$sample_log"
        counts=$(count_with_gffcompare "$output_gtf" "$sample_name")
        transcript_count=$(echo "$counts" | cut -d':' -f1)
        gene_count=$(echo "$counts" | cut -d':' -f2)
        
        # Count exons using grep/awk
        exon_count=$(grep -v "^#" "$output_gtf" | awk '$3=="exon"' | wc -l)
        
        # Get output file sizes
        gtf_size=$(stat -c%s "$output_gtf")
        gtf_size_kb=$(echo "scale=2; $gtf_size / 1024" | bc)
        
        cov_size=0
        abund_size=0
        if [ -f "$coverage_file" ]; then
            cov_size=$(stat -c%s "$coverage_file")
            cov_size_kb=$(echo "scale=2; $cov_size / 1024" | bc)
        fi
        
        if [ -f "$abundance_file" ]; then
            abund_size=$(stat -c%s "$abundance_file")
            abund_size_kb=$(echo "scale=2; $abund_size / 1024" | bc)
        fi
        
        echo "StringTie successfully processed sample: $sample_name in $sample_duration_formatted" >> "$progress_log"
        echo "  Transcripts (GFFCompare): $transcript_count" >> "$progress_log"
        echo "  Genes (GFFCompare): $gene_count" >> "$progress_log"
        echo "  Exons: $exon_count" >> "$progress_log"
        echo "  GTF file size: ${gtf_size_kb} KB" >> "$progress_log"
        
        # Add detailed stats to sample log
        echo "Assembly Statistics:" >> "$sample_log"
        echo "  Processing time: $sample_duration_formatted" >> "$sample_log"
        echo "  Transcripts assembled (GFFCompare): $transcript_count" >> "$sample_log"
        echo "  Genes identified (GFFCompare): $gene_count" >> "$sample_log"
        echo "  Exons detected: $exon_count" >> "$sample_log"
        echo "  GTF file size: ${gtf_size_kb} KB" >> "$sample_log"
        if [ -f "$coverage_file" ]; then
            echo "  Coverage file size: ${cov_size_kb} KB" >> "$sample_log"
        fi
        if [ -f "$abundance_file" ]; then
            echo "  Abundance file size: ${abund_size_kb} KB" >> "$sample_log"
        fi
        
        # Store counts in files for later collection
        echo "$transcript_count" > "${gffcompare_logs_dir}/${sample_name}_transcripts.txt"
        echo "$gene_count" > "${gffcompare_logs_dir}/${sample_name}_genes.txt"
        
        return 0
    fi
}

# Export all functions and variables needed by parallel processes
export -f run_stringtie count_with_gffcompare log_error log_debug
export gtf_file assemblies_dir threads_per_job error_log progress_log sample_logs_dir gffcompare_logs_dir gffcompare_exec debug_log debug_dir ENABLE_NOVEL_DISCOVERY

# --- Step 6: Run StringTie in Parallel ---
log_message "Starting parallel StringTie processing with $max_parallel_jobs jobs"

parallel -u --env gtf_file --env assemblies_dir --env threads_per_job --env error_log --env progress_log --env sample_logs_dir --env gffcompare_logs_dir --env gffcompare_exec --env count_with_gffcompare --env log_error --env debug_log --env debug_dir --env log_debug --env ENABLE_NOVEL_DISCOVERY \
    --jobs "$max_parallel_jobs" --no-notice --joblog "$parallel_joblog" run_stringtie ::: "${bam_files[@]}"

# --- Step 7: Validate Assembly Output ---
log_message "Validating StringTie output files"
total_samples=${#bam_files[@]}
successful_assemblies=0
failed_assemblies=0

for bam_file in "${bam_files[@]}"; do
    sample_name=$(basename "$bam_file" .bam)
    output_gtf="${assemblies_dir}/${sample_name}.gtf"
    
    if [ -f "$output_gtf" ] && [ -s "$output_gtf" ]; then
        # Check if file is a valid GTF with transcript_id entries
        if grep -q "transcript_id" "$output_gtf"; then
            successful_assemblies=$((successful_assemblies + 1))
        else
            log_error "GTF file for $sample_name appears to be invalid (no transcript_id entries)"
            failed_assemblies=$((failed_assemblies + 1))
        fi
    else
        log_error "GTF file for $sample_name is missing or empty"
        failed_assemblies=$((failed_assemblies + 1))
    fi
done

# --- Step 8: Calculate Assembly Statistics from GFFCompare Results ---
log_message "Calculating transcript and gene statistics from GFFCompare results"
# We'll only track min, max, and individual sample stats instead of totals
min_transcripts=999999
max_transcripts=0
min_genes=999999
max_genes=0
sample_stats=()

for bam_file in "${bam_files[@]}"; do
    sample_name=$(basename "$bam_file" .bam)
    transcript_file="${gffcompare_logs_dir}/${sample_name}_transcripts.txt"
    gene_file="${gffcompare_logs_dir}/${sample_name}_genes.txt"
    
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

# Record end time and calculate duration
end_time=$(date +%s)
duration=$((end_time - start_time))
duration_formatted=$(printf '%dh:%dm:%ds' $((duration/3600)) $((duration%3600/60)) $((duration%60)))

# --- Step 9: Generate Summary Report ---
log_message "Generating summary report"
failed_jobs=0
if [ -f "$parallel_joblog" ]; then
    failed_count=$(grep -c 'exit code: [^0]' "$parallel_joblog" 2>/dev/null || echo "0")
    if [[ "$failed_count" =~ ^[0-9]+$ ]]; then
        failed_jobs="$failed_count"
    fi
fi

summary_file="${logs_dir}/11.stringtie_summary.txt"
{
    echo "===== StringTie Assembly Summary ====="
    echo "Date: $(date)"
    echo "Execution time: $duration_formatted"
    echo
    echo "Pipeline Configuration:"
    echo "  Novel transcript discovery: $ENABLE_NOVEL_DISCOVERY"
    echo
    echo "Processing Statistics:"
    echo "  Total samples processed: $total_samples"
    echo "  Successfully assembled: $successful_assemblies"
    echo "  Failed assemblies: $failed_assemblies"
    echo "  Failed jobs (from parallel): $failed_jobs"
    echo
    echo "Reference GTF Statistics (from GFFCompare):"
    echo "  Reference genes: $ref_gene_count"
    echo "  Reference transcripts: $ref_transcript_count"
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
    if [ "$failed_assemblies" -gt 0 ]; then
        echo "Failed assemblies:"
        grep -E 'StringTie failed for sample|GTF file for .* is missing or empty|GTF file for .* appears to be invalid' "$error_log" | \
        sed -E 's/.*StringTie failed for sample: ([^ ]+).*/  - \1/' | \
        sed -E 's/.*GTF file for ([^ ]+) is.*/  - \1/' | \
        sed -E 's/.*GTF file for ([^ ]+) appears.*/  - \1/' | \
        sort | uniq
    fi
    echo
    echo "Output:"
    echo "  Assembly GTF files location: $assemblies_dir"
    echo "  Sample logs location: $sample_logs_dir"
    echo "  GFFCompare logs location: $gffcompare_logs_dir"
    echo "  Debug logs location: $debug_dir"
    echo "==============================================="
} > "$summary_file"

# --- Step 10: Log Completion ---
if [ "$failed_assemblies" -gt 0 ]; then
    log_error "StringTie assembly completed with $failed_assemblies errors. See $error_log for details."
    cat "$summary_file"
    exit 1
else
    log_message "StringTie assembly completed successfully for all $total_samples samples in $duration_formatted"
    cat "$summary_file"
    exit 0
fi
