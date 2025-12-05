#!/bin/bash

#SBATCH --job-name=prepDE                # Job name
#SBATCH --output="16.prepDE_logs/16.prepDE_%j.out"
#SBATCH --error="16.prepDE_logs/16.prepDE_%j.err"
#SBATCH --time=01:00:00                  # Time limit (hh:mm:ss)
#SBATCH --partition=short                # Partition/queue
#SBATCH --ntasks=1                       # Number of tasks
#SBATCH --cpus-per-task=4                # Number of CPU cores per task
#SBATCH --mem=16G                        # Request sufficient memory

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
prepde_script="${scripts_dir}/prepDE.py" 
ballgown_dir="${base_dir}/ballgown"
gene_count_matrix="${ballgown_dir}/gene_count_matrix.csv"
transcript_count_matrix="${ballgown_dir}/transcript_count_matrix.csv"

# Standard log setup
current_dir="$(pwd)"
logs_dir="${current_dir}/16.prepDE_logs"

# Define log files with consistent naming convention
error_log="${logs_dir}/16.prepDE_errors.log"
progress_log="${logs_dir}/16.prepDE_progress.log"
python_log="${logs_dir}/16.prepDE_python.log"

# Clear logs at the start of the script
> "$error_log"
> "$progress_log"
> "$python_log"

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
$MODULE_PYTHON  # Adjust based on your environment

# --- Step 2: Validate Input Files and Directories ---
log_message "Validating input files and directories"

# Check prepDE.py script
if [ ! -f "$prepde_script" ]; then
    log_error "prepDE.py script not found: $prepde_script"
    exit 1
fi

if [ ! -x "$prepde_script" ]; then
    log_message "Making prepDE.py script executable"
    chmod +x "$prepde_script"
    if [ $? -ne 0 ]; then
        log_error "Failed to make prepDE.py script executable"
        exit 1
    fi
fi

# Check ballgown directory
if [ ! -d "$ballgown_dir" ]; then
    log_error "Ballgown directory not found: $ballgown_dir"
    exit 1
fi

# Check if there are sample subdirectories in the ballgown directory
sample_dirs=()
while IFS= read -r -d '' dir; do
    sample_dirs+=("$dir")
done < <(find "$ballgown_dir" -mindepth 1 -maxdepth 1 -type d -not -path "$ballgown_dir/.*" -print0)

if [ ${#sample_dirs[@]} -eq 0 ]; then
    log_error "No sample directories found in ballgown directory: $ballgown_dir"
    exit 1
fi

log_message "Found ${#sample_dirs[@]} sample directories in ballgown directory"

# --- Step 3: Validate Sample Directories ---
log_message "Validating sample directories for prepDE.py input"

valid_samples=0
invalid_samples=0
for sample_dir in "${sample_dirs[@]}"; do
    sample_name=$(basename "$sample_dir")
    sample_gtf="${sample_dir}/${sample_name}.gtf"
    
    if [ ! -f "$sample_gtf" ]; then
        log_error "Sample GTF file not found: $sample_gtf"
        invalid_samples=$((invalid_samples + 1))
    else
        # Check if GTF file is valid (contains transcript_id entries)
        if grep -q "transcript_id" "$sample_gtf"; then
            valid_samples=$((valid_samples + 1))
        else
            log_error "Invalid GTF file for sample $sample_name (no transcript_id entries)"
            invalid_samples=$((invalid_samples + 1))
        fi
    fi
done

if [ $valid_samples -eq 0 ]; then
    log_error "No valid sample GTF files found for prepDE.py"
    exit 1
fi

log_message "Validated $valid_samples sample directories (found $invalid_samples invalid samples)"

# --- Step 4: Create a Sample GTF Map File ---
log_message "Creating sample GTF map file"
gtf_map_file="${logs_dir}/gtf_map.txt"
> "$gtf_map_file"  # Clear the file if it exists

for sample_dir in "${sample_dirs[@]}"; do
    sample_name=$(basename "$sample_dir")
    sample_gtf="${sample_dir}/${sample_name}.gtf"
    
    if [ -f "$sample_gtf" ] && grep -q "transcript_id" "$sample_gtf"; then
        echo "${sample_name} ${sample_gtf}" >> "$gtf_map_file"
    fi
done

# Count lines in the GTF map file
map_entries=$(wc -l < "$gtf_map_file")
log_message "Created GTF map file with $map_entries entries"

if [ $map_entries -eq 0 ]; then
    log_error "No valid GTF files found for any samples"
    exit 1
fi

# --- Step 5: Run prepDE.py Script ---
log_message "Running prepDE.py script"
log_message "Command: python $prepde_script -i $ballgown_dir -g $gene_count_matrix -t $transcript_count_matrix -l 150"

# Run prepDE.py with explicit options
python "$prepde_script" -i "$ballgown_dir" -g "$gene_count_matrix" -t "$transcript_count_matrix" -l 150 &> "$python_log"

# Check if prepDE.py ran successfully
if [ $? -ne 0 ]; then
    log_error "prepDE.py script failed. Check $python_log for details."
    exit 1
else
    log_message "prepDE.py script completed successfully"
fi

# --- Step 6: Validate Output Files ---
log_message "Validating output count matrices"

# Check gene count matrix
if [ ! -f "$gene_count_matrix" ]; then
    log_error "Gene count matrix not created: $gene_count_matrix"
    exit 1
fi

if [ ! -s "$gene_count_matrix" ]; then
    log_error "Gene count matrix is empty: $gene_count_matrix"
    exit 1
fi

# Check transcript count matrix
if [ ! -f "$transcript_count_matrix" ]; then
    log_error "Transcript count matrix not created: $transcript_count_matrix"
    exit 1
fi

if [ ! -s "$transcript_count_matrix" ]; then
    log_error "Transcript count matrix is empty: $transcript_count_matrix"
    exit 1
fi

# Count rows and columns in count matrices
gene_rows=$(wc -l < "$gene_count_matrix")
gene_rows=$((gene_rows - 1))  # Subtract header line
gene_cols=$(head -n 1 "$gene_count_matrix" | tr ',' '\n' | wc -l)
gene_cols=$((gene_cols - 1))  # Subtract gene ID column

transcript_rows=$(wc -l < "$transcript_count_matrix")
transcript_rows=$((transcript_rows - 1))  # Subtract header line
transcript_cols=$(head -n 1 "$transcript_count_matrix" | tr ',' '\n' | wc -l)
transcript_cols=$((transcript_cols - 1))  # Subtract transcript ID column

log_message "Gene count matrix contains:"
log_message "  $gene_rows genes"
log_message "  $gene_cols samples"

log_message "Transcript count matrix contains:"
log_message "  $transcript_rows transcripts"
log_message "  $transcript_cols samples"

# Check for consistent sample counts
if [ $gene_cols -ne $transcript_cols ]; then
    log_error "Inconsistent sample counts between gene ($gene_cols) and transcript ($transcript_cols) matrices"
fi

# Record end time and calculate duration
end_time=$(date +%s)
duration=$((end_time - start_time))
duration_formatted=$(printf '%dm:%ds' $((duration/60)) $((duration%60)))

# --- Step 7: Generate Summary Report ---
summary_file="${logs_dir}/16.prepDE_summary.txt"
{
    echo "===== prepDE.py Summary ====="
    echo "Date: $(date)"
    echo "Execution time: $duration_formatted"
    echo
    echo "Input:"
    echo "  Ballgown directory: $ballgown_dir"
    echo "  Total sample directories: ${#sample_dirs[@]}"
    echo "  Valid sample GTFs: $valid_samples"
    echo "  Invalid sample GTFs: $invalid_samples"
    echo
    echo "Output:"
    echo "  Gene count matrix: $gene_count_matrix"
    echo "    Genes: $gene_rows"
    echo "    Samples: $gene_cols"
    echo "  Transcript count matrix: $transcript_count_matrix"
    echo "    Transcripts: $transcript_rows"
    echo "    Samples: $transcript_cols"
    echo
    if [ $invalid_samples -gt 0 ]; then
        echo "Warning: $invalid_samples sample directories had invalid or missing GTF files"
    fi
    if [ $gene_cols -ne $transcript_cols ]; then
        echo "Warning: Inconsistent sample counts between gene and transcript matrices"
    fi
    echo "======================================="
} > "$summary_file"

# --- Step 8: Log Completion ---
if [ $invalid_samples -gt 0 ]; then
    log_error "prepDE.py completed with warnings. $invalid_samples sample directories had invalid or missing GTF files."
    cat "$summary_file"
    # Continue with exit code 0 as this is not a fatal error
    exit 0
elif [ -s "$error_log" ]; then
    log_error "prepDE.py completed with errors. Check $error_log for details."
    cat "$summary_file"
    exit 1
else
    log_message "prepDE.py completed successfully in $duration_formatted"
    cat "$summary_file"
    exit 0
fi
