#!/bin/bash

# Script to refresh the pipeline configuration file for a new project
# This should be run when starting a new RNA-seq analysis project

# Log file
log_file="refresh_config_$(date +%Y%m%d_%H%M%S).log"
echo "Starting config refresh at $(date)" > "$log_file"

# Target config file
config_file="./00.pipeline_config.sh"

# Check if the config file exists
if [ ! -f "$config_file" ]; then
    echo "Error: $config_file not found! Make sure you're running this from the scripts directory." | tee -a "$log_file"
    exit 1
fi

echo "Processing $config_file..." | tee -a "$log_file"

# Create backup of the original file
backup_file="${config_file}.bak_$(date +%Y%m%d_%H%M%S)"
cp "$config_file" "$backup_file"
echo "Created backup: $backup_file" | tee -a "$log_file"

# Create the refreshed configuration file
refreshed_config=$(cat << 'EOL'
#!/bin/bash
# RNA-seq Pipeline Configuration File
# Customize this file for your specific HPC environment and project

#=====================================================================
# ENVIRONMENT SETTINGS
#=====================================================================

# Base project directory - CHANGE THIS FOR NEW PROJECT
export PROJECT_DIR="/path/to/your/project"

# Reference genome settings - CUSTOMIZE FOR YOUR GENOME
export REFERENCE_DIR="/path/to/genome_references/species"

# Base genome information - customize for your specific genome
export GENOME_VERSION="genome_version"      # Example: "GCF_000500325.1_Ldec_2.0"
export GENOME_PREFIX="genome_prefix"        # Example: "Ldec_2.0"

# Genome filenames - these will be constructed from the above variables
export GENOME_FASTA_FILENAME="${GENOME_VERSION}_genomic.fna"
export GENOME_GTF_FILENAME="${GENOME_VERSION}_genomic.gtf"
export GENOME_FEATURE_FILENAME="${GENOME_VERSION}_feature_table.txt"
export SPLICE_SITES_FILENAME="${GENOME_PREFIX}.ss"
export EXONS_FILENAME="${GENOME_PREFIX}.exon"
export GENOME_INDEX_PREFIX="${GENOME_PREFIX}_index"

# FTP path for genome download - customize for your genome
export GENOME_FTP_BASE="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/path/to/your/genome/"

# SRA Download settings (for step 01) - CUSTOMIZE FOR YOUR PROJECT
export BIOPROJECT="PRJNAXXXXXX"
export NCBI_API_KEY="your_api_key_here"

# HPC-specific settings
export PARTITION="short"  # SLURM partition/queue to use
export DEFAULT_MEM="16G"  # Default memory request
export HIGH_MEM="64G"     # High memory request for demanding steps

#=====================================================================
# PIPELINE BEHAVIOR OPTIONS
#=====================================================================
# Control how the pipeline runs

# Set to true to skip SRA download and use existing FASTQ files
export SKIP_SRA_DOWNLOAD=true

# Choose the transcript assembly merging method: "taco" or "stringtie"
export MERGE_METHOD="stringtie"

# Enable or disable novel transcript discovery in StringTie
export ENABLE_NOVEL_TRANSCRIPT_DISCOVERY=false

#=====================================================================
# MODULE SETTINGS
#=====================================================================
# Customize these based on the modules available in your HPC environment

# Module load commands - set to empty string if not using module system
export MODULE_FASTQC="module load fastqc/0.11.9"
export MODULE_PARALLEL="module load parallel"
export MODULE_SAMTOOLS="module load samtools/1.17"
export MODULE_PYTHON="module load python_3/3.11.1"
export MODULE_HISAT2="module load hisat2/2.2.1"
export MODULE_STRINGTIE="module load stringtie/2.2.0"
export MODULE_MULTIQC="module load multiqc/1.15"
export MODULE_SRA="module load sratoolkit/3.2.0"
export MODULE_EDIRECT="module load edirect"
export MODULE_PIGZ="module load pigz"
export MODULE_MINICONDA="module load miniconda"

# Path to executables that might not be available as modules
export GFFCOMPARE_EXEC="/path/to/gffcompare"
export TACO_EXEC="/path/to/taco_run"
export PREPDE_SCRIPT="${PROJECT_DIR}/scripts/prepDE.py"

#=====================================================================
# CONDA ENVIRONMENTS
#=====================================================================
# Set names of conda environments for tools that require them

export CONDA_ENV_FASTP="fastpEnv"  # Conda environment for fastp

#=====================================================================
# DIRECTORY STRUCTURE
#=====================================================================
# These settings control the directory structure - usually no need to modify

# Main project directories
export RAW_READS_DIR="${PROJECT_DIR}/raw_reads"
export SRA_DIR="${PROJECT_DIR}/sra_files"
export SRA_CACHE="${PROJECT_DIR}/sra_cache"
export RAW_FASTQC_DIR="${PROJECT_DIR}/raw_reads_fastqc_results"
export TRIMMED_READS_DIR="${PROJECT_DIR}/trimmed_reads"
export TRIMMING_SUMMARIES_DIR="${PROJECT_DIR}/trimming_summaries"
export TRIMMED_FASTQC_DIR="${PROJECT_DIR}/trimmed_reads_fastqc_results"
export ALIGNMENTS_SAM_DIR="${PROJECT_DIR}/alignments_sam"
export ALIGNMENTS_BAM_DIR="${PROJECT_DIR}/alignments_bam"
export ALIGNMENTS_BAM_INDEXED_DIR="${PROJECT_DIR}/alignments_bam_indexed"
export ASSEMBLIES_DIR="${PROJECT_DIR}/assemblies"
export TACO_ASSEMBLY_DIR="${ASSEMBLIES_DIR}/taco_assembly_merge"
export STRINGTIE_ASSEMBLY_DIR="${ASSEMBLIES_DIR}/stringtie_merged"
export GFFCOMPARE_OUTPUT_DIR="${PROJECT_DIR}/gffcompare_output"
export BALLGOWN_DIR="${PROJECT_DIR}/ballgown"
export MULTIQC_OUTPUT_DIR="${PROJECT_DIR}/multiqc_report"
export GENOME_INDEX_DIR="${REFERENCE_DIR}/genome_index"

# Constructed paths for reference files
export REFERENCE_GENOME="${REFERENCE_DIR}/${GENOME_FASTA_FILENAME}"
export REFERENCE_GTF="${REFERENCE_DIR}/${GENOME_GTF_FILENAME}"
export REFERENCE_FEATURE="${REFERENCE_DIR}/${GENOME_FEATURE_FILENAME}"
export SPLICE_SITES_FILE="${REFERENCE_DIR}/${SPLICE_SITES_FILENAME}"
export EXONS_FILE="${REFERENCE_DIR}/${EXONS_FILENAME}"
export GENOME_INDEX="${GENOME_INDEX_DIR}/${GENOME_INDEX_PREFIX}"

# Full file paths for FTP downloads
export GENOME_FTP_FASTA="${GENOME_FTP_BASE}${GENOME_FASTA_FILENAME}.gz"
export GENOME_FTP_GTF="${GENOME_FTP_BASE}${GENOME_GTF_FILENAME}.gz"
export GENOME_FTP_FEATURE="${GENOME_FTP_BASE}${GENOME_FEATURE_FILENAME}.gz"

# Define the merged GTF file path based on selected merge method
if [ "$MERGE_METHOD" = "taco" ]; then
    export MERGED_GTF="${TACO_ASSEMBLY_DIR}/assembly.gtf"
else
    export MERGED_GTF="${STRINGTIE_ASSEMBLY_DIR}/stringtie_merged.gtf"
fi

#=====================================================================
# PERFORMANCE SETTINGS
#=====================================================================
# Adjust these based on your HPC's capabilities and your dataset size

# Default threads for processing
export THREADS_PER_FASTQC=2
export THREADS_PER_FASTP=4
export THREADS_PER_HISAT=4
export THREADS_PER_SAMTOOLS=4
export THREADS_PER_STRINGTIE=4
export THREADS_PER_TACO=20

# Maximum parallel jobs 
export MAX_PARALLEL_FASTQC=20
export MAX_PARALLEL_FASTP=10
export MAX_PARALLEL_HISAT=15
export MAX_PARALLEL_SAMTOOLS=15
export MAX_PARALLEL_STRINGTIE=15
EOL
)

# Write the refreshed config to the file
echo "$refreshed_config" > "$config_file"

echo "Refreshed configuration file successfully" | tee -a "$log_file"

# Now, update the config path in all pipeline scripts
echo "Updating config file path in all pipeline scripts..." | tee -a "$log_file"

# Template path to use in all scripts
NEW_PATH="/project/igb_fargo/cpb_diapause_rnaseq_v2/scripts/00.pipeline_config.sh"

# Find all shell scripts in the current directory
SCRIPTS=$(find . -maxdepth 1 -name "*.sh" | sort)

# Count of modified files
MODIFIED=0

# Common patterns to look for in the config file path sections
patterns=(
    "/project/igb_fargo/cpb_diapause_rnaseq_v2/scripts/00.pipeline_config.sh"
    "/90daydata/igb_fargo/cpb_diapause_rnaseq_lebenzon/scripts/00.pipeline_config.sh"
    "/project/igb_fargo/cpb_diapause_rnaseq_v2/scripts/00.pipeline_config.sh"
)

for script in $SCRIPTS; do
    # Skip this refresh script itself if it's in the same directory
    if [[ "$(basename "$script")" == "$(basename "$0")" ]]; then
        continue
    fi
    
    # Skip the config file itself
    if [[ "$(basename "$script")" == "00.pipeline_config.sh" ]]; then
        continue
    fi
    
    echo "Processing $script..." | tee -a "$log_file"
    
    # Create backup of original file
    cp "$script" "${script}.bak"
    echo "  Created backup: ${script}.bak" | tee -a "$log_file"
    
    # Flag to track if we've updated this file
    updated=false
    
    # First try exact pattern matches
    for pattern in "${patterns[@]}"; do
        if grep -q "$pattern" "$script"; then
            sed -i "s|$pattern|$NEW_PATH|g" "$script"
            echo "  Updated config path in $script using pattern: $pattern" | tee -a "$log_file"
            updated=true
            break
        fi
    done
    
    # If none of the specific patterns matched, try a more generic approach
    if [ "$updated" = false ]; then
        # Look for lines that contain both "config_file" and "pipeline_config.sh"
        if grep -q "config_file.*00.pipeline_config.sh" "$script"; then
            # Use sed to replace the line containing config file path with template path
            sed -i "/config_file.*00.pipeline_config.sh/c\config_file=\"$NEW_PATH\"" "$script"
            echo "  Updated config path in $script using generic pattern" | tee -a "$log_file"
            updated=true
        else
            echo "  No matching config path found in $script" | tee -a "$log_file"
        fi
    fi
    
    # Increment counter if updated
    if [ "$updated" = true ]; then
        MODIFIED=$((MODIFIED + 1))
    fi
done

echo "Update completed at $(date)" | tee -a "$log_file"
echo "Modified $MODIFIED script files. See $log_file for details." | tee -a "$log_file"

echo
echo "You now need to edit $config_file and set the following:" | tee -a "$log_file"
echo "1. PROJECT_DIR - Path to your new project directory" | tee -a "$log_file"
echo "2. REFERENCE_DIR - Path to your genome references" | tee -a "$log_file"
echo "3. GENOME_VERSION and GENOME_PREFIX - Set according to your genome" | tee -a "$log_file"
echo "4. GENOME_FTP_BASE - Update with the correct FTP path" | tee -a "$log_file"
echo "5. BIOPROJECT - Set to your project's BioProject accession" | tee -a "$log_file"
echo "6. NCBI_API_KEY - Update with your API key" | tee -a "$log_file"
echo "7. Verify all paths to executables (GFFCOMPARE_EXEC, TACO_EXEC, etc.)" | tee -a "$log_file"