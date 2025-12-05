#!/bin/bash
# RNA-seq Pipeline Configuration File
# Customize this file for your specific HPC environment and project

#=====================================================================
# ENVIRONMENT SETTINGS
#=====================================================================

# Base project directory
export PROJECT_DIR="/project/igb_fargo/cpb_diapause_rnaseq"

# Reference genome settings
export REFERENCE_DIR="/project/igb_fargo/genome_references/Ldec_2.0"

# Base genome information - customize for your specific genome
export GENOME_VERSION="GCF_000500325.1_Ldec_2.0"
export GENOME_PREFIX="Ldec_2.0"  # Used for splice sites and exon files

# Genome filenames - these will be constructed from the above variables
export GENOME_FASTA_FILENAME="${GENOME_VERSION}_genomic.fna"
export GENOME_GTF_FILENAME="${GENOME_VERSION}_genomic.gtf"
export GENOME_FEATURE_FILENAME="${GENOME_VERSION}_feature_table.txt"
export SPLICE_SITES_FILENAME="${GENOME_PREFIX}.ss"
export EXONS_FILENAME="${GENOME_PREFIX}.exon"
export GENOME_INDEX_PREFIX="${GENOME_PREFIX}_index"

# FTP path for genome download - updated with proper protocol
export GENOME_FTP_BASE="ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/500/325/${GENOME_VERSION}/"

# SRA Download settings (for step 01)
export BIOPROJECT="PRJNA553565"
export NCBI_API_KEY="4a00e356c07511cb0bce047d61e116178d08"

#=====================================================================
# PIPELINE BEHAVIOR OPTIONS
#=====================================================================
# Control how the pipeline runs

# Set to true to skip SRA download and use existing FASTQ files
export SKIP_SRA_DOWNLOAD=false

# Set to true to skip genome references download and use existing reference files
export SKIP_GENOME_DOWNLOAD=true

# Choose the transcript assembly merging method: "taco" or "stringtie"
export MERGE_METHOD="taco"

# Enable or disable novel transcript discovery in StringTie
export ENABLE_NOVEL_TRANSCRIPT_DISCOVERY=true

# HPC-specific settings
export PARTITION="short"  # SLURM partition/queue to use
export DEFAULT_MEM="16G"  # Default memory request
export HIGH_MEM="64G"     # High memory request for demanding steps

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
export GFFCOMPARE_EXEC="/project/igb_fargo/programs/gffcompare/gffcompare"
export TACO_EXEC="/project/igb_fargo/programs/taco_run"
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
    export MERGED_GTF="${ASSEMBLIES_DIR}/stringtie_merged.gtf"
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
