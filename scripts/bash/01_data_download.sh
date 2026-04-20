#!/bin/bash

#SBATCH --account=igb_fargo
#SBATCH --job-name="01_data_download"
#SBATCH --output="01_data_download.out"
#SBATCH --error="01_data_download.err"
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --partition=short

# --- Variables ---
BIOPROJECT="PRJNA553565"
RAW_READS_DIR="/project/igb_fargo/cpb_diapause_rnaseq/raw_reads"

# --- Modules ---
module load sratoolkit/3.2.1
module load edirect
module load pigz

mkdir -p "$RAW_READS_DIR"

# --- Get SRR accessions and sample names from NCBI ---
# efetch runinfo CSV: column 1 = SRR accession, column 12 = library name
runinfo=$(esearch -db sra -query "$BIOPROJECT" | efetch -format runinfo)

# --- Download, convert, and rename each sample ---
while IFS=',' read -r srr sample_name; do

    # Strip hidden carriage return characters from the sample name
    sample_name=$(echo "$sample_name" | tr -d '\r')

    # Download the SRA file
    prefetch -O "$RAW_READS_DIR" "$srr"

    # Convert to paired FASTQ files
    fasterq-dump \
        --outdir "$RAW_READS_DIR" \
        --split-files \
        --threads 4 \
        "${RAW_READS_DIR}/${srr}/${srr}.sra"

    # Compress
    pigz -p 4 "${RAW_READS_DIR}/${srr}_1.fastq"
    pigz -p 4 "${RAW_READS_DIR}/${srr}_2.fastq"

    # Rename from SRR number to meaningful sample name
    mv "${RAW_READS_DIR}/${srr}_1.fastq.gz" "${RAW_READS_DIR}/${sample_name}_R1.fastq.gz"
    mv "${RAW_READS_DIR}/${srr}_2.fastq.gz" "${RAW_READS_DIR}/${sample_name}_R2.fastq.gz"

done < <(echo "$runinfo" | cut -d',' -f1,12 | grep '^SRR')
