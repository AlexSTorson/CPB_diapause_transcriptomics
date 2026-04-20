#!/bin/bash

#SBATCH --account=igb_fargo
#SBATCH --job-name="02_rename_reads"
#SBATCH --output="02_rename_reads.out"
#SBATCH --error="02_rename_reads.err"
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=short

# --- Variables ---
RAW_READS_DIR="/project/igb_fargo/cpb_diapause_rnaseq/raw_reads"
RENAME_CSV="/project/igb_fargo/cpb_diapause_rnaseq/scripts/cpb_sample_rename.csv"

# --- Rename each sample ---
# CSV format: LibraryName,SampleName
# tr strips Windows carriage returns and the UTF-8 BOM from the file
while IFS=',' read -r library_name sample_name; do

    # Skip the header row
    [ "$library_name" = "LibraryName" ] && continue

    mv "${RAW_READS_DIR}/${library_name}_R1.fastq.gz" \
       "${RAW_READS_DIR}/${sample_name}_R1.fastq.gz"

    mv "${RAW_READS_DIR}/${library_name}_R2.fastq.gz" \
       "${RAW_READS_DIR}/${sample_name}_R2.fastq.gz"

    echo "Renamed ${library_name} -> ${sample_name}"

done < <(tr -d '\r\xef\xbb\xbf' < "$RENAME_CSV")
