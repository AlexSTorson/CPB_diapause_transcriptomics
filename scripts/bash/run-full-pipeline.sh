#!/bin/bash
# RNA-seq master pipeline controller script with enhanced flexibility

# --- Display help function (defined early) ---
display_help() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  --start=STEP    Start pipeline from specified step number (1-17)"
    echo "  --help          Display this help message"
    echo ""
    echo "Step numbers:"
    echo "  1: Create log directories (01.create_log_dirs.sh)"
    echo "  2: Download SRA data (02.data_download.sh) - skipped if SKIP_SRA_DOWNLOAD=true"
    echo "  3: Download genome references (03.genome_references_download.sh) - skipped if SKIP_GENOME_DOWNLOAD=true"
    echo "  4: FastQC on raw reads (04.fastqc_pre_trimming.sh)"
    echo "  5: Quality trimming (05.fastp_quality_trim.sh)"
    echo "  6: FastQC on trimmed reads (06.fastqc_post_trimming.sh)"
    echo "  7: Extract exons and splice sites (07.exon_and_splice_sites.sh)"
    echo "  8: Build genome index (08.genome_index.sh)"
    echo "  9: Map reads to genome (09.read_mapping.sh)"
    echo "  10: Convert SAM to BAM (10.sam_to_bam.sh)"
    echo "  11: Run StringTie (11.stringtie.sh)"
    echo "  12: Create mergelist (12.create_mergelist.sh)"
    echo "  13: Run assembly merge (13.taco_merge.sh or 13.stringtie_merge.sh)"
    echo "  14: Run gffcompare (14.gffcompare.sh)"
    echo "  15: Run count estimation (15.count_estimation.sh)"
    echo "  16: Run prepDE (16.prepDE.sh)"
    echo "  17: Run MultiQC (17.multiqc.sh)"
    echo ""
    # Important: exit after displaying help
    exit 0
}

# --- Process arguments first before loading configuration ---
# Default is to start from the beginning (step 1)
start_step=1

# Parse command line options
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --start=*)
            start_step="${1#*=}"
            # Validate that start_step is a number between 1 and 17
            if ! [[ "$start_step" =~ ^[0-9]+$ ]] || [ "$start_step" -lt 1 ] || [ "$start_step" -gt 17 ]; then
                echo "ERROR: --start value must be a number between 1 and 17"
                exit 1
            fi
            ;;
        --help)
            display_help
            ;;
        *)
            echo "Unknown parameter: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
    shift
done

# --- Load the Pipeline Configuration ---
script_dir="$(dirname "$(readlink -f "$0")")"
config_file="${script_dir}/00.pipeline_config.sh"

# Verify the config file exists
if [ ! -f "$config_file" ]; then
    echo "ERROR: Configuration file not found at $config_file"
    exit 1
fi

# Load the configuration
source "$config_file"

# --- Global settings ---
pipeline_dir="$(pwd)"  # Use the current directory as the pipeline directory
log_file="${pipeline_dir}/pipeline_execution.log"
job_ids_file="${pipeline_dir}/job_ids.txt"
timestamp=$(date +"%Y%m%d_%H%M%S")

# Clear log files
> "$log_file"
> "$job_ids_file"

# --- Simple Logging function ---
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$log_file"
}

# Print pipeline configuration information
log_message "Pipeline configuration:"
log_message "  - Skip SRA download: $SKIP_SRA_DOWNLOAD"
log_message "  - Skip genome download: $SKIP_GENOME_DOWNLOAD"
log_message "  - Assembly merge method: $MERGE_METHOD"
log_message "  - Novel transcript discovery: $ENABLE_NOVEL_TRANSCRIPT_DISCOVERY"
log_message "  - Merged GTF will be: $MERGED_GTF"
log_message "  - Starting from step: $start_step"
log_message ""

# Function to check if a step should be run
should_run_step() {
    local step_num=$1
    if [ $step_num -ge $start_step ]; then
        return 0 # True in bash
    else
        return 1 # False in bash
    fi
}

# --- Step 1: Create log directories (first step) ---
if should_run_step 1; then
    log_message "Submitting job: 01_create_log_dirs"
    chmod +x "${pipeline_dir}/01.create_log_dirs.sh"
    job_info_1=$(sbatch "${pipeline_dir}/01.create_log_dirs.sh")
    job_id_1=$(echo "$job_info_1" | awk '{print $4}')
    echo "01_create_log_dirs:$job_id_1" >> "$job_ids_file"
    log_message "Job submitted: 01_create_log_dirs (Job ID: $job_id_1)"
    last_dependency="afterok:$job_id_1"
else
    log_message "Skipping job: 01_create_log_dirs (starting from step $start_step)"
    # Set a dummy dependency value for later use
    last_dependency=""
fi

# --- Step 2: Download data (depends on log directory creation) ---
# Only submit if SKIP_SRA_DOWNLOAD is false and start_step <= 2
if [ "$SKIP_SRA_DOWNLOAD" = "false" ] && should_run_step 2; then
    log_message "Submitting job: 02_data_download"
    chmod +x "${pipeline_dir}/02.data_download.sh"
    
    # Set dependency if not starting from this step
    if [ -n "$last_dependency" ]; then
        job_info_2=$(sbatch --dependency=$last_dependency "${pipeline_dir}/02.data_download.sh")
    else
        job_info_2=$(sbatch "${pipeline_dir}/02.data_download.sh")
    fi
    
    job_id_2=$(echo "$job_info_2" | awk '{print $4}')
    echo "02_data_download:$job_id_2" >> "$job_ids_file"
    log_message "Job submitted: 02_data_download (Job ID: $job_id_2)"
    # Set the dependency for later steps that need data
    data_dependency="afterok:$job_id_2"
    last_dependency=$data_dependency
else
    log_message "Skipping SRA download (SKIP_SRA_DOWNLOAD=$SKIP_SRA_DOWNLOAD, start_step=$start_step)"
    # If skipping SRA download, data dependency is either log directory or none
    if [ -n "$last_dependency" ]; then
        data_dependency=$last_dependency
    else
        data_dependency=""
    fi
fi

# --- Step 3: Download genome references (depends on log directory creation) ---
if [ "$SKIP_GENOME_DOWNLOAD" = "false" ] && should_run_step 3; then
    log_message "Submitting job: 03_genome_references_download"
    chmod +x "${pipeline_dir}/03.genome_references_download.sh"
    
    # Set dependency if not starting from this step
    if [ -n "$last_dependency" ]; then
        job_info_3=$(sbatch --dependency=$last_dependency "${pipeline_dir}/03.genome_references_download.sh")
    else
        job_info_3=$(sbatch "${pipeline_dir}/03.genome_references_download.sh")
    fi
    
    job_id_3=$(echo "$job_info_3" | awk '{print $4}')
    echo "03_genome_references_download:$job_id_3" >> "$job_ids_file"
    log_message "Job submitted: 03_genome_references_download (Job ID: $job_id_3)"
    genome_dependency="afterok:$job_id_3"
    last_dependency=$genome_dependency
else
    log_message "Skipping genome references download (SKIP_GENOME_DOWNLOAD=$SKIP_GENOME_DOWNLOAD, start_step=$start_step)"
    if [ -n "$last_dependency" ]; then
        genome_dependency=$last_dependency
    else
        genome_dependency=""
    fi
fi

# --- Step 4: Run FastQC on raw reads (depends on data availability) ---
if should_run_step 4; then
    log_message "Submitting job: 04_fastqc_pre_trimming"
    chmod +x "${pipeline_dir}/04.fastqc_pre_trimming.sh"
    
    # Set dependency if not starting from this step and if data dependency exists
    if [ -n "$data_dependency" ]; then
        job_info_4=$(sbatch --dependency=$data_dependency "${pipeline_dir}/04.fastqc_pre_trimming.sh")
    else
        job_info_4=$(sbatch "${pipeline_dir}/04.fastqc_pre_trimming.sh")
    fi
    
    job_id_4=$(echo "$job_info_4" | awk '{print $4}')
    echo "04_fastqc_pre_trimming:$job_id_4" >> "$job_ids_file"
    log_message "Job submitted: 04_fastqc_pre_trimming (Job ID: $job_id_4)"
    fastqc_pre_dependency="afterok:$job_id_4"
    last_dependency=$fastqc_pre_dependency
else
    log_message "Skipping job: 04_fastqc_pre_trimming (starting from step $start_step)"
    fastqc_pre_dependency=""
fi

# --- Step 5: Run fastp quality trimming (depends on data availability) ---
if should_run_step 5; then
    log_message "Submitting job: 05_fastp_quality_trim"
    chmod +x "${pipeline_dir}/05.fastp_quality_trim.sh"
    
    # Set dependency if not starting from this step and if data dependency exists
    if [ -n "$data_dependency" ]; then
        job_info_5=$(sbatch --dependency=$data_dependency "${pipeline_dir}/05.fastp_quality_trim.sh")
    else
        job_info_5=$(sbatch "${pipeline_dir}/05.fastp_quality_trim.sh")
    fi
    
    job_id_5=$(echo "$job_info_5" | awk '{print $4}')
    echo "05_fastp_quality_trim:$job_id_5" >> "$job_ids_file"
    log_message "Job submitted: 05_fastp_quality_trim (Job ID: $job_id_5)"
    fastp_dependency="afterok:$job_id_5"
    last_dependency=$fastp_dependency
else
    log_message "Skipping job: 05_fastp_quality_trim (starting from step $start_step)"
    fastp_dependency=""
fi

# --- Step 6: Run FastQC on trimmed reads (depends on fastp) ---
if should_run_step 6; then
    log_message "Submitting job: 06_fastqc_post_trimming"
    chmod +x "${pipeline_dir}/06.fastqc_post_trimming.sh"
    
    # Set dependency if not starting from this step and if fastp dependency exists
    if [ -n "$fastp_dependency" ]; then
        job_info_6=$(sbatch --dependency=$fastp_dependency "${pipeline_dir}/06.fastqc_post_trimming.sh")
    else
        job_info_6=$(sbatch "${pipeline_dir}/06.fastqc_post_trimming.sh")
    fi
    
    job_id_6=$(echo "$job_info_6" | awk '{print $4}')
    echo "06_fastqc_post_trimming:$job_id_6" >> "$job_ids_file"
    log_message "Job submitted: 06_fastqc_post_trimming (Job ID: $job_id_6)"
    fastqc_post_dependency="afterok:$job_id_6"
    last_dependency=$fastqc_post_dependency
else
    log_message "Skipping job: 06_fastqc_post_trimming (starting from step $start_step)"
    fastqc_post_dependency=""
fi

# --- Step 7: Extract exons and splice sites (depends on genome download) ---
if should_run_step 7; then
    log_message "Submitting job: 07_exon_splice_extraction"
    chmod +x "${pipeline_dir}/07.exon_and_splice_sites.sh"
    
    # Set dependency if not starting from this step and if genome dependency exists
    if [ -n "$genome_dependency" ]; then
        job_info_7=$(sbatch --dependency=$genome_dependency "${pipeline_dir}/07.exon_and_splice_sites.sh")
    else
        job_info_7=$(sbatch "${pipeline_dir}/07.exon_and_splice_sites.sh")
    fi
    
    job_id_7=$(echo "$job_info_7" | awk '{print $4}')
    echo "07_exon_splice_extraction:$job_id_7" >> "$job_ids_file"
    log_message "Job submitted: 07_exon_splice_extraction (Job ID: $job_id_7)"
    exon_splice_dependency="afterok:$job_id_7"
    last_dependency=$exon_splice_dependency
else
    log_message "Skipping job: 07_exon_splice_extraction (starting from step $start_step)"
    exon_splice_dependency=""
fi

# --- Step 8: Build genome index (depends on exon/splice extraction) ---
if should_run_step 8; then
    log_message "Submitting job: 08_genome_index"
    chmod +x "${pipeline_dir}/08.genome_index.sh"
    
    # Set dependency if not starting from this step and if exon/splice dependency exists
    if [ -n "$exon_splice_dependency" ]; then
        job_info_8=$(sbatch --dependency=$exon_splice_dependency "${pipeline_dir}/08.genome_index.sh")
    else
        job_info_8=$(sbatch "${pipeline_dir}/08.genome_index.sh")
    fi
    
    job_id_8=$(echo "$job_info_8" | awk '{print $4}')
    echo "08_genome_index:$job_id_8" >> "$job_ids_file"
    log_message "Job submitted: 08_genome_index (Job ID: $job_id_8)"
    index_dependency="afterok:$job_id_8"
    last_dependency=$index_dependency
else
    log_message "Skipping job: 08_genome_index (starting from step $start_step)"
    index_dependency=""
fi

# --- Step 9: Map reads to genome (depends on genome index and trimmed reads) ---
if should_run_step 9; then
    log_message "Submitting job: 09_read_mapping"
    chmod +x "${pipeline_dir}/09.read_mapping.sh"
    
    # Build dependency string based on what steps were run
    mapping_dependency=""
    if [ -n "$fastp_dependency" ] && [ -n "$index_dependency" ]; then
        mapping_dependency="${fastp_dependency},${index_dependency}"
    elif [ -n "$fastp_dependency" ]; then
        mapping_dependency="$fastp_dependency"
    elif [ -n "$index_dependency" ]; then
        mapping_dependency="$index_dependency"
    fi
    
    # Submit with or without dependency
    if [ -n "$mapping_dependency" ]; then
        job_info_9=$(sbatch --dependency=$mapping_dependency "${pipeline_dir}/09.read_mapping.sh")
    else
        job_info_9=$(sbatch "${pipeline_dir}/09.read_mapping.sh")
    fi
    
    job_id_9=$(echo "$job_info_9" | awk '{print $4}')
    echo "09_read_mapping:$job_id_9" >> "$job_ids_file"
    log_message "Job submitted: 09_read_mapping (Job ID: $job_id_9)"
    mapping_job_dependency="afterok:$job_id_9"
    last_dependency=$mapping_job_dependency
else
    log_message "Skipping job: 09_read_mapping (starting from step $start_step)"
    mapping_job_dependency=""
fi

# --- Step 10: Convert SAM to BAM (depends on read mapping) ---
if should_run_step 10; then
    log_message "Submitting job: 10_sam_to_bam"
    chmod +x "${pipeline_dir}/10.sam_to_bam.sh"
    
    # Set dependency if needed
    if [ -n "$mapping_job_dependency" ]; then
        job_info_10=$(sbatch --dependency=$mapping_job_dependency "${pipeline_dir}/10.sam_to_bam.sh")
    else
        job_info_10=$(sbatch "${pipeline_dir}/10.sam_to_bam.sh")
    fi
    
    job_id_10=$(echo "$job_info_10" | awk '{print $4}')
    echo "10_sam_to_bam:$job_id_10" >> "$job_ids_file"
    log_message "Job submitted: 10_sam_to_bam (Job ID: $job_id_10)"
    bam_dependency="afterok:$job_id_10"
    last_dependency=$bam_dependency
else
    log_message "Skipping job: 10_sam_to_bam (starting from step $start_step)"
    bam_dependency=""
fi

# --- Step 11: Run StringTie on BAM files (depends on SAM to BAM) ---
# Modify the StringTie script execution based on novel transcript discovery setting
if should_run_step 11; then
    log_message "Submitting job: 11_stringtie"
    chmod +x "${pipeline_dir}/11.stringtie.sh"
    
    # Add the novel transcript discovery flag as an environment variable for the job
    export_option=""
    if [ "$ENABLE_NOVEL_TRANSCRIPT_DISCOVERY" = "true" ]; then
        export_option="--export=ENABLE_NOVEL_DISCOVERY=true"
    else
        export_option="--export=ENABLE_NOVEL_DISCOVERY=false"
    fi
    
    # Set dependency if needed
    if [ -n "$bam_dependency" ]; then
        job_info_11=$(sbatch --dependency=$bam_dependency $export_option "${pipeline_dir}/11.stringtie.sh")
    else
        job_info_11=$(sbatch $export_option "${pipeline_dir}/11.stringtie.sh")
    fi
    
    job_id_11=$(echo "$job_info_11" | awk '{print $4}')
    echo "11_stringtie:$job_id_11" >> "$job_ids_file"
    log_message "Job submitted: 11_stringtie (Job ID: $job_id_11) with novel transcript discovery: $ENABLE_NOVEL_TRANSCRIPT_DISCOVERY"
    stringtie_dependency="afterok:$job_id_11"
    last_dependency=$stringtie_dependency
else
    log_message "Skipping job: 11_stringtie (starting from step $start_step)"
    stringtie_dependency=""
fi

# --- Step 12: Create mergelist (depends on StringTie) ---
if should_run_step 12; then
    log_message "Submitting job: 12_create_mergelist"
    chmod +x "${pipeline_dir}/12.create_mergelist.sh"
    
    # Set dependency if needed
    if [ -n "$stringtie_dependency" ]; then
        job_info_12=$(sbatch --dependency=$stringtie_dependency "${pipeline_dir}/12.create_mergelist.sh")
    else
        job_info_12=$(sbatch "${pipeline_dir}/12.create_mergelist.sh")
    fi
    
    job_id_12=$(echo "$job_info_12" | awk '{print $4}')
    echo "12_create_mergelist:$job_id_12" >> "$job_ids_file"
    log_message "Job submitted: 12_create_mergelist (Job ID: $job_id_12)"
    mergelist_dependency="afterok:$job_id_12"
    last_dependency=$mergelist_dependency
else
    log_message "Skipping job: 12_create_mergelist (starting from step $start_step)"
    mergelist_dependency=""
fi

# --- Step 13: Run assembly merge based on selected method ---
if should_run_step 13; then
    if [ "$MERGE_METHOD" = "taco" ]; then
        # Run TACO merge
        log_message "Submitting job: 13_taco_merge (selected merge method: TACO)"
        chmod +x "${pipeline_dir}/13.taco_merge.sh"
        
        # Set dependency if needed
        if [ -n "$mergelist_dependency" ]; then
            job_info_13=$(sbatch --dependency=$mergelist_dependency "${pipeline_dir}/13.taco_merge.sh")
        else
            job_info_13=$(sbatch "${pipeline_dir}/13.taco_merge.sh")
        fi
        
        job_id_13=$(echo "$job_info_13" | awk '{print $4}')
        echo "13_taco_merge:$job_id_13" >> "$job_ids_file"
        log_message "Job submitted: 13_taco_merge (Job ID: $job_id_13)"
        # Set merge job ID for later dependencies
        merge_job_id=$job_id_13
    else
        # Run StringTie merge
        log_message "Submitting job: 13_stringtie_merge (selected merge method: StringTie)"
        chmod +x "${pipeline_dir}/13.stringtie_merge.sh"
        
        # Set dependency if needed
        if [ -n "$mergelist_dependency" ]; then
            job_info_13=$(sbatch --dependency=$mergelist_dependency "${pipeline_dir}/13.stringtie_merge.sh")
        else
            job_info_13=$(sbatch "${pipeline_dir}/13.stringtie_merge.sh")
        fi
        
        job_id_13=$(echo "$job_info_13" | awk '{print $4}')
        echo "13_stringtie_merge:$job_id_13" >> "$job_ids_file"
        log_message "Job submitted: 13_stringtie_merge (Job ID: $job_id_13)"
        # Set merge job ID for later dependencies
        merge_job_id=$job_id_13
    fi
    merge_dependency="afterok:$merge_job_id"
    last_dependency=$merge_dependency
else
    log_message "Skipping job: 13_merge (starting from step $start_step)"
    merge_dependency=""
    merge_job_id=""
fi

# --- Step 14: Run gffcompare (depends on merge) ---
if should_run_step 14; then
    log_message "Submitting job: 14_gffcompare"
    chmod +x "${pipeline_dir}/14.gffcompare.sh"
    
    # Set dependency if needed
    if [ -n "$merge_dependency" ]; then
        job_info_14=$(sbatch --dependency=$merge_dependency "${pipeline_dir}/14.gffcompare.sh")
    else
        job_info_14=$(sbatch "${pipeline_dir}/14.gffcompare.sh")
    fi
    
    job_id_14=$(echo "$job_info_14" | awk '{print $4}')
    echo "14_gffcompare:$job_id_14" >> "$job_ids_file"
    log_message "Job submitted: 14_gffcompare (Job ID: $job_id_14)"
    gffcompare_dependency="afterok:$job_id_14"
    last_dependency=$gffcompare_dependency
else
    log_message "Skipping job: 14_gffcompare (starting from step $start_step)"
    gffcompare_dependency=""
fi

# --- Step 15: Run count estimation (depends on merge and BAM files) ---
if should_run_step 15; then
    log_message "Submitting job: 15_count_estimation"
    chmod +x "${pipeline_dir}/15.count_estimation.sh"
    
    # Build dependency string based on what steps were run
    count_dependency=""
    if [ -n "$bam_dependency" ] && [ -n "$merge_dependency" ]; then
        count_dependency="${bam_dependency},${merge_dependency}"
    elif [ -n "$bam_dependency" ]; then
        count_dependency="$bam_dependency"
    elif [ -n "$merge_dependency" ]; then
        count_dependency="$merge_dependency"
    fi
    
    # Submit with or without dependency
    if [ -n "$count_dependency" ]; then
        job_info_15=$(sbatch --dependency=$count_dependency "${pipeline_dir}/15.count_estimation.sh")
    else
        job_info_15=$(sbatch "${pipeline_dir}/15.count_estimation.sh")
    fi
    
    job_id_15=$(echo "$job_info_15" | awk '{print $4}')
    echo "15_count_estimation:$job_id_15" >> "$job_ids_file"
    log_message "Job submitted: 15_count_estimation (Job ID: $job_id_15)"
    count_dependency="afterok:$job_id_15"
    last_dependency=$count_dependency
else
    log_message "Skipping job: 15_count_estimation (starting from step $start_step)"
    count_dependency=""
fi

# --- Step 16: Run prepDE (depends on count estimation) ---
if should_run_step 16; then
    log_message "Submitting job: 16_prepDE"
    chmod +x "${pipeline_dir}/16.prepDE.sh"
    
    # Set dependency if needed
    if [ -n "$count_dependency" ]; then
        job_info_16=$(sbatch --dependency=$count_dependency "${pipeline_dir}/16.prepDE.sh")
    else
        job_info_16=$(sbatch "${pipeline_dir}/16.prepDE.sh")
    fi
    
    job_id_16=$(echo "$job_info_16" | awk '{print $4}')
    echo "16_prepDE:$job_id_16" >> "$job_ids_file"
    log_message "Job submitted: 16_prepDE (Job ID: $job_id_16)"
    prepde_dependency="afterok:$job_id_16"
    last_dependency=$prepde_dependency
else
    log_message "Skipping job: 16_prepDE (starting from step $start_step)"
    prepde_dependency=""
fi

# --- Step 17: Run MultiQC (depends on both FastQC steps and mapping) ---
if should_run_step 17; then
    log_message "Submitting job: 17_multiqc"
    chmod +x "${pipeline_dir}/17.multiqc.sh"
    
    # Build dependency string based on what steps were run
    multiqc_dependency=""
    dependencies=()
    
    [ -n "$fastqc_pre_dependency" ] && dependencies+=("$fastqc_pre_dependency")
    [ -n "$fastqc_post_dependency" ] && dependencies+=("$fastqc_post_dependency")
    [ -n "$mapping_job_dependency" ] && dependencies+=("$mapping_job_dependency")
    
    if [ ${#dependencies[@]} -gt 0 ]; then
        multiqc_dependency=$(IFS=, ; echo "${dependencies[*]}")
    fi
    
    # Submit with or without dependency
    if [ -n "$multiqc_dependency" ]; then
        job_info_17=$(sbatch --dependency=$multiqc_dependency "${pipeline_dir}/17.multiqc.sh")
    else
        job_info_17=$(sbatch "${pipeline_dir}/17.multiqc.sh")
    fi
    
    job_id_17=$(echo "$job_info_17" | awk '{print $4}')
    echo "17_multiqc:$job_id_17" >> "$job_ids_file"
    log_message "Job submitted: 17_multiqc (Job ID: $job_id_17)"
else
    log_message "Skipping job: 17_multiqc (starting from step $start_step)"
fi

# --- Create a monitoring script for the pipeline ---
monitor_script="${pipeline_dir}/monitor_pipeline_${timestamp}.sh"
{
    echo '#!/bin/bash'
    echo "# Pipeline monitoring script"
    echo "# Created on: $(date)"
    echo
    echo 'echo "===== RNA-Seq Pipeline Status ====="'
    echo 'echo "Status as of: $(date)"'
    echo 'echo'
    echo 'echo "Pipeline configuration:"'
    echo "echo \"  - Skip SRA download: $SKIP_SRA_DOWNLOAD\""
    echo "echo \"  - Skip genome download: $SKIP_GENOME_DOWNLOAD\""
    echo "echo \"  - Assembly merge method: $MERGE_METHOD\""
    echo "echo \"  - Novel transcript discovery: $ENABLE_NOVEL_TRANSCRIPT_DISCOVERY\""
    echo "echo \"  - Starting from step: $start_step\""
    echo 'echo'
    echo 'echo "Job Status:"'
    
    while IFS=: read -r job_name job_id; do
        echo "echo \"  - $job_name (Job ID: $job_id): \$(sacct -j $job_id --format=State --noheader | head -1)\""
    done < "$job_ids_file"
    
    echo
    echo 'echo "To view detailed job information, run: sacct -j JOBID --format=JobID,JobName,State,Elapsed,ExitCode"'
    echo 'echo "To cancel all remaining jobs, run: scancel $(cat '"$job_ids_file"' | cut -d\":\" -f2 | tr \"\\n\" \" \")"'
    echo 'echo "To restart the pipeline from a specific step, run:"'
    echo 'echo "  ./run-full-pipeline.sh --start=STEP_NUMBER"'
    echo 'echo "============================================"'
} > "$monitor_script"

chmod +x "$monitor_script"
log_message "Created pipeline monitoring script: $monitor_script"

# --- Final message ---
log_message "Pipeline has been submitted successfully!"
log_message "To monitor the pipeline status, run: $monitor_script"