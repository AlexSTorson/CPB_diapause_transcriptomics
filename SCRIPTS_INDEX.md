# Scripts Index

This file provides an overview of all scripts in this repository.

## Bash Pipeline Scripts (17 scripts)

Located in `scripts/bash/`

### Configuration and Setup
- **00_pipeline_config.sh** - Master configuration file with all pipeline parameters
- **01_create_log_dirs.sh** - Creates log directories for all pipeline steps
- **refresh_config.sh** - Helper script to refresh configuration for new projects
- **rel_config_paths.sh** - Updates configuration file paths across all scripts  
- **run-full-pipeline.sh** - Master script to run the complete bash pipeline

### Data Acquisition
- **02_data_download.sh** - Downloads raw FASTQ files from NCBI SRA
- **03_genome_references_download.sh** - Downloads reference genome and annotations

### Quality Control and Preprocessing
- **04_fastqc_pre_trimming.sh** - Quality assessment of raw reads
- **05_fastp_quality_trim.sh** - Quality trimming and adapter removal
- **06_fastqc_post_trimming.sh** - Quality assessment of trimmed reads

### Read Mapping
- **07_exon_and_splice_sites.sh** - Extracts exons and splice sites from genome
- **08_genome_index.sh** - Builds HISAT2 genome index
- **09_read_mapping.sh** - Maps trimmed reads to reference genome
- **10_sam_to_bam.sh** - Converts SAM files to sorted, indexed BAM files

### Transcriptome Assembly
- **11_stringtie.sh** - Assembles transcriptomes for each sample
- **12_create_mergelist.sh** - Creates list of GTF files for merging
- **13_stringtie_merge.sh** - Merges transcriptomes using StringTie (alternative to TACO)
- **13_taco_merge.sh** - Merges transcriptomes using TACO (alternative to StringTie merge)
- **14_gffcompare.sh** - Compares assembled transcriptome to reference annotation

### Count Matrix Generation
- **15_count_estimation.sh** - Re-estimates transcript abundances with merged GTF
- **16_prepDE.sh** - Generates count matrices for DESeq2 analysis
- **17_multiqc.sh** - Aggregates QC reports from all pipeline steps

---

## R Analysis Scripts (14 scripts)

Located in `scripts/R/`

### Core Analyses
- **00_CBP_Diapause_Annotation.R** - Prepares annotation data for downstream analyses
- **01_CBP_Diapause_DESeq2.R** - Differential expression analysis with interaction model
- **02_CPB_Diapause_Multivariate_Analysis.R** - PCA and K-means clustering
- **09_CPB_Diapause_WGCNA.R** - Co-expression network analysis

### Functional Enrichment
- **03_CBP_Diapause_DE_topGO.R** - Gene Ontology enrichment for DE categories
- **04_CBP_Diapause_DE_topGO_Summary_Plots.R** - Visualization of GO enrichment results
- **05_CPB_Diapause_DE_KEGG_Enrichment.R** - KEGG pathway enrichment for DE categories
- **06_CPB_Diapause_DE_KEGG_Summary.R** - Visualization of KEGG enrichment results
- **12_CPB_Diapause_WGCNA_GO_Enrichment_Analysis.R** - GO enrichment for WGCNA modules
- **13_CPB_Diapause_WGCNA_KEGG_Enrichment_Analysis.R** - KEGG enrichment for WGCNA modules

### Comparative and Overlap Analyses
- **07_CBP_Diapause_DE_Overlap_Venn.R** - Venn diagrams of DE overlaps between sexes
- **08_CBP_Diapause_DE_Overlap_Enrichment_Summary_Plots.R** - Enrichment plots for overlaps

### WGCNA Visualization and Analysis
- **10_CPB_Diapause_WGCNA_Dimorphism_Analysis.R** - Analysis of sexual dimorphism in modules
- **11_CPB_Diapause_WGCNA_Figures.R** - Network visualization and module plots

### Special Analyses
- **14_CBP_Diapause_Couchpotato.R** - Detailed analysis of couchpotato isoforms

---

## Python Helper Scripts

Located in `scripts/bash/` (used by bash pipeline)

- **hisat2_extract_exons.py** - Extracts exon coordinates from GTF (used by script 07)
- **hisat2_extract_splice_sites.py** - Extracts splice site coordinates (used by script 07)
- **prepDE.py** - Extracts count matrices from StringTie output (used by script 16)

---

## Usage Notes

### Running the Bash Pipeline
The bash pipeline is designed to run on HPC systems with SLURM job scheduling. Each script can be run independently or as part of the full pipeline using `run-full-pipeline.sh`.

**To run the complete pipeline:**
```bash
cd scripts/bash
./run-full-pipeline.sh
```

**To run from a specific step:**
```bash
./run-full-pipeline.sh --start=5  # Start from step 5 (quality trimming)
```

### Running R Analyses
R scripts should be run sequentially after the bash pipeline completes. Scripts are numbered to indicate the recommended execution order.

**Typical workflow:**
```R
# 1. Prepare annotations
source("scripts/R/00_CBP_Diapause_Annotation.R")

# 2. Differential expression
source("scripts/R/01_CBP_Diapause_DESeq2.R")

# 3. Multivariate analysis
source("scripts/R/02_CPB_Diapause_Multivariate_Analysis.R")

# 4-8. Functional enrichment
# (Run sequentially)

# 9-11. WGCNA
# (Run sequentially)
```

---

## Dependencies

See README.md for complete software requirements.
