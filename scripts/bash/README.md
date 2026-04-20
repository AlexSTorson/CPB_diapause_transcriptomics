# Bash Pipeline Scripts

This directory contains 15 bash scripts and 3 Python helper scripts for processing RNA-seq data from raw SRA downloads through to count matrices ready for differential expression analysis.

## Complete Pipeline Overview

The pipeline follows this sequential workflow:

```
SRA BioProject PRJNA553565
    ↓
01_data_download.sh          # Download raw FASTQ files from SRA
    ↓
02_rename_reads.sh           # Rename files using sample metadata
    ↓
03_genome_download.sh        # Download L. decemlineata reference genome
    ↓
04_fastqc_pre.sh            # Quality control (pre-trimming)
    ↓
05_fastp.sh                 # Quality trimming and adapter removal
    ↓
06_fastqc_post.sh           # Quality control (post-trimming)
    ↓
07_extract_splice_exons.sh  # Extract splice sites/exons from reference GTF
    ↓
08_genome_index.sh          # Build HISAT2 splice-aware genome index
    ↓
09_read_mapping.sh          # Map reads to genome with HISAT2
    ↓
10_sam_to_bam.sh           # Convert SAM to sorted, indexed BAM
    ↓
11_stringtie.sh            # Transcript assembly per sample
    ↓
12_create_mergelist.sh     # Prepare file list for merging
    ↓
13_taco_merge.sh           # Merge assemblies across samples (TACO)
    ↓
14_count_estimation.sh     # Re-quantify against merged transcriptome
    ↓
15_prepDE.sh               # Generate count matrices for DESeq2
    ↓
Count matrices ready for downstream analysis in R
```

## Script Details

### Data Acquisition (01-03)
- **01_data_download.sh**: Downloads all 22 RNA-seq samples from NCBI SRA BioProject PRJNA553565
- **02_rename_reads.sh**: Renames files from SRR accessions to meaningful sample names using `cpb_sample_rename.csv`
- **03_genome_download.sh**: Downloads *Leptinotarsa decemlineata* RefSeq assembly GCF_000500325.1_Ldec_2.0

### Quality Control (04-06)
- **04_fastqc_pre.sh**: FastQC quality reports on raw reads
- **05_fastp.sh**: Adapter removal, quality trimming, and filtering
- **06_fastqc_post.sh**: FastQC quality reports on trimmed reads

### Genome Preparation (07-08)
- **07_extract_splice_exons.sh**: Extracts splice sites and exons from reference GTF using HISAT2 helper scripts
- **08_genome_index.sh**: Builds splice-aware HISAT2 genome index

### Read Mapping (09-10)
- **09_read_mapping.sh**: Maps trimmed reads to genome using HISAT2 with `--dta` flag for transcript assembly
- **10_sam_to_bam.sh**: Converts SAM to sorted, indexed BAM files

### Transcript Assembly (11-13)
- **11_stringtie.sh**: Assembles transcripts per sample guided by reference GTF
- **12_create_mergelist.sh**: Creates file list of individual assemblies
- **13_taco_merge.sh**: Merges assemblies across samples using TACO

### Quantification (14-15)
- **14_count_estimation.sh**: Re-quantifies all samples against merged transcriptome using StringTie `-e` mode
- **15_prepDE.sh**: Converts StringTie output to gene and transcript count matrices for DESeq2

## Python Helper Scripts

- **prepDE.py3**: Converts StringTie Ballgown output to DESeq2-compatible count matrices
- **hisat2_extract_splice_sites.py**: Extracts splice junction coordinates from GTF files
- **hisat2_extract_exons.py**: Extracts exon coordinates from GTF files

## Configuration Files

- **cpb_sample_rename.csv**: Maps SRA library names to meaningful sample names (LibraryName,SampleName format)

## Key Software Requirements

- **SLURM**: Job scheduler (all scripts use SLURM directives)
- **FastQC** 0.11.9: Quality control reports
- **fastp**: Adapter trimming and quality filtering
- **HISAT2** 2.2.1: Splice-aware read mapping
- **StringTie** 2.2.0: Transcript assembly and quantification
- **TACO**: Multi-sample transcript merging
- **SAMtools** 1.17: SAM/BAM manipulation
- **Python** 3.11.1: For helper scripts

## System Configuration

All scripts are configured for the USDA-ARS SCINet computing environment:
- Account: `igb_fargo`
- Partitions: `short` (most scripts), with appropriate time/memory limits
- Module system for software loading
- Project directory: `/project/igb_fargo/cpb_diapause_rnaseq/`

## Output Directory Structure

```
/project/igb_fargo/cpb_diapause_rnaseq/
├── raw_reads/                    # Downloaded FASTQ files
├── raw_reads_fastqc/             # Pre-trimming QC reports
├── trimmed_reads/                # Quality-trimmed FASTQ files
├── trimmed_reads_fastqc/         # Post-trimming QC reports
├── trimming_summaries/           # fastp HTML/JSON reports
├── alignments_sam/               # HISAT2 SAM output
├── alignments_bam/               # Sorted, indexed BAM files
├── assemblies/                   # Per-sample StringTie assemblies
│   ├── taco_assembly_merge/      # TACO merged assembly
│   └── mergelist.txt            # List of input assemblies
├── ballgown/                     # StringTie quantification output
│   └── [sample_dirs]/           # Per-sample Ballgown directories
└── scripts/                      # This directory
```

## Running the Pipeline

### Sequential Execution
Run scripts in numerical order (01-15), checking for successful completion before proceeding:

```bash
# Submit jobs individually
sbatch 01_data_download.sh
# Wait for completion, then:
sbatch 02_rename_reads.sh
# Continue sequentially...
```

### Dependencies
- Scripts 04-06 can run in parallel after script 02
- Script 07 requires script 03 completion
- Scripts 08 requires script 07 completion  
- Scripts 09+ must run sequentially as written

## Expected Runtime

Estimated wallclock times on SCINet (22 samples, ~150M read pairs):
- Data download: ~12 hours
- Quality trimming: ~8 hours  
- Read mapping: ~12 hours
- Transcript assembly: ~8 hours
- TACO merging: ~4 hours
- Count estimation: ~10 hours

**Total pipeline runtime**: ~2-3 days

## Final Output

The pipeline produces two key files for downstream analysis:
- `ballgown/gene_count_matrix.csv`: Gene-level counts for DESeq2
- `ballgown/transcript_count_matrix.csv`: Transcript-level counts for DESeq2

These matrices contain raw read counts for all 22 samples across detected genes/transcripts, ready for import into R for differential expression analysis.

## Troubleshooting

- **Job failures**: Check `.out` and `.err` files for each script
- **Path issues**: Verify project directory structure matches expected layout
- **Module loading**: Ensure all required software modules are available on your system
- **Memory errors**: Increase `--mem` values for memory-intensive steps
- **Storage space**: Monitor disk usage, especially during BAM file generation

For additional details on the complete analysis pipeline, see the R scripts in `../R/` and the main repository documentation.
