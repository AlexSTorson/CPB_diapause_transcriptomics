# Bash Pipeline Scripts

This directory contains 17 bash scripts for processing RNA-seq data from raw FASTQ files through to count matrices.

## Scripts Overview

For a complete index of all scripts and their functions, see `/SCRIPTS_INDEX.md` at the repository root.

## Pipeline Workflow

```
Raw FASTQ files (SRA)
    ↓
[02] Data Download
    ↓
[04] FastQC (pre-trimming)
    ↓
[05] fastp (quality trimming)
    ↓
[06] FastQC (post-trimming)
    ↓
[07-08] Genome indexing
    ↓
[09] HISAT2 (read mapping)
    ↓
[10] SAM to BAM conversion
    ↓
[11] StringTie (transcriptome assembly)
    ↓
[12-13] Merge assemblies (TACO or StringTie)
    ↓
[14] gffcompare (assembly evaluation)
    ↓
[15-16] Count matrix generation
    ↓
[17] MultiQC (aggregate QC report)
```

## Configuration

The pipeline is configured through `00_pipeline_config.sh`. Key settings include:

- Project directories
- Reference genome paths
- Software module paths (for HPC)
- Performance parameters (threads, memory)
- Pipeline behavior flags

**Before running:** Edit `00_pipeline_config.sh` to match your environment.

## Running the Pipeline

### Option 1: Full Pipeline
```bash
./run-full-pipeline.sh
```

### Option 2: Start from Specific Step
```bash
./run-full-pipeline.sh --start=5  # Start from quality trimming
```

### Option 3: Individual Scripts
```bash
sbatch 05_fastp_quality_trim.sh
```

## Requirements

- SLURM job scheduler
- Environment modules or conda
- Software: fastp, FastQC, HISAT2, StringTie, TACO, SAMtools, MultiQC, GNU Parallel

## Output Structure

```
project_dir/
├── raw_reads/
├── trimmed_reads/
├── alignments_bam/
├── assemblies/
├── merged_assembly/
└── count_matrices/
```

## Troubleshooting

- **Job failures**: Check individual log directories (e.g., `05_fastp_quality_trim_logs/`)
- **Path issues**: Verify all paths in `00_pipeline_config.sh`
- **Memory errors**: Adjust memory requests in individual SLURM scripts
- **Module not found**: Update `MODULE_*` variables in config file

## Contact

For questions about the pipeline, contact Alex.Torson@usda.gov
