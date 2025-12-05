# R Analysis Scripts

This directory contains 14 R scripts for statistical analysis of RNA-seq data, including differential expression, co-expression networks, and functional enrichment.

## Scripts Overview

For a complete index of all scripts and their functions, see `/SCRIPTS_INDEX.md` at the repository root.

## Analysis Workflow

```
Count Matrices (from bash pipeline)
    ↓
[00] Annotation Preparation
    ↓
[01] DESeq2 Differential Expression
    ├── Males: Diapause vs Non-diapause
    ├── Females: Diapause vs Non-diapause  
    └── Sex × Diapause Interaction
    ↓
[02] Multivariate Analysis (PCA, K-means)
    ↓
[03-08] Functional Enrichment
    ├── [03-04] GO enrichment & visualization
    ├── [05-06] KEGG enrichment & visualization
    ├── [07] Overlap analysis (Venn diagrams)
    └── [08] Overlap enrichment plots
    ↓
[09-13] WGCNA Co-expression Networks
    ├── [09] Network construction
    ├── [10] Sexual dimorphism analysis
    ├── [11] Network visualization
    ├── [12] Module GO enrichment
    └── [13] Module KEGG enrichment
    ↓
[14] Special Analyses (e.g., couchpotato isoforms)
```

## Requirements

### R Version
- R ≥ 4.4.2

### Required Packages
```R
# Differential Expression
install.packages("BiocManager")
BiocManager::install("DESeq2")

# Network Analysis
install.packages("WGCNA")

# Functional Enrichment
BiocManager::install("topGO")
BiocManager::install("clusterProfiler")

# Data Manipulation & Visualization
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("igraph")
install.packages("ggraph")
```

## Input Files

Scripts expect the following directory structure:

```
project_dir/
├── 00_Annotation/
│   ├── cpb_diapause_annotation.csv
│   ├── cpb_diapause_go_omicsbox_export.txt
│   └── cpb_diapause_kaas.txt
├── 01_DESeq2/
│   ├── transcript_count_matrix.csv
│   └── cpb_sample_information.csv
└── (output directories created by scripts)
```

## Running the Analyses

### Sequential Execution (Recommended)
```R
# Set working directory to project root
setwd("/path/to/your/project")

# Run scripts in order
source("scripts/R/00_CBP_Diapause_Annotation.R")
source("scripts/R/01_CBP_Diapause_DESeq2.R")
source("scripts/R/02_CPB_Diapause_Multivariate_Analysis.R")
# ... continue with remaining scripts
```

### Individual Script Execution
Each script is self-contained and can be run independently after required dependencies are met.

## Output Files

### Differential Expression (01)
- `Diapause_vs_NonDiapause_Males_DE_Transcripts.csv`
- `Diapause_vs_NonDiapause_Females_DE_Transcripts.csv`
- `Int_DE_Transcripts.csv`
- MA plots for each comparison

### Multivariate Analysis (02)
- PCA plots
- K-means clustering results
- Variance contribution analyses

### Functional Enrichment (03-08)
- GO enrichment tables (BP, MF, CC)
- KEGG pathway enrichment tables
- Enrichment visualization plots
- Venn diagrams

### WGCNA (09-13)
- Module assignments
- Module-trait relationships
- Network visualization plots
- Module enrichment results

## Key Analysis Parameters

### DESeq2
- Interaction model: `~ sex + status + sex:status`
- Significance threshold: padj < 0.05
- Log2 fold change shrinkage: ashr method

### WGCNA
- Soft threshold power: 6
- Minimum module size: 30
- Merge cut height: 0.25 (correlation > 0.75)
- Network type: unsigned

### Enrichment Analyses
- GO: Fisher's exact test, parentChild algorithm, p < 0.01
- KEGG: Hypergeometric test, q < 0.01

## Troubleshooting

- **Memory issues**: WGCNA requires substantial RAM (16+ GB recommended)
- **Missing packages**: Install using BiocManager for Bioconductor packages
- **File path errors**: Verify working directory and input file locations
- **Plot rendering**: Ensure graphics device is properly configured

## Custom Modifications

### Adjusting Significance Thresholds
```R
# In DESeq2 analysis (script 01)
results(dds, alpha = 0.01)  # More stringent

# In functional enrichment (scripts 03, 05)
fisher.test(p.value.cutoff = 0.001)  # More stringent
```

### Modifying WGCNA Parameters
```R
# In script 09
# Adjust soft threshold
power = 8  # Increase for stricter network

# Adjust module merging
mergeCutHeight = 0.15  # Lower for more modules
```

## Contact

For questions about the R analyses, contact Alex.Torson@usda.gov
