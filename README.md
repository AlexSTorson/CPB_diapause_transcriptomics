# Sex-Specific Transcriptomics of Diapause in Colorado Potato Beetle

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

## Overview

This repository contains analysis code and supplemental data for:

**Diapause in the Colorado potato beetle is characterized by sex-specific gene expression**

*Alex S. Torson, Dacotah Melicher, George D. Yocum, Joseph P. Rinehart*

USDA-ARS Edward T. Schafer Agricultural Research Center, Biosciences Research Laboratory, Fargo, ND 58102, USA

Published in *Journal of Insect Physiology* (2025)

## Abstract

This study characterizes sex-specific transcriptomic profiles of diapausing and non-diapausing adult Colorado potato beetles (*Leptinotarsa decemlineata*) using RNA-seq. We identified 10,060 differentially expressed transcripts, with distinct patterns of shared responses (2,366 transcripts), female-specific responses (4,653 transcripts), male-specific responses (2,813 transcripts), and sex×diapause interactions (228 transcripts). Our findings reveal that diapause involves differential regulation of core sex-independent processes alongside robust sex-specific regulatory mechanisms.

## Repository Structure

```
CPB_diapause_transcriptomics/
├── scripts/
│   ├── bash/              # Shell scripts for read processing and mapping (17 scripts)
│   └── R/                 # R scripts for statistical analyses (14 scripts)
├── data/
│   └── supplemental_tables/  # Supplemental data tables (S1-S18)
├── docs/
│   ├── WORKFLOW.md        # Detailed workflow documentation
│   └── BASH_SCRIPTS.md    # Bash pipeline documentation
│   └── R_SCRIPTS.md       # R analysis documentation
├── README.md              # This file
└── LICENSE                # MIT License
```

## Data Availability

### Raw Sequencing Data
- **NCBI BioProject**: PRJNA553565
- **SRA Accessions**: SRX6428284-SRX6428305
- 21 RNA-seq libraries (150 bp paired-end, Illumina NextSeq 500)

### Reference Genome
- **Species**: *Leptinotarsa decemlineata*
- **RefSeq Assembly**: GCF_000500325.1
- **BioProject**: PRJNA171749

## Quick Start

### 1. Clone Repository
```bash
git clone https://github.com/YOUR-USERNAME/CPB_diapause_transcriptomics.git
cd CPB_diapause_transcriptomics
```

### 2. Bash Pipeline (Read Processing)
See `docs/BASH_SCRIPTS.md` for detailed documentation of all 17 bash scripts:
- Data download and quality control
- Read mapping with HISAT2
- Transcriptome assembly with StringTie/TACO
- Count matrix generation

### 3. R Analysis Pipeline
See `docs/R_SCRIPTS.md` for detailed documentation of all 14 R scripts:
- Differential expression analysis (DESeq2)
- Co-expression network analysis (WGCNA)
- Functional enrichment (GO, KEGG)
- Data visualization

## Key Findings

### Core Diapause Responses (Shared between sexes)
- Suppression of oxidative phosphorylation (49 transcripts)
- Enhanced stress tolerance (HIF-1α: 24.3-fold; L(2)efl: 31.7-fold)
- Metabolic reorganization and substrate conservation
- Cytoskeletal stabilization

### Female-Specific Responses
- Enhanced insulin signaling pathway regulation
- Spliceosome pathway enrichment (28 transcripts)
- Increased lipid mobilization capacity
- Suppression of circadian clock components

### Male-Specific Responses
- mRNA surveillance pathway activation (15 transcripts)
- Enhanced DNA repair mechanisms
- Post-transcriptional control emphasis
- Genomic maintenance systems

### Sex×Diapause Interactions
- Opposing growth control strategies (Hippo, TGF-β)
- Differential signaling network modulation
- Sex-specific precision in regulatory responses

## Software Requirements

### Bash Pipeline
- fastp ≥ 1.0.1
- FastQC ≥ 0.11.9
- HISAT2 ≥ 2.2.1
- StringTie ≥ 1.3.4
- TACO ≥ 0.7.3
- SAMtools ≥ 1.x
- MultiQC ≥ 1.7
- GNU Parallel

### R Analysis
- R ≥ 4.4.2
- DESeq2 ≥ 1.34.0
- WGCNA ≥ 1.73
- topGO ≥ 2.58.0
- clusterProfiler ≥ 4.14.6
- tidyverse
- ggplot2
- igraph
- ggraph

## Supplemental Tables

All supplemental tables (S1-S18) are provided in `data/supplemental_tables/`:

| Table | Description | Size |
|-------|-------------|------|
| **S1** | Manuscript transcript reference table | 6 KB |
| **S2** | Complete CPB diapause annotation | 5.8 MB |
| **S3** | GO annotations (OmicsBox export) | 1.2 MB |
| **S4** | KEGG annotations (KAAS) | 474 KB |
| **S5** | K-means clustering results | 19 KB |
| **S6** | Male diapause vs non-diapause DE transcripts | 1.2 MB |
| **S7** | Female diapause vs non-diapause DE transcripts | 1.6 MB |
| **S8** | Sex×diapause interaction DE transcripts | 535 KB |
| **S9** | Shared DE (no interaction) | 398 KB |
| **S10** | Shared DE (with interaction) | 42 KB |
| **S11** | Male-unique DE transcripts | 428 KB |
| **S12** | Female-unique DE transcripts | 701 KB |
| **S13** | WGCNA module overview | 6 KB |
| **S14** | Module DE enrichment results | 62 KB |
| **S15** | DE overlap GO enrichment | 61 KB |
| **S16** | DE overlap KEGG enrichment | 41 KB |
| **S17** | All modules GO enrichment | 701 KB |
| **S18** | All modules KEGG enrichment | 129 KB |

## Citation

If you use this code or data, please cite:

```
Torson AS, Melicher D, Yocum GD, Rinehart JP (2025). 
Diapause in the Colorado potato beetle is characterized by sex-specific gene expression.
Journal of Insect Physiology. DOI: XX.XXXX/XXXXX
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

**Corresponding Author**: Alex S. Torson  
Email: Alex.Torson@usda.gov  
USDA-ARS Edward T. Schafer Agricultural Research Center  
Biosciences Research Laboratory, Fargo, ND 58102, USA

## Acknowledgments

This research used resources provided by the SCINet project of the USDA Agricultural Research Service (project number 0500-00093-001-00-D). Funded by USDA-ARS project 3060-21220-03200D.

---

**Repository maintained by**: Alex Torson  
**Last updated**: December 2025
