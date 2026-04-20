# CPB Diapause Transcriptomics

RNA-seq analysis of sex-specific gene expression during diapause in Colorado potato beetle (*Leptinotarsa decemlineata*).

## Citation

Torson, A. S., Melicher, D., Yocum, G. D., & Rinehart, J. P. (2026). Diapause in the Colorado potato beetle is characterized by sex-specific gene expression. Journal of Insect Physiology, 104990. DOI: https://doi.org/10.1016/j.jinsphys.2026.104990

## Data

- **Raw sequencing data**: NCBI BioProject PRJNA553565
- **Reference genome**: *L. decemlineata* RefSeq GCF_000500325.1

## Scripts

- **Bash pipeline**: `scripts/bash/` - RNA-seq processing (SRA → count matrices)
- **R analysis**: `scripts/R/` - Differential expression and network analysis
- **Overview**: See `SCRIPTS_INDEX.md`

## Software Requirements

**Bash**: fastp, FastQC, HISAT2, StringTie, TACO, SAMtools  
**R**: DESeq2, WGCNA, topGO, clusterProfiler, tidyverse

## License

US Government Work (Public Domain)

## Contact

Alex Torson - Alex.Torson@usda.gov  
USDA-ARS, Fargo, ND
