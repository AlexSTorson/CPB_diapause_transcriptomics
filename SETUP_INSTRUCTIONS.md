# GitHub Repository Setup Instructions

This document provides step-by-step instructions for completing your GitHub repository setup and archiving it on Zenodo.

## Current Status ✓

The following files are already prepared:
- ✓ README.md
- ✓ LICENSE (MIT)
- ✓ .gitignore
- ✓ SCRIPTS_INDEX.md  
- ✓ Supplemental tables (all 18 tables in `data/supplemental_tables/`)
- ✓ Directory structure created
- ✓ Documentation READMEs for script directories

## What You Need to Do

### Step 1: Add Your Scripts to the Repository

Your bash and R scripts are in the Claude project but need to be copied into the repository directories.

**For bash scripts:**
The following 17 bash scripts need to be placed in `scripts/bash/`:
- 00_pipeline_config.sh
- 01_create_log_dirs.sh
- 02_data_download.sh
- 03_genome_references_download.sh
- 04_fastqc_pre_trimming.sh
- 05_fastp_quality_trim.sh
- 06_fastqc_post_trimming.sh
- 07_exon_and_splice_sites.sh
- 08_genome_index.sh
- 09_read_mapping.sh
- 10_sam_to_bam.sh
- 11_stringtie.sh
- 12_create_mergelist.sh
- 13_stringtie_merge.sh
- 13_taco_merge.sh
- 14_gffcompare.sh
- 15_count_estimation.sh
- 16_prepDE.sh
- 17_multiqc.sh
- run-full-pipeline.sh
- refresh_config.sh
- rel_config_paths.sh

Plus these Python helper scripts:
- hisat2_extract_exons.py
- hisat2_extract_splice_sites.py
- prepDE.py

**For R scripts:**
The following 14 R scripts need to be placed in `scripts/R/`:
- 00_CBP_Diapause_Annotation.R
- 01_CBP_Diapause_DESeq2.R
- 02_CPB_Diapause_Multivariate_Analysis.R
- 03_CBP_Diapause_DE_topGO.R
- 04_CBP_Diapause_DE_topGO_Summary_Plots.R
- 05_CPB_Diapause_DE_KEGG_Enrichment.R
- 06_CPB_Diapause_DE_KEGG_Summary.R
- 07_CBP_Diapause_DE_Overlap_Venn.R
- 08_CBP_Diapause_DE_Overlap_Enrichment_Summary_Plots.R
- 09_CPB_Diapause_WGCNA.R
- 10_CPB_Diapause_WGCNA_Dimorphism_Analysis.R
- 11_CPB_Diapause_WGCNA_Figures.R
- 12_CPB_Diapause_WGCNA_GO_Enrichment_Analysis.R
- 13_CPB_Diapause_WGCNA_KEGG_Enrichment_Analysis.R
- 14_CBP_Diapause_Couchpotato.R

### Step 2: Initialize Git Repository

```bash
cd /home/claude/CPB_diapause_transcriptomics

# Initialize git repository
git init

# Add all files
git add .

# Create initial commit
git commit -m "Initial commit: RNA-seq analysis pipeline and supplemental data"
```

### Step 3: Create GitHub Repository

1. Go to https://github.com
2. Click "New repository"
3. Name it: `CPB_diapause_transcriptomics`
4. Description: "Sex-specific transcriptomics of diapause in Colorado potato beetle"
5. Keep it Public (required for Zenodo)
6. Do NOT initialize with README (you already have one)
7. Click "Create repository"

### Step 4: Push to GitHub

GitHub will show you commands. They'll look like:

```bash
git remote add origin https://github.com/YOUR-USERNAME/CPB_diapause_transcriptomics.git
git branch -M main
git push -u origin main
```

### Step 5: Create a GitHub Release

Before archiving on Zenodo, create a release:

1. Go to your repository on GitHub
2. Click "Releases" (right sidebar)
3. Click "Create a new release"
4. Tag version: `v1.0.0`
5. Release title: `Initial Release - Publication Version`
6. Description: Add publication info and DOI once available
7. Click "Publish release"

### Step 6: Archive on Zenodo

1. Go to https://zenodo.org
2. Sign in (create account if needed)
3. Click your name → "GitHub" (connect GitHub account)
4. Find `CPB_diapause_transcriptomics` in the repository list
5. Toggle the switch to ON to enable archiving
6. Go back to GitHub and create a new release (or use existing v1.0.0)
7. Zenodo will automatically create an archive with a DOI

### Step 7: Update README with DOI

Once Zenodo creates your DOI:

1. Copy the DOI (will look like: `10.5281/zenodo.XXXXXXX`)
2. Update the README.md badge at the top:
   ```markdown
   [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
   ```
3. Commit and push the change:
   ```bash
   git add README.md
   git commit -m "Add Zenodo DOI"
   git push
   ```

### Step 8: Update Your Manuscript

In your manuscript's Data Availability section, add:

```
Bash and R scripts for all analyses are available at GitHub 
(https://github.com/YOUR-USERNAME/CPB_diapause_transcriptomics) and permanently 
archived at Zenodo (DOI: 10.5281/zenodo.XXXXXXX).
```

## Verification Checklist

Before considering the repository complete, verify:

- [ ] All 17 bash scripts are in `scripts/bash/`
- [ ] All 14 R scripts are in `scripts/R/`
- [ ] All 3 Python helper scripts are in `scripts/bash/`
- [ ] All 18 supplemental tables are in `data/supplemental_tables/`
- [ ] README.md is complete and formatted correctly
- [ ] LICENSE file is present
- [ ] Repository is pushed to GitHub
- [ ] Release v1.0.0 is created
- [ ] Zenodo archive is created and DOI obtained
- [ ] README updated with correct DOI
- [ ] Manuscript updated with repository and DOI links

## Troubleshooting

**Problem**: Git is not installed
**Solution**: `sudo apt-get install git` or download from https://git-scm.com

**Problem**: Can't push to GitHub (authentication)
**Solution**: Use Personal Access Token instead of password. Generate at: GitHub → Settings → Developer settings → Personal access tokens

**Problem**: Files too large for GitHub
**Solution**: The .gitignore already excludes large data files. Only scripts and supplemental tables will be uploaded.

**Problem**: Zenodo won't connect to GitHub
**Solution**: Make sure repository is Public and you've authorized Zenodo to access your GitHub account

## Need Help?

If you encounter issues:
1. Check GitHub's documentation: https://docs.github.com
2. Check Zenodo's GitHub guide: https://help.zenodo.org/#github
3. Contact: Alex.Torson@usda.gov

---

**Note**: Once you've completed all steps, this file can be deleted or moved out of the repository before creating the final release.
