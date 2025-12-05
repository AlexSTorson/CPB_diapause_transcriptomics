# Quick Reference Card

## Repository Location
```
/mnt/user-data/outputs/CPB_diapause_transcriptomics/
```

## What's Complete ✓
- Repository structure
- README.md with full documentation
- All 18 supplemental tables
- LICENSE file (MIT)
- .gitignore
- Complete documentation

## What You Need to Add ⚠
- 17 bash scripts → `scripts/bash/`
- 14 R scripts → `scripts/R/`
- 3 Python helpers → `scripts/bash/`

## GitHub Setup (5 Steps)

1. **Add scripts** to repository directories
2. **Initialize git**: `git init && git add . && git commit -m "Initial commit"`
3. **Create GitHub repo** at github.com (keep public!)
4. **Push code**: `git remote add origin URL && git push -u origin main`
5. **Create release**: v1.0.0 on GitHub

## Zenodo Archiving (3 Steps)

1. **Connect GitHub** at zenodo.org → GitHub settings
2. **Enable repository** toggle switch for your repo
3. **Get DOI** automatically created from release

## Update Manuscript

Replace this in your Data Availability section:

**Before**:
> "Bash and R scripts are available as supplemental material."

**After**:
> "Bash and R scripts are available at GitHub (https://github.com/USERNAME/CPB_diapause_transcriptomics) and permanently archived at Zenodo (DOI: 10.5281/zenodo.XXXXXXX)."

## Key Files to Review

- `REPOSITORY_SUMMARY.md` - This summary
- `SETUP_INSTRUCTIONS.md` - Complete step-by-step guide  
- `SCRIPTS_INDEX.md` - All 31 scripts listed
- `README.md` - Main repository page

## Scripts to Copy

### Bash (scripts/bash/)
00, 01, 02, 03, 04, 05, 06, 07, 08, 09, 10, 11, 12, 13 (two versions), 14, 15, 16, 17
Plus: run-full-pipeline.sh, refresh_config.sh, rel_config_paths.sh
Plus: hisat2_extract_*.py, prepDE.py

### R (scripts/R/)
00-14 (all R scripts from annotation through couchpotato analysis)

## Important Notes

- Keep repository PUBLIC (required for Zenodo)
- Create release BEFORE archiving on Zenodo
- Update README with DOI after Zenodo creates it
- Scripts referenced in manuscript as "File S1" and "File S2" are now in GitHub

## Timeline Estimate

- Adding scripts: 15 minutes
- GitHub setup: 10 minutes
- Zenodo setup: 5 minutes
- **Total: ~30 minutes**

---

**Need detailed help?** → Open `SETUP_INSTRUCTIONS.md`
