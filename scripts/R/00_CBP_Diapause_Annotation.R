################################################################################
#
# Title: Transcriptome Annotation Pipeline
#        Colorado Potato Beetle (CPB) Diapause RNA-seq
#
# Author: Alex Torson
# Institution: USDA-ARS
# Email: Alex.Torson@usda.gov
#
# Description: Annotation of Colorado potato beetle diapause transcriptome sequences
#              combining RefSeq annotations with novel transcript annotations from
#              BLAST analysis conduced in Omicsbox. Creates full annotation data set for
#              downstream analyses.
#
# Dependencies: Output files from transcriptome assembly pipeline
#
# Input files:
#   - tacocomp.assembly.gtf.tmap (gffcompare tracking file)
#   - GCF_000500325.1_Ldec_2.0_feature_table.txt (RefSeq feature table)
#   - CPB_Diapause_blast2go_topblast.txt (OmicsBox BLAST results)
#
# Output files:
#   - 00_annotation/cpb_diapause_annotation.csv (comprehensive annotation file)
#
################################################################################

# Package Loading ---------------------------------------------------------

library(tidyverse)
library(seqinr)

# Working Directory and Output Setup --------------------------------------

setwd("/Users/alextorson/Library/CloudStorage/OneDrive-USDA/Torson_Lab/CPB_Diapause_RNA-seq/00_annotation/")

# Data Loading and Processing ----------------------------------------------

# Load gffcompare tracking file
taco_annot <- read.delim(file = "tacocomp.assembly.gtf.tmap", 
                         header = TRUE, sep = "\t", na.strings = c("-", ""))

# Clean and rename columns
taco_annot <- taco_annot %>%
  dplyr::rename(ref_transcript_id = ref_id, qry_transcript_id = qry_id) %>%
  dplyr::select(ref_gene_id, qry_gene_id, ref_transcript_id, qry_transcript_id, class_code) 

# Load RefSeq feature table
refseq_annot <- read.delim(file = "GCF_000500325.1_Ldec_2.0_feature_table.txt", 
                           header = TRUE, sep = "\t", dec = ".")

# Filter for RNA features and select relevant columns
refseq_annot <- refseq_annot %>%
  filter(feature != "CDS" | feature != "gene") %>%
  dplyr::select(product_accession, name, symbol) %>%
  dplyr::rename(ref_transcript_id = product_accession, 
                ref_transcript_name = name,
                ref_gene_id = symbol)

# Merge RefSeq annotations with transcript mappings
diapause_annot <- taco_annot %>%
  left_join(refseq_annot, by = c("ref_transcript_id", "ref_gene_id"))

# Create gene name by removing variant information
diapause_annot$ref_gene_name <- sub(",.*", "", diapause_annot$ref_transcript_name)

# Create column that denotes whether the gene was prev. existing or IDed as novel by stringtie
diapause_annot$gene_annotation <- ifelse(is.na(diapause_annot$ref_gene_id), "novel", "refseq")

# Create list of transcript IDs that do not have a reference annotation
no_annot_trans <- diapause_annot %>%
  filter(is.na(ref_gene_id)) %>%
  dplyr::select(qry_transcript_id)

# Load BLAST2GO annotation result to add novel transcript annotations
blast_annot <- read.delim(file = "CPB_Diapause_blast2go_topblast.txt",
                          header = TRUE, sep = "\t", quote = "", na.strings = "")

# Process BLAST annotations and merge with unannotated transcripts
novel_omics_blast <- blast_annot %>%
  dplyr::rename(qry_transcript_id = Sequence.name,
                ref_transcript_name = Sequence.desc.,
                ref_transcript_id = Hit.ACC) %>%
  dplyr::select(qry_transcript_id, ref_transcript_name, ref_transcript_id) %>%
  right_join(no_annot_trans)

# Clean up extra accession numbers
novel_omics_blast$ref_transcript_id <- sub(",.*", "", as.character(novel_omics_blast$ref_transcript_id))

# Create gene name by removing isoform info
novel_omics_blast$ref_gene_name <- sub("isoform X[0-9]*", "", as.character(novel_omics_blast$ref_transcript_name))

# Add asterisk to novel gene names to mark them as novel
novel_omics_blast$ref_gene_name <- sub("^", "* ", as.character(novel_omics_blast$ref_gene_name))

# Merge novel annotations with existing annotations
final_annotation <- diapause_annot %>%
  left_join(novel_omics_blast, by = "qry_transcript_id") %>%
  mutate(ref_transcript_name = coalesce(ref_transcript_name.x, ref_transcript_name.y)) %>%
  mutate(ref_transcript_id = coalesce(ref_transcript_id.x, ref_transcript_id.y)) %>%
  mutate(ref_gene_name = coalesce(ref_gene_name.x, ref_gene_name.y)) %>%
  dplyr::select(ref_gene_id, qry_gene_id, ref_gene_name, gene_annotation, ref_transcript_id, qry_transcript_id, ref_transcript_name, class_code)

# Save final annotation file
write_csv(final_annotation, "cpb_diapause_annotation.csv")

# Session Info ------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "./00_Annotation/00_session_info.txt")
