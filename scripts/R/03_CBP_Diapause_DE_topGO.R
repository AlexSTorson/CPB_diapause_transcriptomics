################################################################################
#
# Title: Gene Ontology Enrichment Analysis of Differential Expression
#        Colorado Potato Beetle (CPB) Diapause RNA-seq
#
# Author: Alex Torson
# Institution: USDA-ARS
# Email: Alex.Torson@usda.gov
#
# Description: Performs Gene Ontology enrichment analysis on differentially
#              expressed gene sets from various comparisons using topGO. Tests
#              for enrichment across Biological Process, Molecular Function,
#              and Cellular Component ontologies for each DE gene set.
#
# Dependencies: Differential expression results from DESeq2 
#
# Input files:
#   - 00_Annotation/cpb_diapause_annotation.csv
#   - 00_Annotation/cpb_diapause_go_omicsbox_export_one_seq_per_row.txt
#   - 01_DESeq2/*_DE_Transcripts.csv (differential expression results)
#
# Output files:
#   - 03_GO_Analyses/*_BP_Enrichment.csv (Biological Process enrichment)
#   - 03_GO_Analyses/*_MF_Enrichment.csv (Molecular Function enrichment)
#   - 03_GO_Analyses/*_CC_Enrichment.csv (Cellular Component enrichment)
#
################################################################################

# Package Loading ---------------------------------------------------------

library(tidyverse)
library(topGO)

# Working Directory and Output Setup --------------------------------------

setwd('/Users/alextorson/Library/CloudStorage/OneDrive-USDA/Torson_Lab/CPB_Diapause_RNA-seq/')

if (!dir.exists("./03_GO_Analyses/")) {
  dir.create("./03_GO_Analyses/", recursive = TRUE)
}

# Gene Ontology Data Preparation ------------------------------------------

# Load transcript-level GO annotations from OmicsBox
transID2GO <- read.delim(
  file = "./00_Annotation/cpb_diapause_go_omicsbox_export_one_seq_per_row.txt",
  sep = "\t",
  header = TRUE,
  fill = TRUE,
  na.strings = "",
  stringsAsFactors = FALSE
) %>%
  dplyr::rename(qry_transcript_id = Sequence.Name,
                go_id = Annotation.GO.ID) %>%
  filter(!is.na(go_id), go_id != "")

# Create gene-to-GO mapping list (required format for topGO)
transID2GOList <- strsplit(as.character(transID2GO$go_id), split = ",")
names(transID2GOList) <- transID2GO$qry_transcript_id

# Define gene universe (all genes with GO annotations)
geneUniverse <- as.character(transID2GO$qry_transcript_id)

# Load Differential Expression Results ------------------------------------

# Diapause effects in females (higher abundance in diapause)
d_vs_nd_female_up <- read.csv("./01_DESeq2/Diapause_vs_NonDiapause_Females_DE_Transcripts.csv") %>%
  filter(log2FoldChange > 0) %>%
  pull(qry_transcript_id)

# Diapause effects in females (higher abundance in non-diapause)
d_vs_nd_female_down <- read.csv("./01_DESeq2/Diapause_vs_NonDiapause_Females_DE_Transcripts.csv") %>%
  filter(log2FoldChange < 0) %>%
  pull(qry_transcript_id)

# Diapause effects in males (higher abundance in diapause)
d_vs_nd_male_up <- read.csv("./01_DESeq2/Diapause_vs_NonDiapause_Males_DE_Transcripts.csv") %>%
  filter(log2FoldChange > 0) %>%
  pull(qry_transcript_id)

# Diapause effects in males (higher abundance in non-diapause)
d_vs_nd_male_down <- read.csv("./01_DESeq2/Diapause_vs_NonDiapause_Males_DE_Transcripts.csv") %>%
  filter(log2FoldChange < 0) %>%
  pull(qry_transcript_id)

# Sex-specific diapause responses (interaction effects)
int_up <- read.csv("./01_DESeq2/Int_DE_Transcripts.csv") %>%
  filter(log2FoldChange > 0) %>%
  pull(qry_transcript_id)

int_down <- read.csv("./01_DESeq2/Int_DE_Transcripts.csv") %>%
  filter(log2FoldChange < 0) %>%
  pull(qry_transcript_id)

# Create named list of gene sets for enrichment analysis (excluding sex-controlled for now)
de_gene_sets <- list(
  d_vs_nd_female_up = d_vs_nd_female_up,
  d_vs_nd_female_down = d_vs_nd_female_down,
  d_vs_nd_male_up = d_vs_nd_male_up,
  d_vs_nd_male_down = d_vs_nd_male_down,
  int_up = int_up,
  int_down = int_down
)

# GO Enrichment Function --------------------------------------------------

perform_go_enrichment <- function(genes_of_interest, gene_set_name) {
  # Create binary gene list for topGO (1 = interesting, 0 = background)
  gene_list <- factor(as.integer(geneUniverse %in% genes_of_interest))
  names(gene_list) <- geneUniverse
  
  # Biological Process Enrichment -----------------------------------------
  
  # Create topGO object for Biological Process
  go_data_bp <- new("topGOdata",
                    description = "BP",
                    ontology = "BP",
                    allGenes = gene_list,
                    annot = annFUN.gene2GO,
                    gene2GO = transID2GOList,
                    nodeSize = 1)
  
  # Run Fisher's exact test
  result_fisher_bp <- runTest(go_data_bp, algorithm = "parentchild", statistic = "fisher")
  
  # Extract results
  all_go_bp <- usedGO(go_data_bp)
  results_bp <- GenTable(go_data_bp,
                         p_value = result_fisher_bp,
                         topNodes = length(all_go_bp),
                         numChar = 1000)
  
  # Add gene IDs for each enriched term
  results_bp$gene_id <- sapply(results_bp$GO.ID, function(x) {
    genes <- genesInTerm(go_data_bp, x)
    paste(genes[[1]][genes[[1]] %in% genes_of_interest], collapse = ", ")
  })
  
  # Filter results
  results_bp <- results_bp %>%
    filter(as.numeric(p_value) <= 0.01)
  
  # Save Biological Process results
  write.csv(results_bp,
            file = paste0(gene_set_name, "_BP_Enrichment.csv"),
            row.names = FALSE)
  
  # Molecular Function Enrichment -----------------------------------------
  
  # Create topGO object for Molecular Function
  go_data_mf <- new("topGOdata",
                    description = "MF",
                    ontology = "MF",
                    allGenes = gene_list,
                    annot = annFUN.gene2GO,
                    gene2GO = transID2GOList,
                    nodeSize = 1)
  
  # Run Fisher's exact test
  result_fisher_mf <- runTest(go_data_mf, algorithm = "parentchild", statistic = "fisher")
  
  # Extract results
  all_go_mf <- usedGO(go_data_mf)
  results_mf <- GenTable(go_data_mf,
                         p_value = result_fisher_mf,
                         topNodes = length(all_go_mf),
                         numChar = 1000)
  
  # Add gene IDs for each enriched term
  results_mf$gene_id <- sapply(results_mf$GO.ID, function(x) {
    genes <- genesInTerm(go_data_mf, x)
    paste(genes[[1]][genes[[1]] %in% genes_of_interest], collapse = ", ")
  })
  
  # Filter results
  results_mf <- results_mf %>%
    filter(as.numeric(p_value) <= 0.01)
  
  # Save Molecular Function results
  write.csv(results_mf,
            file = paste0(gene_set_name, "_MF_Enrichment.csv"),
            row.names = FALSE)
  
  # Cellular Component Enrichment -----------------------------------------
  
  # Create topGO object for Cellular Component
  go_data_cc <- new("topGOdata",
                    description = "CC",
                    ontology = "CC",
                    allGenes = gene_list,
                    annot = annFUN.gene2GO,
                    gene2GO = transID2GOList,
                    nodeSize = 1)
  
  # Run Fisher's exact test
  result_fisher_cc <- runTest(go_data_cc, algorithm = "parentchild", statistic = "fisher")
  
  # Extract results
  all_go_cc <- usedGO(go_data_cc)
  results_cc <- GenTable(go_data_cc,
                         p_value = result_fisher_cc,
                         topNodes = length(all_go_cc),
                         numChar = 1000)
  
  # Add gene IDs for each enriched term
  results_cc$gene_id <- sapply(results_cc$GO.ID, function(x) {
    genes <- genesInTerm(go_data_cc, x)
    paste(genes[[1]][genes[[1]] %in% genes_of_interest], collapse = ", ")
  })
  
  # Filter results
  results_cc <- results_cc %>%
    filter(as.numeric(p_value) <= 0.01)
  
  # Save Cellular Component results
  write.csv(results_cc,
            file = paste0(gene_set_name, "_CC_Enrichment.csv"),
            row.names = FALSE)
  
  return(list(BP = nrow(results_bp), MF = nrow(results_mf), CC = nrow(results_cc)))
}

# Run GO Enrichment Analysis ----------------------------------------------

# Change to GO_Analyses directory for output
setwd("./03_GO_Analyses/")

# Apply GO enrichment to all gene sets
enrichment_results <- map2(de_gene_sets, names(de_gene_sets), perform_go_enrichment)
names(enrichment_results) <- names(de_gene_sets)

# Session Info ------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "03_session_info.txt")
