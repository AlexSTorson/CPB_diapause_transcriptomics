################################################################################
#
# Title: KEGG Pathway Enrichment Analysis and Visualization
#        Colorado Potato Beetle (CPB) Diapause RNA-seq
#
# Author: Alex Torson
# Institution: USDA-ARS
# Email: Alex.Torson@usda.gov
#
# Description: Performs KEGG pathway enrichment analysis on differentially
#              expressed gene sets and creates pathway visualizations with
#              fold change data overlaid. Tests for enriched metabolic and 
#              signaling pathways and generates visual pathway maps.
#
# Dependencies: Differential expression results from DESeq2
#
# Input files:
#   - 00_Annotation/cpb_diapause_KAAS.txt (KEGG ortholog mappings)
#   - 00_Annotation/cpb_diapause_annotation.csv (gene annotations)
#   - 01_DESeq2/*_DE_Transcripts.csv (differential expression results)
#
# Output files:
#   - 05_KEGG_Enrichment/*_enrichment_results.csv (enrichment statistics)
#   - 05_KEGG_Enrichment/Mapping_Plots_*/*.pathview.png (pathway visualizations)
#
################################################################################

# Package Loading ---------------------------------------------------------

library(tidyverse)
library(pathview)
library(gage)
library(clusterProfiler)

# Working Directory and Output Setup --------------------------------------

setwd("/Users/alextorson/Library/CloudStorage/OneDrive-USDA/Torson_Lab/CPB_Diapause_RNA-seq/")

if (!dir.exists("./05_KEGG_Enrichment/")) {
  dir.create("./05_KEGG_Enrichment/", recursive = TRUE)
}

# Change to KEGG enrichment directory for output
setwd("./05_KEGG_Enrichment/")

# Create output directories
kegg_dirs <- c(
  "Mapping_Plots_KO_IDs",
  "Mapping_Plots_Enzymes", 
  "Lists"
)

walk(kegg_dirs, ~ if (!dir.exists(.x)) dir.create(.x, recursive = TRUE))

# KEGG Data Preparation ---------------------------------------------------

# Load KEGG ortholog mappings from KAAS annotation
cpb_kegg <- read.table(
  "../00_Annotation/cpb_diapause_KAAS.txt",
  header = FALSE,
  colClasses = c("character", "factor"),
  sep = "",
  fill = TRUE,
  col.names = c("qry_transcript_id", "kegg_id"),
  na.strings = ""
) %>%
  drop_na()

# Load gene annotations
cpb_annot <- read.csv('../00_Annotation/cpb_diapause_annotation.csv')

# Merge KEGG IDs with gene annotations
cpb_annot <- cpb_annot %>%
  left_join(cpb_kegg, by = "qry_transcript_id")

# Create gene-to-KEGG mapping
kegg_genes <- cpb_annot %>%
  dplyr::select(qry_transcript_id, kegg_id) %>%
  drop_na() %>%
  distinct()

# Get current KEGG pathway gene sets
httr::set_config(httr::config(ssl_verifypeer = FALSE))
kg.ko <- kegg.gsets("ko")
kegg.gs <- kg.ko$kg.sets[kg.ko$sigmet.idx]  # Remove disease pathways

# Extract pathway IDs for filtering
kegg.gs_names <- names(kegg.gs)
kegg.gs_names <- as.data.frame(gsub(" .*$", "", kegg.gs_names))
names(kegg.gs_names) <- "ID"

# Load Differential Expression Results ------------------------------------

# Diapause effects in females
d_vs_nd_female <- read.csv("../01_DESeq2/Diapause_vs_NonDiapause_Females_DE_Transcripts.csv") %>%
  dplyr::select(qry_transcript_id, log2FoldChange, padj)

# Diapause effects in males
d_vs_nd_male <- read.csv("../01_DESeq2/Diapause_vs_NonDiapause_Males_DE_Transcripts.csv") %>%
  dplyr::select(qry_transcript_id, log2FoldChange, padj)

# Sex-specific diapause responses (interaction effects)
int <- read.csv("../01_DESeq2/Int_DE_Transcripts.csv") %>%
  dplyr::select(qry_transcript_id, log2FoldChange, padj)

# Create list of all differential expression results
de_list <- list(
  d_vs_nd_female = d_vs_nd_female,
  d_vs_nd_male = d_vs_nd_male,
  int = int
)

# KEGG Enrichment and Visualization Function ------------------------------

perform_kegg_analysis <- function(de_df, label) {
  # Store the original working directory
  original_wd <- getwd()
  
  # Join differential expression data with KEGG annotations
  de_kegg <- de_df %>%
    left_join(kegg_genes, by = "qry_transcript_id") %>%
    drop_na() %>%
    mutate(log2FoldChange = as.numeric(log2FoldChange))
  
  # Create character vector of KEGG gene IDs
  de_kegg_chr <- as.character(de_kegg$kegg_id)
  
  # Create named vector of fold changes with KEGG IDs as names
  de_kegg_fc <- de_kegg %>%
    group_by(kegg_id) %>%
    summarize(log2FoldChange = mean(log2FoldChange, na.rm = TRUE), .groups = 'drop') %>%
    {setNames(.$log2FoldChange, .$kegg_id)}
  
  # Create named vector of transcript IDs with KEGG IDs as names
  de_kegg_transcripts <- with(de_kegg, setNames(qry_transcript_id, kegg_id))
  
  # Match genes to pathways
  de_pathway_genes <- lapply(kegg.gs, function(x)
    as.vector(na.omit(de_kegg_transcripts[x])))
  
  # Remove pathway IDs without differentially expressed genes
  de_pathway_genes <- de_pathway_genes[lengths(de_pathway_genes) > 0]
  
  # Create character vector of pathway IDs
  de_pathway_names_full <- as.character(names(de_pathway_genes))
  
  # Remove problematic pathways from plotting
  excluded_pathways <- c(
    "ko01100 Metabolic pathways",
    "ko01110 Biosynthesis of secondary metabolites",
    "ko04723 Retrograde endocannabinoid signaling",
    "ko01120 Microbial metabolism in diverse environments",
    "ko04215 Apoptosis - multiple species",
    "ko00510 N-Glycan biosynthesis",
    "ko00511 Other glycan degradation",
    "ko00512 Mucin type O-glycan biosynthesis",
    "ko00513 Various types of N-glycan biosynthesis",
    "ko00514 Other types of O-glycan biosynthesis",
    "ko00515 Mannose type O-glycan biosynthesis"
  )
  
  de_pathway_names_full <- de_pathway_names_full[!de_pathway_names_full %in% excluded_pathways]
  
  # Extract pathway IDs and names
  de_pathway_ids <- gsub(" .*$", "", de_pathway_names_full)
  de_pathway_names <- gsub("^.*? ", "", de_pathway_names_full)
  
  # Calculate fold change limit for visualization
  all_fold_changes <- de_kegg$log2FoldChange
  max_fc <- quantile(abs(all_fold_changes), 0.95, na.rm = TRUE)
  fc_limit <- ceiling(max_fc)
  
  # Create absolute paths for output directories
  plot_dir_ko <- file.path(original_wd, "Mapping_Plots_KO_IDs", label)
  plot_dir_enzymes <- file.path(original_wd, "Mapping_Plots_Enzymes", label)
  list_dir <- file.path(original_wd, "Lists", label)
  
  walk(c(plot_dir_ko, plot_dir_enzymes, list_dir), ~ {
    if (!dir.exists(.x)) dir.create(.x, recursive = TRUE)
  })
  
  # Generate KO ID pathway visualizations
  setwd(plot_dir_ko)
  
  ko_results <- map(seq_along(de_pathway_ids), function(i) {
    tryCatch({
      path_id <- de_pathway_ids[i]
      path_name <- de_pathway_names[i]
      path_name_clean <- gsub("[^a-zA-Z0-9_]", "_", path_name)
      
      out_suffix <- paste0("KO_", label, "_", path_name_clean)
      
      pathview_result <- pathview(
        gene.data = de_kegg_fc,
        pathway.id = path_id,
        species = "ko",
        out.suffix = out_suffix,
        kegg.native = TRUE,
        multi.state = FALSE,
        same.layer = FALSE,  # KO IDs use different layer
        limit = list(gene = c(-fc_limit, fc_limit)),
        low = "red",
        mid = "grey90",
        high = "blue",
        plot.col.key = TRUE,
        key.pos = "topright",
        bins = list(gene = 20)
      )
      
      return(out_suffix)
    }, error = function(e) {
      return(NULL)
    })
  })
  
  # Reset to KO directory and generate enzyme visualizations
  setwd(plot_dir_enzymes)
  
  enzyme_results <- map(seq_along(de_pathway_ids), function(i) {
    tryCatch({
      path_id <- de_pathway_ids[i]
      path_name <- de_pathway_names[i]
      path_name_clean <- gsub("[^a-zA-Z0-9_]", "_", path_name)
      
      out_suffix <- paste0("Enzyme_", label, "_", path_name_clean)
      
      pathview_result <- pathview(
        gene.data = de_kegg_fc,
        pathway.id = path_id,
        species = "ko",
        out.suffix = out_suffix,
        kegg.native = TRUE,
        multi.state = FALSE,
        same.layer = TRUE,  # Enzyme names use same layer
        limit = list(gene = c(-fc_limit, fc_limit)),
        low = "red",
        mid = "grey90",
        high = "blue",
        plot.col.key = TRUE,
        key.pos = "topright",
        bins = list(gene = 20)
      )
      
      return(out_suffix)
    }, error = function(e) {
      return(NULL)
    })
  })
  
  # Perform enrichment analysis
  enrich <- enrichKEGG(
    gene = de_kegg_chr,
    organism = "ko",
    keyType = 'kegg',
    pvalueCutoff = 0.01
  )
  
  # Extract and filter enrichment results
  enrich_df <- data.frame(enrich@result) %>%
    filter(qvalue <= 0.01) %>%
    subset(ID %in% kegg.gs_names$ID) %>%  # Remove disease pathways
    unite(pathway, ID, Description, remove = FALSE, sep = " ")
  
  # Save enrichment results
  setwd(list_dir)
  write.csv(enrich_df,
            file = paste0(label, "_enrichment_results.csv"),
            row.names = FALSE)
  
  # Create pathway-gene mapping file for easy reference
  if (nrow(enrich_df) > 0) {
    pathway_gene_mapping <- map_dfr(1:nrow(enrich_df), function(i) {
      pathway_id <- enrich_df$ID[i]
      pathway_name <- enrich_df$Description[i]
      
      # Get genes in this pathway
      pathway_genes <- strsplit(enrich_df$geneID[i], "/")[[1]]
      
      # Create mapping dataframe
      data.frame(
        pathway_id = pathway_id,
        pathway_name = pathway_name,
        kegg_id = pathway_genes,
        stringsAsFactors = FALSE
      )
    }) %>%
      # Add transcript IDs by joining with KEGG annotations
      left_join(kegg_genes, by = "kegg_id") %>%
      dplyr::select(pathway_id, pathway_name, qry_transcript_id, kegg_id) %>%
      filter(!is.na(qry_transcript_id)) %>%
      arrange(pathway_id, qry_transcript_id)
    
    # Save pathway-gene mapping
    write.csv(pathway_gene_mapping,
              file = paste0(label, "_pathway_gene_ids.csv"),
              row.names = FALSE)
  }
  
  # Always return to original working directory
  setwd(original_wd)
  
  return(list(
    n_pathways = length(de_pathway_ids),
    n_enriched = nrow(enrich_df),
    ko_plots = sum(!sapply(ko_results, is.null)),
    enzyme_plots = sum(!sapply(enzyme_results, is.null))
  ))
}

# Run KEGG Analysis on All Comparisons ------------------------------------

# Apply KEGG analysis to all differential expression datasets
kegg_results <- imap(de_list, perform_kegg_analysis)

# Session Info ------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "05_session_info.txt")
