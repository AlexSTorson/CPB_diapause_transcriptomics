################################################################################
#
# Title: Core WGCNA Analysis - Module Detection and Hub Gene Identification
#        Colorado Potato Beetle (CPB) Diapause Gene Co-expression Networks
#
# Author: Alex Torson
# Institution: USDA-ARS
# Email: Alex.Torson@usda.gov
#
# Description: Core WGCNA analysis for module detection, network 
#              construction, and hub gene identification. 
#
# Dependencies: DESeq2 object from previous differential expression analysis
#
# Input files:
#   - ./DESeq2/cpb_diapause_deseq2_object
#   - ./DESeq2/cpb_sample_information.csv
#   - ./Annotation/cpb_diapause_annotation.csv
#
# Output files:
#   - 09_core_results/*.csv (analysis summaries and module assignments)
#   - 09_transcript_lists/*.txt (gene lists by module)
#   - 09_network_relationships/*.csv (module connectivity data)
#   - 09_intermediate_data/core_wgcna_results.rds (R objects for downstream)
#
################################################################################

# Package Loading and Setup -----------------------------------------------

library(WGCNA)
cor <- WGCNA::cor
enableWGCNAThreads(nThreads = 8)

library(DESeq2)
library(tidyverse)
library(pheatmap)
library(corrplot)
library(RColorBrewer)

allowWGCNAThreads()
set.seed(12345)

# Working Directory and Output Setup --------------------------------------

setwd("/project/igb_fargo/cpb_diapause_rnaseq/wgcna/")

output_dirs <- c("09_core_results", "09_transcript_lists", "09_network_relationships", "09_intermediate_data")
walk(output_dirs, ~ if (!dir.exists(.x)) dir.create(.x, recursive = TRUE))

# Data Preparation ---------------------------------------------------------

# Load DESeq object
load("./01_DESeq2/cpb_diapause_deseq2_object")

# Extract normalized counts and sample information
norm_counts <- counts(dds, normalized = TRUE)
sample_df <- read_csv("./01_DESeq2/cpb_sample_information.csv") %>%
  dplyr::select(Sample_ID, status, sex)

# Create ordered sample info matching norm_counts columns
sample_order <- match(colnames(norm_counts), sample_df$Sample_ID)
sample_info <- sample_df[sample_order, ] %>%
  mutate(condition = str_c(sex, status, sep = "_"))

# Load gene annotations
raw_annotations <- read_csv("./00_Annotation/cpb_diapause_annotation.csv")
gene_annotations <- raw_annotations %>%
  dplyr::select(any_of(c("qry_transcript_id", "ref_transcript_id", "ref_transcript_name", "gene_annotation"))) %>%
  distinct(qry_transcript_id, .keep_all = TRUE)

# Data Filtering and Preprocessing ----------------------------------------

# Step 1: Filter low-expression genes (>10 counts in ≥3 samples)
keep_genes <- rowSums(norm_counts >= 10) >= 3
filtered_counts <- norm_counts[keep_genes, ]

# Step 2: Additional variance filtering
variance_filter <- apply(filtered_counts, 1, var) > 0
expr_data_filtered <- filtered_counts[variance_filter, ]

# Log2 transform with pseudocount
expr_data <- log2(expr_data_filtered + 1)

# Transpose for WGCNA (samples as rows, genes as columns)
expr_data <- t(expr_data)

# Create gene filtering summary
gene_filter_summary <- tibble(
  Step = c("Initial", "Expression_Filter", "Variance_Filter", "Final"),
  N_Genes = c(nrow(norm_counts), sum(keep_genes), sum(variance_filter), ncol(expr_data)),
  Percentage_Retained = round(N_Genes / nrow(norm_counts) * 100, 1)
)

# Network Construction and Module Detection -------------------------------

# Choose soft-thresholding power
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
sft <- pickSoftThreshold(expr_data, powerVector = powers, verbose = 5)
soft_power_all <- sft$powerEstimate

if (is.na(soft_power_all)) {
  soft_power_all <- 6
}

# Construct network and detect modules
net_all <- blockwiseModules(expr_data, 
                            power = soft_power_all,
                            TOMType = "unsigned", 
                            minModuleSize = 30,
                            reassignThreshold = 0, 
                            mergeCutHeight = 0.25,
                            numericLabels = TRUE, 
                            pamRespectsDendro = FALSE,
                            saveTOMs = FALSE,
                            verbose = 3)

# Convert numeric labels to colors for visualization
module_colors <- labels2colors(net_all$colors)

# Calculate module eigengenes (first principal component of each module)
MEs <- moduleEigengenes(expr_data, module_colors)$eigengenes

# Clean Module Numbering System -------------------------------------------

# Create systematic module numbering: Module_00 for grey, Module_01+ by size
module_sizes <- table(module_colors)
module_sizes_sorted <- sort(module_sizes, decreasing = TRUE)

# Create color-to-number mapping
color_to_number_map_clean <- setNames(
  sprintf("%02d", seq(0, length(module_sizes_sorted) - 1)), 
  names(module_sizes_sorted)
)

# Ensure grey module gets 00 regardless of size
if ("grey" %in% names(color_to_number_map_clean)) {
  color_to_number_map_clean["grey"] <- "00"
  other_modules <- names(module_sizes_sorted)[names(module_sizes_sorted) != "grey"]
  color_to_number_map_clean[other_modules] <- sprintf("%02d", seq(1, length(other_modules)))
}

# Hub Gene Analysis -------------------------------------------------------

# Calculate module membership (kME) for all genes
module_membership <- cor(expr_data, MEs, use = "pairwise.complete.obs")
rownames(module_membership) <- colnames(expr_data)
colnames(module_membership) <- names(MEs)

# Calculate intramodular connectivity for all genes
connectivity <- intramodularConnectivity.fromExpr(
  datExpr = expr_data,
  colors = module_colors,
  power = soft_power_all,
  scaleByMax = FALSE
)

# Identify hub genes using both kWithin and kME
hub_genes <- connectivity %>%
  as_tibble(rownames = "gene_index") %>%
  mutate(
    gene_index = as.numeric(gene_index),
    gene = colnames(expr_data)[gene_index],
    module_color = module_colors[gene_index]
  ) %>%
  # Add clean module numbering
  mutate(
    module_number = color_to_number_map_clean[module_color],
    module_id = paste0("Module_", module_number)
  ) %>%
  # Exclude grey module genes
  filter(module_color != "grey") %>%
  # Add kME values
  mutate(
    kME = map2_dbl(gene, module_color, 
                   ~module_membership[.x, paste0("ME", .y)]),
    kWithin = pmax(0, kWithin, na.rm = TRUE)
  ) %>%
  # Calculate hub metrics within each module
  group_by(module_color) %>%
  mutate(
    kWithin_percentile = percent_rank(kWithin),
    kME_percentile = percent_rank(abs(kME)),
    combined_percentile = (kWithin_percentile + kME_percentile) / 2
  ) %>%
  ungroup() %>%
  # Define hub categories
  mutate(
    hub_category = case_when(
      kWithin_percentile >= 0.95 & kME_percentile >= 0.95 ~ "key_regulator",
      kWithin_percentile >= 0.90 & kME_percentile >= 0.90 ~ "major_hub",
      kWithin_percentile >= 0.80 & kME_percentile >= 0.80 ~ "hub_gene",
      TRUE ~ "peripheral"
    )
  ) %>%
  arrange(desc(combined_percentile))

# Module Connectivity Analysis --------------------------------------------

# Calculate module-module correlations
ME_cor <- cor(MEs, use = "pairwise.complete.obs")

# Analyze module connectivity patterns
module_connectivity <- abs(ME_cor)
diag(module_connectivity) <- 0

module_connectivity_summary <- module_connectivity %>%
  as_tibble(rownames = "module") %>%
  pivot_longer(-module, names_to = "connected_module", values_to = "correlation") %>%
  group_by(module) %>%
  summarise(
    mean_connectivity = mean(correlation, na.rm = TRUE),
    max_connectivity = max(correlation, na.rm = TRUE),
    strong_connections = sum(correlation > 0.50, na.rm = TRUE),
    total_connectivity = sum(correlation, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    connectivity_rank = rank(desc(total_connectivity)),
    hub_module = strong_connections > 20 & mean_connectivity > 0.30
  ) %>%
  arrange(desc(total_connectivity))

# Comprehensive Module Summary --------------------------------------------

# Create comprehensive summary with annotations
gene_to_module_mapping <- tibble(
  qry_transcript_id = colnames(expr_data),
  module_color = module_colors,
  module_number = color_to_number_map_clean[module_colors],
  module_id = paste0("Module_", module_number)
) %>%
  mutate(
    module_type = case_when(
      module_color == "grey" ~ "unassigned",
      TRUE ~ "assigned"
    )
  )

# Merge with gene annotations
comprehensive_module_summary <- gene_to_module_mapping %>%
  left_join(gene_annotations, by = "qry_transcript_id") %>%
  group_by(module_color, module_id) %>%
  mutate(genes_in_module = n()) %>%
  ungroup() %>%
  dplyr::select(
    Module_ID = module_id,
    Module_Number = module_number, 
    Module_Color = module_color,
    Module_Type = module_type,
    Transcript_ID = qry_transcript_id,
    Reference_Transcript_ID = ref_transcript_id,
    Reference_Transcript_Name = ref_transcript_name,
    Gene_Annotation = gene_annotation,
    Genes_in_Module = genes_in_module
  ) %>%
  arrange(as.numeric(Module_Number), Transcript_ID)

# Create module overview table
module_overview_table <- comprehensive_module_summary %>%
  group_by(Module_ID, Module_Number, Module_Color, Module_Type) %>%
  summarise(
    Total_Genes = n(),
    Annotated_Genes = sum(!is.na(Gene_Annotation)),
    Annotation_Rate = round(Annotated_Genes / Total_Genes * 100, 1),
    .groups = 'drop'
  ) %>%
  arrange(as.numeric(Module_Number))

# Export Gene Lists by Module ---------------------------------------------

# Create gene-module annotation for exports
gene_module_annotations <- comprehensive_module_summary %>%
  rename(
    qry_transcript_id = Transcript_ID,
    module_id = Module_ID,
    module_color = Module_Color,
    ref_transcript_id = Reference_Transcript_ID,
    ref_transcript_name = Reference_Transcript_Name,
    gene_annotation = Gene_Annotation
  )

# Export gene lists for each colored module (exclude Module_00)
gene_module_annotations %>%
  filter(module_id != "Module_00") %>%
  group_by(module_id, module_color) %>%
  group_walk(~ {
    # Simple transcript list
    transcript_names <- .x %>%
      mutate(export_name = coalesce(ref_transcript_name, qry_transcript_id)) %>%
      pull(export_name)
    
    write_lines(transcript_names, 
                file.path("09_transcript_lists", str_c(.y$module_id, "_transcripts.txt")))
    
    # Annotated transcript list
    .x %>%
      dplyr::select(qry_transcript_id, ref_transcript_id, ref_transcript_name, gene_annotation) %>%
      write_csv(file.path("09_transcript_lists", str_c(.y$module_id, "_transcripts_annotated.csv")))
  })

# Export Module_00 (grey) separately if it exists
module_00_data <- gene_module_annotations %>% filter(module_id == "Module_00")

if (nrow(module_00_data) > 0) {
  module_00_transcript_names <- module_00_data %>%
    mutate(export_name = coalesce(ref_transcript_name, qry_transcript_id)) %>%
    pull(export_name)
  
  write_lines(module_00_transcript_names, 
              file.path("09_transcript_lists", "Module_00_transcripts.txt"))
  
  module_00_data %>%
    dplyr::select(qry_transcript_id, ref_transcript_id, ref_transcript_name, gene_annotation) %>%
    write_csv(file.path("09_transcript_lists", "Module_00_transcripts_annotated.csv"))
}

# Save Results -------------------------------------------------------------

# Main results tables
write_csv(comprehensive_module_summary, file.path("09_core_results", "Comprehensive_Module_Summary.csv"))
write_csv(module_overview_table, file.path("09_core_results", "Module_Overview_Table.csv"))
write_csv(gene_module_annotations, file.path("09_core_results", "transcript_module_assignments.csv"))
write_csv(hub_genes, file.path("09_core_results", "hub_genes_analysis.csv"))
write_csv(module_connectivity_summary, file.path("09_core_results", "module_connectivity_analysis.csv"))
write_csv(gene_filter_summary, file.path("09_core_results", "Gene_Filtering_Summary.csv"))

# Network summary statistics
network_summary <- tibble(
  Analysis = "Overall",
  N_Modules = length(unique(module_colors)) - 1,
  N_Genes = ncol(expr_data),
  Soft_Power = soft_power_all,
  Module_00_Genes = sum(module_colors == "grey"),
  Module_00_Percentage = round(sum(module_colors == "grey") / length(module_colors) * 100, 1),
  Hub_Genes = nrow(hub_genes),
  Hub_Gene_Percentage = round(nrow(hub_genes) / sum(module_colors != "grey") * 100, 1),
  Key_Regulators = sum(hub_genes$hub_category == "key_regulator"),
  Major_Hubs = sum(hub_genes$hub_category == "major_hub"),
  Hub_Modules = sum(module_connectivity_summary$hub_module)
)

write_csv(network_summary, file.path("09_core_results", "Network_Summary_Statistics.csv"))

# Save all objects for downstream analysis
saveRDS(list(
  expr_data = expr_data,
  sample_info = sample_info,
  gene_annotations = gene_annotations,
  module_colors = module_colors,
  module_eigengenes = MEs,
  module_membership = module_membership,
  hub_genes = hub_genes,
  color_to_number_map = color_to_number_map_clean,
  comprehensive_summary = comprehensive_module_summary,
  module_overview = module_overview_table,
  module_connectivity = module_connectivity_summary,
  ME_correlations = ME_cor,
  network_object = net_all,
  soft_power = soft_power_all,
  gene_filter_summary = gene_filter_summary
), file.path("09_intermediate_data", "core_wgcna_results.rds"))

# Session Info ------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "./09_core_results/09_session_info.txt")
