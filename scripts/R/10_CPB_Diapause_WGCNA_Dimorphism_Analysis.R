################################################################################
#
# Title: Sexual Dimorphism Analysis for WGCNA Modules
#        Colorado Potato Beetle (CPB) Diapause Gene Co-expression Networks
#
# Author: Alex Torson
# Institution: USDA-ARS
# Email: Alex.Torson@usda.gov
#
# Description: Analysis of sexual dimorphism in module expression patterns
#              and DE enrichment testing for WGCNA modules. Focuses on 
#              manuscript analyses: sex differences, sex×diapause interactions,
#              and DE category enrichment.
#
# Dependencies: Must run 09_core_wgcna_analysis.R first
#
# Input files:
#   - 09_intermediate_data/core_wgcna_results.rds
#   - DE_Overlaps/Annotated_Overlap_Lists/*.csv
#
# Output files:
#   - 10_dimorphism_analysis/sex_difference_tests.csv
#   - 10_dimorphism_analysis/interaction_effect_tests.csv
#   - 10_dimorphism_analysis/de_module_enrichments.csv
#   - 10_dimorphism_analysis/module_connectivity_analysis.csv
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

dimorphism_dirs <- c("10_dimorphism_analysis", "10_module_profiles", "10_intermediate_data")
walk(dimorphism_dirs, ~ if (!dir.exists(.x)) dir.create(.x, recursive = TRUE))

# Load Core WGCNA Results --------------------------------------------------

# Load core results from script 09
core_results <- readRDS(file.path("09_intermediate_data", "core_wgcna_results.rds"))

# Extract key objects
expr_data <- core_results$expr_data
sample_info <- core_results$sample_info
gene_annotations <- core_results$gene_annotations
module_colors <- core_results$module_colors
MEs <- core_results$module_eigengenes
hub_genes <- core_results$hub_genes
color_to_number_map_clean <- core_results$color_to_number_map
comprehensive_summary <- core_results$comprehensive_summary

# Analysis 1: Module Expression Patterns and Sexual Dimorphism -----------

# Calculate module expression profiles across all samples
module_profiles <- MEs %>%
  as_tibble(rownames = "sample_id") %>%
  left_join(sample_info %>% mutate(sample_id = rownames(expr_data)), by = "sample_id") %>%
  pivot_longer(cols = starts_with("ME"), names_to = "module_eigen", values_to = "expression") %>%
  mutate(
    module_color = str_remove(module_eigen, "ME"),
    module_number = color_to_number_map_clean[module_color],
    module_id = paste0("Module_", module_number)
  ) %>%
  filter(module_color != "grey")

# Sex difference tests for each module within each diapause status
sex_difference_tests <- module_profiles %>%
  group_by(module_id, status) %>%
  filter(n_distinct(sex) == 2, n() >= 6) %>%
  summarise(
    n_samples = n(),
    male_mean = mean(expression[sex == "Male"], na.rm = TRUE),
    female_mean = mean(expression[sex == "Female"], na.rm = TRUE),
    t_test_result = list(t.test(expression ~ sex)),
    t_statistic = t_test_result[[1]]$statistic,
    p_value = t_test_result[[1]]$p.value,
    male_sd = sd(expression[sex == "Male"], na.rm = TRUE),
    female_sd = sd(expression[sex == "Female"], na.rm = TRUE),
    pooled_sd = sqrt(((sum(sex == "Male") - 1) * male_sd^2 + 
                        (sum(sex == "Female") - 1) * female_sd^2) / 
                       (n() - 2)),
    cohens_d = (male_mean - female_mean) / pooled_sd,
    .groups = 'drop'
  ) %>%
  dplyr::select(-t_test_result) %>%
  mutate(
    p_value_adj = p.adjust(p_value, method = "BH"),
    significant = p_value_adj < 0.05,
    effect_size_category = case_when(
      abs(cohens_d) >= 0.8 ~ "large",
      abs(cohens_d) >= 0.5 ~ "medium", 
      abs(cohens_d) >= 0.2 ~ "small",
      TRUE ~ "negligible"
    )
  )

# Interaction effect tests (sex × diapause status)
interaction_tests <- module_profiles %>%
  group_by(module_id) %>%
  filter(n_distinct(sex) == 2, n_distinct(status) == 2, n() >= 8) %>%
  summarise(
    n_samples = n(),
    anova_result = list(aov(expression ~ sex * status)),
    interaction_p = summary(anova_result[[1]])[[1]]["sex:status", "Pr(>F)"],
    total_ss = sum(summary(anova_result[[1]])[[1]][["Sum Sq"]]),
    interaction_ss = summary(anova_result[[1]])[[1]]["sex:status", "Sum Sq"],
    eta_squared = interaction_ss / total_ss,
    .groups = 'drop'
  ) %>%
  dplyr::select(-anova_result) %>%
  mutate(
    interaction_p_adj = p.adjust(interaction_p, method = "BH"),
    significant_interaction = interaction_p_adj < 0.05,
    effect_size_category = case_when(
      eta_squared >= 0.14 ~ "large",
      eta_squared >= 0.06 ~ "medium",
      eta_squared >= 0.01 ~ "small", 
      TRUE ~ "negligible"
    )
  )

# Analysis 2: DE-Module Enrichment (Hypergeometric tests) ----------------

# Test which modules are enriched for different categories of DE genes
de_overlap_dir <- file.path(getwd(), "07_DE_Overlaps", "Annotated_Overlap_Lists")

# Load DE gene lists if they exist
if (dir.exists(de_overlap_dir)) {
  
  de_files <- list.files(de_overlap_dir, pattern = "*.csv$", full.names = TRUE)
  
  if (length(de_files) > 0) {
    
    # Load each DE category
    de_categories <- map(de_files, ~ {
      file_name <- basename(.x)
      category_name <- str_remove(file_name, ".csv$")
      
      de_genes <- read_csv(.x, show_col_types = FALSE) %>%
        pull(qry_transcript_id)
      
      tibble(category = category_name, genes = list(de_genes))
    }) %>%
      bind_rows()
    
    # Test enrichment for each module and DE category
    de_enrichments <- map_dfr(1:nrow(de_categories), function(i) {
      
      category_name <- de_categories$category[i]
      de_genes <- de_categories$genes[[i]]
      
      # Test each module (exclude Module_00)
      comprehensive_summary %>%
        filter(Module_ID != "Module_00") %>%
        group_by(Module_ID) %>%
        summarise(
          module_size = n(),
          de_in_module = sum(Transcript_ID %in% de_genes),
          .groups = 'drop'
        ) %>%
        filter(module_size >= 10) %>%
        rowwise() %>%
        mutate(
          # Total genes in background
          total_genes = nrow(comprehensive_summary),
          total_de = length(de_genes),
          # Hypergeometric test
          hypergeom_p = phyper(
            q = de_in_module - 1,
            m = total_de,
            n = total_genes - total_de,
            k = module_size,
            lower.tail = FALSE
          ),
          # Enrichment ratio
          expected_overlap = (module_size * total_de) / total_genes,
          enrichment_ratio = de_in_module / expected_overlap,
          # Category information
          de_category = category_name
        ) %>%
        ungroup()
    })
    
    # Multiple testing correction
    de_enrichments <- de_enrichments %>%
      mutate(
        hypergeom_p_adj = p.adjust(hypergeom_p, method = "BH"),
        significant = hypergeom_p_adj < 0.05 & de_in_module >= 3 & enrichment_ratio >= 1.5
      ) %>%
      arrange(hypergeom_p)
    
  } else {
    de_enrichments <- NULL
  }
} else {
  de_enrichments <- NULL
}

# Module Connectivity Analysis --------------------------------------------

# Load module correlation matrix
module_connectivity <- core_results$ME_correlations
module_connectivity <- abs(module_connectivity)
diag(module_connectivity) <- 0

# Summarize connectivity patterns
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

# Save Results -------------------------------------------------------------

# Analysis results
write_csv(sex_difference_tests, file.path("10_dimorphism_analysis", "sex_difference_tests.csv"))
write_csv(interaction_tests, file.path("10_dimorphism_analysis", "interaction_effect_tests.csv"))
write_csv(module_connectivity_summary, file.path("10_dimorphism_analysis", "module_connectivity_analysis.csv"))

# Save DE enrichments if available
if (!is.null(de_enrichments)) {
  write_csv(de_enrichments, file.path("10_dimorphism_analysis", "de_module_enrichments.csv"))
}

# Module expression profiles for plotting
write_csv(module_profiles, file.path("10_module_profiles", "module_expression_profiles.csv"))

# Save module connectivity matrix
write_csv(as_tibble(module_connectivity, rownames = "module"), 
          file.path("10_dimorphism_analysis", "module_correlation_matrix.csv"))

# Summary statistics
analysis_summary <- tibble(
  Analysis = c("Sex Differences", "Interactions", "Module Connectivity"),
  N_Tests = c(nrow(sex_difference_tests), nrow(interaction_tests), 
              nrow(module_connectivity_summary)),
  N_Significant = c(sum(sex_difference_tests$significant), 
                    sum(interaction_tests$significant_interaction),
                    sum(module_connectivity_summary$hub_module))
)

if (!is.null(de_enrichments)) {
  analysis_summary <- analysis_summary %>%
    add_row(Analysis = "DE Enrichments", 
            N_Tests = nrow(de_enrichments),
            N_Significant = sum(de_enrichments$significant))
}

write_csv(analysis_summary, file.path("10_dimorphism_analysis", "analysis_summary.csv"))

# Save objects for downstream analysis
saveRDS(list(
  module_profiles = module_profiles,
  sex_difference_tests = sex_difference_tests,
  interaction_tests = interaction_tests,
  de_enrichments = de_enrichments,
  module_connectivity = module_connectivity,
  module_connectivity_summary = module_connectivity_summary,
  analysis_summary = analysis_summary
), file.path("10_intermediate_data", "dimorphism_analysis_results.rds"))

# Session Info ------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "./10_dimorphism_analysis/10_session_info.txt")