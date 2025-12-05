################################################################################
#
# Title: Publication Figures for WGCNA Sexual Dimorphism Analysis
#        Colorado Potato Beetle (CPB) Diapause Gene Co-expression Networks
#
# Author: Alex Torson
# Institution: USDA-ARS
# Email: Alex.Torson@usda.gov
#
# Description: Create figures for WGCNA analysis 
#              focusing on manuscript analyses: sexual dimorphism, interactions,
#              DE enrichment, and module connectivity.
#
# Dependencies: Must run 01_core_wgcna_analysis.R and 02_dimorphism_analysis.R first
#
# Input files:
#   - 01_intermediate_data/core_wgcna_results.rds
#   - 02_intermediate_data/dimorphism_analysis_results.rds
#
# Output files:
#   - 11_figures_PDFs/*.pdf (publication figures)
#
################################################################################

# Package Loading and Setup -----------------------------------------------

library(WGCNA)
library(tidyverse)
library(pheatmap)
library(corrplot)
library(RColorBrewer)
library(patchwork)
library(igraph)
library(ggraph)

set.seed(12345)

# Working Directory and Output Setup --------------------------------------

setwd("/project/igb_fargo/cpb_diapause_rnaseq/wgcna/")

graphics_dirs <- c("11_figures_PDFs")
walk(graphics_dirs, ~ if (!dir.exists(.x)) dir.create(.x, recursive = TRUE))

# Standardized Color Palette and Theme ------------------------------------

# Diapause status colors (primary scheme)
diapause_colors <- c("Diapause" = "#0072B2", "NonDiapause" = "#E69F00")

# Sex colors (assigned from variable clusters)
sex_colors <- c("Male" = "#00AFBB", "Female" = "#009E73")

# K-means cluster colors (for PCA variable plots)
cluster_colors <- c("#E69F00", "#00AFBB", "#009E73", "#0072B2")

# Significance colors for MA plots
sig_colors <- c("TRUE" = "#D55E00", "FALSE" = "black")

# Standardized theme
theme_cpb <- theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'plain', size = 12),
    axis.text = element_text(face = "plain", size = 10),
    axis.title = element_text(face = "plain", size = 10),
    axis.title.x = element_text(margin = margin(t = 8, r = 20, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 6, b = 0, l = 0)),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# Publication theme with legend
theme_cpb_legend <- theme_cpb + theme(legend.position = "bottom")

# Color palettes for publication
pub_colors <- list(
  sex = c("female" = "#009E73", "male" = "#00AFBB"),
  status = c("Diapause" = "#0072B2", "NonDiapause" = "#E69F00"),
  significance = c("significant" = "#D55E00", "not_significant" = "#BDC3C7")
)

# Load Previous Results ---------------------------------------------------

# Load core and dimorphism results
core_results <- readRDS(file.path("09_intermediate_data", "core_wgcna_results.rds"))
dimorphism_results <- readRDS(file.path("10_intermediate_data", "dimorphism_analysis_results.rds"))

# Extract key objects
expr_data <- core_results$expr_data
sample_info <- core_results$sample_info
module_colors <- core_results$module_colors
MEs <- core_results$module_eigengenes
hub_genes <- core_results$hub_genes
color_to_number_map_clean <- core_results$color_to_number_map
ME_cor <- core_results$ME_correlations
comprehensive_summary <- core_results$comprehensive_summary

module_profiles <- dimorphism_results$module_profiles
sex_difference_tests <- dimorphism_results$sex_difference_tests
interaction_tests <- dimorphism_results$interaction_tests
de_enrichments <- dimorphism_results$de_enrichments
module_connectivity_summary <- dimorphism_results$module_connectivity_summary

# No utility function needed - using ggsave directly

# Experimental Overview ---------------------------------------------------

experimental_overview <- sample_info %>%
  count(sex, status) %>%
  ggplot(aes(x = status, y = n, fill = sex)) +
  geom_col(position = "dodge", color = "black", alpha = 0.8, width = 0.7) +
  geom_text(aes(label = n), position = position_dodge(width = 0.7), 
            vjust = -0.3, fontface = "bold", size = 5) +
  scale_fill_manual(values = pub_colors$sex, name = "Sex") +
  labs(x = "Diapause Status", y = "Number of Samples") +
  theme_cpb_legend

module_distribution <- table(module_colors) %>%
  as.data.frame() %>%
  filter(module_colors != "grey") %>%
  ggplot(aes(x = Freq)) +
  geom_histogram(bins = 12, fill = pub_colors$sex[["male"]], alpha = 0.7, color = "black") +
  labs(x = "Number of Transcripts per Module", y = "Number of Modules") +
  theme_cpb

combined_overview <- experimental_overview | module_distribution
ggsave(file.path("11_figures_PDFs", "Experimental_Overview.pdf"), plot = combined_overview, width = 12, height = 6, dpi = 300)



# Sexually Dimorphic Module Expression Profiles --------------------------

# Get modules with significant sex differences
sex_dimorphic_modules <- sex_difference_tests %>%
  filter(significant) %>%
  pull(module_id) %>%
  unique()

if (length(sex_dimorphic_modules) > 0) {
  sex_dimorphic_data <- module_profiles %>%
    filter(module_id %in% sex_dimorphic_modules) %>%
    group_by(module_id, sex, status) %>%
    summarise(
      mean_expr = mean(expression, na.rm = TRUE),
      se_expr = sd(expression, na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    )
  
  # Add effect size information to facet labels
  module_effect_sizes <- sex_difference_tests %>%
    filter(significant) %>%
    group_by(module_id) %>%
    summarise(max_effect = max(abs(cohens_d)), .groups = 'drop') %>%
    arrange(desc(max_effect))
  
  # Plot top 15 sexually dimorphic modules (or all if fewer than 15)
  top_sex_modules <- module_effect_sizes %>%
    slice_head(n = 15) %>%
    pull(module_id)
  
  p_sex_dimorphic <- sex_dimorphic_data %>%
    filter(module_id %in% top_sex_modules) %>%
    left_join(module_effect_sizes, by = "module_id") %>%
    mutate(module_label = paste0(module_id, "\n(d = ", round(max_effect, 2), ")")) %>%
    ggplot(aes(x = status, y = mean_expr, color = sex, group = sex)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = mean_expr - se_expr, ymax = mean_expr + se_expr), 
                  width = 0.1, linewidth = 0.8) +
    facet_wrap(~reorder(module_label, -max_effect), scales = "free_y", ncol = 5) +
    scale_color_manual(values = pub_colors$sex, name = "Sex") +
    labs(x = "Diapause Status", y = "Module Expression (eigengene)") +
    theme_cpb_legend
  
  ggsave(file.path("11_figures_PDFs", "Sexually_Dimorphic_Module_Profiles.pdf"), plot = p_sex_dimorphic, width = 14, height = 10, dpi = 300)
}

# Sexual Dimorphism Analysis Heatmap --------------------------------------

sig_sex_differences <- sex_difference_tests %>%
  filter(significant, !is.na(cohens_d)) %>%
  arrange(desc(abs(cohens_d))) %>%
  slice_head(n = 20)

if (nrow(sig_sex_differences) > 0) {
  p_sex_dimorphism <- ggplot(sig_sex_differences, aes(x = status, y = reorder(module_id, abs(cohens_d)), 
                                                      fill = cohens_d)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", cohens_d)), 
              color = "black", size = 3.5, fontweight = "bold") +
    scale_fill_gradient2(
      low = pub_colors$sex[["male"]], mid = "white", high = pub_colors$sex[["female"]],
      midpoint = 0, name = "Cohen's d\n(Effect Size)",
      limits = c(-max(abs(sig_sex_differences$cohens_d)), max(abs(sig_sex_differences$cohens_d)))
    ) +
    labs(x = "Diapause Status", y = "Module") +
    theme_cpb_legend
  
  ggsave(file.path("11_figures_PDFs", "Sexual_Dimorphism_Heatmap.pdf"), plot = p_sex_dimorphism, width = 10, height = 8, dpi = 300)
}

# Interaction Effects ------------------------------------------------------

sig_interactions <- interaction_tests %>%
  filter(significant_interaction, !is.na(eta_squared)) %>%
  arrange(desc(eta_squared)) %>%
  slice_head(n = 15)

if (nrow(sig_interactions) > 0) {
  p_interactions <- ggplot(sig_interactions, aes(x = reorder(module_id, eta_squared), y = eta_squared)) +
    geom_col(fill = pub_colors$significance[["significant"]], alpha = 0.7, color = "black") +
    geom_text(aes(label = sprintf("%.3f", eta_squared)), 
              hjust = -0.1, size = 3.5) +
    coord_flip() +
    labs(x = "Module", y = "Effect Size (η²)") +
    theme_cpb
  
  ggsave(file.path("11_figures_PDFs", "Interaction_Effects.pdf"), plot = p_interactions, width = 10, height = 6, dpi = 300)
}

# Module Correlation Network (Top 20) -------------------------------------

top_20_modules <- module_profiles %>%
  filter(module_id != "Module_00") %>%
  group_by(module_id) %>%
  summarise(n_samples = n(), .groups = 'drop') %>%
  arrange(desc(n_samples)) %>%
  slice_head(n = 20) %>%
  pull(module_id)

top_20_colors <- module_profiles %>%
  filter(module_id %in% top_20_modules) %>%
  distinct(module_id, module_color) %>%
  mutate(me_name = paste0("ME", module_color)) %>%
  pull(me_name)

ME_cor_top20 <- ME_cor[top_20_colors, top_20_colors]

clean_names_top20 <- str_remove(rownames(ME_cor_top20), "ME")
clean_module_names_top20 <- map_chr(clean_names_top20, function(color) {
  if (color %in% names(color_to_number_map_clean)) {
    paste0("M", color_to_number_map_clean[color])
  } else {
    color
  }
})

rownames(ME_cor_top20) <- clean_module_names_top20
colnames(ME_cor_top20) <- clean_module_names_top20

pdf(file.path("11_figures_PDFs", "Module_Correlation_Network.pdf"), width = 10, height = 10)
corrplot(ME_cor_top20,
         method = "color",
         type = "upper", 
         order = "hclust",
         tl.cex = 1.2,
         tl.col = "black",
         tl.srt = 45,
         col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(200),
         mar = c(0, 0, 1, 0),
         cl.cex = 1.0)
dev.off()

# DE-Module Enrichment -----------------------------------------------------

if (!is.null(de_enrichments) && nrow(de_enrichments) > 0) {
  sig_enrichments <- de_enrichments %>%
    filter(significant) %>%
    arrange(desc(enrichment_ratio)) %>%
    slice_head(n = 30)
  
  plot_enrichments <- sig_enrichments %>%
    mutate(
      de_category_clean = case_when(
        str_detect(de_category, "Female.*Unique") ~ "Female\nUnique",
        str_detect(de_category, "Male.*Unique") ~ "Male\nUnique",
        str_detect(de_category, "Shared.*No.*Interaction") ~ "Shared\n(No Interaction)",
        str_detect(de_category, "Shared.*With.*Interaction") ~ "Shared\n(With Interaction)",
        str_detect(de_category, "Three.*Way.*Overlap") ~ "Shared\n(With Interaction)",
        str_detect(de_category, "Female.*Biased") ~ "Female\nBiased",
        str_detect(de_category, "Male.*Biased") ~ "Male\nBiased",
        TRUE ~ str_replace_all(de_category, "_", "\n")
      ),
      module_clean = str_remove(Module_ID, "Module_"),
      module_numeric = as.numeric(module_clean)
    )
  
  p_de_enrichment <- ggplot(plot_enrichments, aes(x = de_category_clean, y = reorder(module_clean, module_numeric), fill = enrichment_ratio)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f", enrichment_ratio)), 
              color = "black", size = 10/.pt) +
    scale_fill_gradient(low = "grey95", high = pub_colors$sex[["male"]], name = "Enrichment\nRatio") +
    labs(x = "DE Category", y = "Module") +
    theme_cpb_legend +
    theme(panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
          panel.grid.minor = element_line(color = "grey95", linewidth = 0.25),
          legend.position = "right")
  
  ggsave(file.path("11_figures_PDFs", "DE_Module_Enrichment.pdf"), plot = p_de_enrichment, width = 180/25.4, height = 6, dpi = 300)
}

# Module Network Visualization --------------------------------------------

top_size_modules <- comprehensive_summary %>%
  filter(Module_ID != "Module_00") %>%
  group_by(Module_ID) %>%
  summarise(module_size = n(), .groups = 'drop') %>%
  arrange(desc(module_size)) %>%
  slice_head(n = 20) %>%
  pull(Module_ID)

enriched_modules <- if (!is.null(de_enrichments)) {
  de_enrichments %>%
    filter(significant) %>%
    pull(Module_ID) %>%
    unique()
} else {
  character(0)
}

network_modules <- unique(c(top_size_modules, enriched_modules))

network_colors <- comprehensive_summary %>%
  filter(Module_ID %in% network_modules) %>%
  distinct(Module_ID, Module_Color) %>%
  mutate(me_name = paste0("ME", Module_Color))

ME_cor_network <- ME_cor[network_colors$me_name, network_colors$me_name]

adj_matrix <- abs(ME_cor_network)
adj_matrix[adj_matrix < 0.5] <- 0

node_data <- network_colors %>%
  left_join(
    comprehensive_summary %>%
      group_by(Module_ID) %>%
      summarise(module_size = n(), .groups = 'drop'),
    by = "Module_ID"
  ) %>%
  left_join(
    sex_difference_tests %>%
      group_by(module_id) %>%
      summarise(has_sex_diff = any(significant), .groups = 'drop'),
    by = c("Module_ID" = "module_id")
  ) %>%
  replace_na(list(has_sex_diff = FALSE)) %>%
  mutate(
    module_number = str_remove(Module_ID, "Module_"),
    module_size_capped = pmin(module_size, 1000),
    node_category = case_when(
      Module_ID %in% enriched_modules & has_sex_diff ~ "DE Enriched + Sex Diff",
      Module_ID %in% enriched_modules ~ "DE Enriched",
      has_sex_diff ~ "Sexually Dimorphic",
      TRUE ~ "Standard"
    )
  )

g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)

V(g)$module_id <- node_data$Module_ID
V(g)$module_size <- node_data$module_size_capped
V(g)$node_category <- node_data$node_category
V(g)$module_number <- node_data$module_number

set.seed(123)
p_network <- ggraph(g, layout = "fr") +
  # Grey edges with width scaling only (no alpha)
  geom_edge_link(aes(width = weight), color = "grey60") +
  geom_node_point(aes(size = module_size, fill = node_category), 
                  shape = 21, color = "black", alpha = 0.8) +
  geom_node_text(aes(label = module_number), size = 4, fontface = "plain", color = "black") +
  scale_size_continuous(name = "Module Size\n(transcripts)", range = c(8, 18)) +
  # Increased width range for better visibility
  scale_edge_width_continuous(name = "Correlation", range = c(0.5, 4)) +
  scale_fill_manual(
    name = "Module Type",
    values = c(
      "DE Enriched + Sex Diff" = pub_colors$sex[["female"]],
      "DE Enriched" = "#FF7F00",
      "Sexually Dimorphic" = pub_colors$sex[["male"]],
      "Standard" = "grey70"
    )
  ) +
  labs(x = "", y = "") +
  theme_cpb_legend +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin(t = 5),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  ) +
  guides(
    size = guide_legend(override.aes = list(fill = "grey70"), order = 1),
    fill = guide_legend(override.aes = list(size = 5), order = 2),
    edge_width = guide_legend(order = 3)
  )

ggsave(file.path("11_figures_PDFs", "Module_Network_Visualization.pdf"), plot = p_network, width = 180/25.4, height = 6.5, dpi = 300)

# Enhanced Network Visualization ------------------------------------------

peripheral_modules <- c("Module_59", "Module_60")
hub_modules <- c("Module_01", "Module_03", "Module_12", "Module_29", "Module_125")
interaction_modules <- c("Module_101", "Module_118", "Module_60", "Module_59", "Module_16")

node_data_enhanced <- node_data %>%
  mutate(
    cluster = case_when(
      Module_ID %in% peripheral_modules ~ "Peripheral Cluster",
      Module_ID %in% hub_modules ~ "Hub Module", 
      TRUE ~ "Core Network"
    )
  )

V(g)$cluster <- node_data_enhanced$cluster

set.seed(123)
p_network_enhanced <- ggraph(g, layout = "fr") +
  # Grey edges with width scaling only (no alpha)
  geom_edge_link(aes(width = weight), color = "grey60") +
  geom_node_point(aes(size = module_size, fill = node_category, shape = cluster), 
                  color = "black", alpha = 0.8, stroke = 1.5) +
  geom_node_text(aes(label = module_number), size = 4, fontface = "plain") +
  scale_shape_manual(
    name = "Network Region",
    values = c("Hub Module" = 23, "Core Network" = 21, "Peripheral Cluster" = 22),
    drop = FALSE
  ) +
  scale_size_continuous(name = "Module Size\n(transcripts)", range = c(8, 18)) +
  # Increased width range for better visibility
  scale_edge_width_continuous(name = "Correlation", range = c(0.5, 4)) +
  scale_fill_manual(
    name = "Module Type",
    values = c(
      "DE Enriched + Sex Diff" = pub_colors$sex[["female"]],
      "DE Enriched" = "#FF7F00",
      "Sexually Dimorphic" = pub_colors$sex[["male"]],
      "Standard" = "grey70"
    )
  ) +
  labs(x = "", y = "") +
  theme_cpb_legend +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin(t = 5),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  ) +
  guides(
    size = guide_legend(order = 1, override.aes = list(fill = "grey70", shape = 21)),
    fill = guide_legend(order = 2, override.aes = list(size = 5, shape = 21)),
    shape = guide_legend(order = 3, override.aes = list(size = 5, fill = "grey70")),
    edge_width = guide_legend(order = 4)
  )

ggsave(file.path("11_figures_PDFs", "Module_Network_Enhanced.pdf"), plot = p_network_enhanced, width = 10, height = 8, dpi = 300)

# Complete Network Visualization ------------------------------------------

all_modules <- comprehensive_summary %>%
  filter(Module_ID != "Module_00") %>%
  distinct(Module_ID, Module_Color) %>%
  mutate(me_name = paste0("ME", Module_Color))

ME_cor_full <- ME_cor[all_modules$me_name, all_modules$me_name]

adj_matrix_full <- abs(ME_cor_full)
adj_matrix_full[adj_matrix_full < 0.3] <- 0  # Manuscript threshold |r| > 0.3

node_data_full <- all_modules %>%
  left_join(
    comprehensive_summary %>%
      group_by(Module_ID) %>%
      summarise(module_size = n(), .groups = 'drop'),
    by = "Module_ID"
  ) %>%
  left_join(
    sex_difference_tests %>%
      group_by(module_id) %>%
      summarise(has_sex_diff = any(significant), .groups = 'drop'),
    by = c("Module_ID" = "module_id")
  ) %>%
  replace_na(list(has_sex_diff = FALSE)) %>%
  mutate(
    module_number = str_remove(Module_ID, "Module_"),
    is_hub = Module_ID %in% hub_modules,
    is_peripheral = Module_ID %in% peripheral_modules,
    is_interaction_enriched = Module_ID %in% interaction_modules,
    # Manuscript color scheme
    module_type = case_when(
      is_hub ~ "Hub", 
      is_peripheral ~ "Peripheral",
      is_interaction_enriched ~ "Interaction",
      TRUE ~ "Standard"
    )
  )

g_full <- graph_from_adjacency_matrix(adj_matrix_full, mode = "undirected", weighted = TRUE, diag = FALSE)

V(g_full)$module_id <- node_data_full$Module_ID
V(g_full)$module_size <- node_data_full$module_size
V(g_full)$module_type <- node_data_full$module_type
V(g_full)$module_number <- node_data_full$module_number
V(g_full)$is_peripheral <- node_data_full$is_peripheral

# Use Kamada-Kawai layout as specified in manuscript methods
layout_kk <- create_layout(g_full, layout = "kk")
set.seed(789)

p_network_complete <- ggraph(layout_kk) +
  # Draw all edges with manuscript threshold
  geom_edge_link(aes(width = weight), alpha = 0.3, color = "grey80") +
  
  # Draw nodes with manuscript color scheme
  geom_node_point(aes(size = module_size, fill = module_type, alpha = module_type), 
                  shape = 21, color = "black", stroke = 0.5) +
  
  # Add labels for key modules (hub, peripheral, interaction-enriched)
  geom_node_text(aes(label = ifelse(module_type %in% c("Hub", "Peripheral", "Interaction"), 
                                    module_number, "")), 
                 size = 3, fontface = "plain") +
  
  # Module size scaling (manuscript: n = 30-5,370)
  scale_size_continuous(
    name = "Module Size\n(transcripts)", 
    range = c(6, 16),
    trans = "sqrt"
  ) +
  
  scale_edge_width_continuous(
    name = "Correlation",
    range = c(0.1, 2)
  ) +
  
  # Manuscript color scheme
  scale_fill_manual(
    name = "Module Type",
    values = c(
      "Hub" = "#1F78B4",           # Blue for hub modules (01, 03, 12, 29, 125)
      "Peripheral" = "#E74C3C",    # Red for peripheral modules (59, 60)
      "Interaction" = "#FF7F00",   # Orange for other interaction modules (101, 118, 16)
      "Standard" = "grey85"        # Grey for standard modules
    )
  ) +
  
  scale_alpha_manual(
    values = c(
      "Hub" = 1,
      "Peripheral" = 1,
      "Interaction" = 0.9,
      "Standard" = 0.5
    ),
    guide = "none"
  ) +
  
  labs(x = "", y = "") +
  theme_cpb_legend +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.position = "right",
    plot.margin = margin(20, 20, 20, 20)
  ) +
  guides(
    size = guide_legend(order = 1, override.aes = list(fill = "grey70", shape = 21)),
    fill = guide_legend(
      title = "Module Type",
      order = 2, 
      override.aes = list(size = 6, shape = 21)
    ),
    edge_width = guide_legend(order = 3)
  )

ggsave(file.path("11_figures_PDFs", "Module_Network_Complete_KK.pdf"), plot = p_network_complete, width = 12, height = 10, dpi = 300)

# Calculate connectivity statistics for peripheral modules
peripheral_me_names <- all_modules %>%
  filter(Module_ID %in% peripheral_modules) %>%
  pull(me_name)

# Use the adjacency matrix from the complete network
peripheral_connectivity <- adj_matrix_full[peripheral_me_names, ] %>%
  as.data.frame() %>%
  rownames_to_column("module_me") %>%
  pivot_longer(-module_me, names_to = "connected_to", values_to = "correlation") %>%
  filter(module_me != connected_to, correlation > 0) %>%
  group_by(module_me) %>%
  summarise(
    n_connections = n(),
    mean_correlation = mean(correlation),
    max_correlation = max(correlation),
    .groups = 'drop'
  ) %>%
  left_join(
    all_modules %>% dplyr::select(me_name, Module_ID),
    by = c("module_me" = "me_name")
  ) %>%
  dplyr::select(Module_ID, n_connections, mean_correlation, max_correlation)

# Calculate average connectivity for non-peripheral modules
non_peripheral_me_names <- all_modules %>%
  filter(!Module_ID %in% peripheral_modules) %>%
  pull(me_name)

avg_connectivity <- map_dfr(non_peripheral_me_names, function(mod_me) {
  connections <- adj_matrix_full[mod_me, ]
  connections <- connections[names(connections) != mod_me & connections > 0]
  tibble(
    n_connections = length(connections),
    mean_correlation = ifelse(length(connections) > 0, mean(connections), 0)
  )
}) %>%
  summarise(
    avg_n_connections = mean(n_connections),
    avg_mean_correlation = mean(mean_correlation)
  )

# Save connectivity comparison
connectivity_comparison <- tibble(
  Module_Type = c(peripheral_connectivity$Module_ID, "Average (Non-Peripheral)"),
  Connections = c(peripheral_connectivity$n_connections, round(avg_connectivity$avg_n_connections, 1)),
  Mean_Correlation = c(peripheral_connectivity$mean_correlation, avg_connectivity$avg_mean_correlation)
)

write_csv(connectivity_comparison, file.path("11_figures_PDFs", "Module_Connectivity_Comparison.csv"))