################################################################################
#
# Title: KEGG Pathway Enrichment Summary Plots
#        Colorado Potato Beetle (CPB) Diapause RNA-seq
#
# Author: Alex Torson
# Institution: USDA-ARS
# Email: Alex.Torson@usda.gov
#
# Description: Creates bubble plots summarizing KEGG pathway enrichments
#              across different differential expression comparisons including
#              sex, diapause status, and interaction effects.
#
# Dependencies: Completed KEGG enrichment analyses from CPB_Diapause_KEGG_Enrichment.R
#
# Input files:
#   - ./05_KEGG_Enrichment/Lists/d_vs_nd_female/d_vs_nd_female_enrichment_results.csv
#   - ./05_KEGG_Enrichment/Lists/d_vs_nd_male/d_vs_nd_male_enrichment_results.csv
#   - ./05_KEGG_Enrichment/Lists/int/int_enrichment_results.csv
#
# Output files:
#   - ./06_KEGG_Summary/combined_KEGG_Enrichment_Plot.png
#
################################################################################

# Load required packages
library(tidyverse)

# Standardized color palettes
comparison_colors <- c(
  "Females" = "#009E73",
  "Males" = "#00AFBB", 
  "Interaction" = "#CC79A7"
)

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
    legend.position = "bottom"
  )

# Working Directory and Output Setup --------------------------------------

setwd('/Users/alextorson/Library/CloudStorage/OneDrive-USDA/Torson_Lab/CPB_Diapause_RNA-seq/')

if (!dir.exists("./06_KEGG_Summary/")) {
  dir.create("./06_KEGG_Summary/", recursive = TRUE)
}

# Change to KEGG summary directory for output
setwd("./06_KEGG_Summary/")

# Load and combine all KEGG enrichment results ---------------------------

# Load each dataset with appropriate comparison labels
d_vs_nd_female_kegg <- read.csv(file = "../05_KEGG_Enrichment/Lists/d_vs_nd_female/d_vs_nd_female_enrichment_results.csv") %>%
  mutate(comparison = "Females")

d_vs_nd_male_kegg <- read.csv(file = "../05_KEGG_Enrichment/Lists/d_vs_nd_male/d_vs_nd_male_enrichment_results.csv") %>%
  mutate(comparison = "Males")

int_kegg <- read.csv(file = "../05_KEGG_Enrichment/Lists/int/int_enrichment_results.csv") %>%
  mutate(comparison = "Interaction")

# Combine all datasets
combined_kegg <- bind_rows(d_vs_nd_female_kegg, d_vs_nd_male_kegg, int_kegg)

# Get all unique pathways across all comparisons
all_pathways <- unique(combined_kegg$Description)

# Create a complete grid of all pathways × all comparisons
complete_grid <- expand_grid(
  Description = all_pathways,
  comparison = c("Females", "Males", "Interaction")
)

# Left join to fill in missing combinations with NA
combined_kegg_complete <- left_join(complete_grid, combined_kegg, 
                                    by = c("Description", "comparison"))

# Set factor levels for consistent ordering
combined_kegg_complete$comparison <- factor(combined_kegg_complete$comparison, 
                                            levels = c("Females", "Males", "Interaction"))

# Create combined plot
(combined_plot <- ggplot(combined_kegg_complete, aes(x = comparison, y = Description)) +
    theme_cpb +
    geom_point(aes(size = -log(qvalue), color = comparison), alpha = 0.75) +
    scale_y_discrete(
      labels = function(x) {
        is_long <- nchar(x) > 60
        x[is_long] <- paste0(substr(x[is_long], 1, 60), "...")
        x
      }
    ) +
    scale_color_manual(values = comparison_colors, na.value = "transparent") +
    scale_size_continuous(name = "-log(q-value)") +
    labs(title = "KEGG Pathway Enrichments", 
         y = "KEGG Pathway", 
         x = "Comparison") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    guides(color = "none"))

# Save the combined plot
ggsave("combined_KEGG_Enrichment_Plot.png", 
       plot = combined_plot,
       width = 8, height = 20)

# Session Info ------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "06_session_info.txt")