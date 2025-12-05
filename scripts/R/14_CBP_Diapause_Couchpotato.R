################################################################################
#
# Title: Corrected Couchpotato Expression Analysis with DESeq2 Fold Changes
#        Colorado Potato beetle (CPB) Diapause RNA-seq
#
# Author: Alex Torson
# Institution: USDA ARS
# Email: Alex.Torson@usda.gov
#
# Description: Creates corrected visualization using actual DESeq2 fold changes
#              with mean ± SD points and individual sample points
#
################################################################################

# Package Loading ---------------------------------------------------------

library(tidyverse)
library(ggplot2)
library(scales)
library(ggh4x)  # For nested axis labels

# Set working directory
setwd("/Users/alextorson/Library/CloudStorage/OneDrive-USDA/Torson_Lab/CPB_Diapause_RNA-seq/")

# Create output directory
dir.create("./14_Couchpotato_Analysis/", showWarnings = FALSE, recursive = TRUE)

# Alex's Theme with legend moved to right
Alex_Theme <- theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, vjust = 0, face = 'plain', size = 14),
    plot.subtitle = element_text(size = 12, hjust = 0.5, face = "italic", color = "black"),
    axis.text = element_text(face = "plain", size = 11),
    axis.title = element_text(face = "plain", size = 12),
    axis.title.x = element_text(margin = margin(t = 8, r = 20, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 6, b = 0, l = 0)),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black", size = 0.5)
  )

# Define standard colorblind-friendly palette
cbPalette <- c("#000000", "#E69F00", "#00AFBB", "#009E73", 
               "#0072B2", "#D55E00", "#CC79A7")

# Data Loading and Preparation --------------------------------------------

# Load normalized counts
norm_counts <- read_csv("./01_DESeq2/cpb_diapause_normalized_counts.csv") %>%
  rename(transcript_id = names(.)[1])

# Load sample information  
sample_info <- read_csv("./01_DESeq2/cpb_sample_information.csv")

# Load DESeq2 results for actual fold changes
male_deseq <- read_csv("./01_DESeq2/Diapause_vs_NonDiapause_Males_All_Transcripts.csv")
female_deseq <- read_csv("./01_DESeq2/Diapause_vs_NonDiapause_Females_All_Transcripts.csv")

# Extract couchpotato data
couchpotato_transcripts <- c("TU8769", "TU8773", "TU8772")

# Get actual fold changes from DESeq2 results
male_fc <- male_deseq %>%
  filter(qry_transcript_id %in% couchpotato_transcripts) %>%
  dplyr::select(qry_transcript_id, log2FoldChange, padj) %>%
  mutate(sex = "male")

female_fc <- female_deseq %>%
  filter(qry_transcript_id %in% couchpotato_transcripts) %>%
  dplyr::select(qry_transcript_id, log2FoldChange, padj) %>%
  mutate(sex = "female")

# Combine fold changes
all_fc <- bind_rows(male_fc, female_fc) %>%
  mutate(
    significant = padj < 0.05,
    fc_text = case_when(
      is.na(log2FoldChange) ~ "NA",
      significant ~ paste0(ifelse(log2FoldChange > 0, "+", ""), 
                           round(log2FoldChange, 1), " log2FC"),
      TRUE ~ paste0(round(log2FoldChange, 1), " log2FC (n.s.)")
    )
  )

# Prepare expression data for plotting with nested structure
couchpotato_expr <- norm_counts %>%
  filter(transcript_id %in% couchpotato_transcripts) %>%
  pivot_longer(cols = -transcript_id, names_to = "Sample_ID", values_to = "normalized_counts") %>%
  left_join(sample_info, by = "Sample_ID") %>%
  mutate(
    condition = paste(sex, status, sep = "_"),
    # Create nested factor for axis labels
    status_nested = factor(status, levels = c("Diapause", "NonDiapause")),
    sex_nested = factor(sex, levels = c("female", "male"))
  )

# Calculate mean and SD for error bars with nested structure
summary_stats <- couchpotato_expr %>%
  group_by(transcript_id, sex, status) %>%
  summarise(
    mean_counts = mean(normalized_counts, na.rm = TRUE),
    sd_counts = sd(normalized_counts, na.rm = TRUE),
    n = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    condition = paste(sex, status, sep = "_"),
    se_counts = sd_counts / sqrt(n),
    # Create nested factors
    status_nested = factor(status, levels = c("Diapause", "NonDiapause")),
    sex_nested = factor(sex, levels = c("female", "male"))
  )

# Create the main plot with points for mean ± SD
p_main <- ggplot() +
  # Mean points with error bars (SD) - BEHIND individual points
  geom_errorbar(data = summary_stats,
                aes(x = interaction(status_nested, sex_nested), ymin = pmax(0, mean_counts - sd_counts), 
                    ymax = mean_counts + sd_counts),
                color = "grey70", width = 0.3, size = 1, alpha = 0.3) +
  geom_point(data = filter(summary_stats, mean_counts > 0),
             aes(x = interaction(status_nested, sex_nested), y = mean_counts, color = status),
             size = 4, shape = 18, alpha = 0.3) +  # Diamond shape for means - only when mean > 0
  # Individual sample points - IN FRONT
  geom_point(data = couchpotato_expr, 
             aes(x = interaction(status_nested, sex_nested), y = normalized_counts, fill = status),
             position = position_jitter(width = 0.2, seed = 123),
             size = 2, alpha = 0.6, shape = 21, stroke = 0.5) +
  facet_wrap(~ transcript_id, scales = "fixed", ncol = 3,
             labeller = labeller(transcript_id = function(x) x)) +
  scale_fill_manual(
    values = c("Diapause" = "#E69F00", "NonDiapause" = "#0072B2"),
    name = "Diapause Status:",
    labels = c("Diapause", "Non-diapause")
  ) +
  scale_color_manual(
    values = c("Diapause" = "#E69F00", "NonDiapause" = "#0072B2"),
    name = "Diapause Status:",
    labels = c("Diapause", "Non-diapause")
  ) +
  scale_y_continuous(
    labels = comma_format(),
    breaks = pretty_breaks(n = 5)
  ) +
  scale_x_discrete(
    labels = c("Diapause.female" = "Diapause", "NonDiapause.female" = "Non-diapause",
               "Diapause.male" = "Diapause", "NonDiapause.male" = "Non-diapause")
  ) +
  # Add vertical line to separate male vs female
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "grey60", alpha = 0.7) +
  # Add manual sex labels at the top of each panel
  annotate("text", x = 1.5, y = 2600, label = "Female", fontface = "plain", size = 4) +
  annotate("text", x = 3.5, y = 2600, label = "Male", fontface = "plain", size = 4) +
  labs(
    x = NULL,
    y = "Normalized Counts"
  ) +
  Alex_Theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    strip.text = element_text(size = 10),
    strip.background = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9)
  ) +
  guides(fill = guide_legend(override.aes = list(alpha = 1, size = 3)),
         color = "none")  # Hide color legend to avoid duplication

# Add fold change annotations from DESeq2 results
annotation_data <- all_fc %>%
  mutate(
    x_pos = case_when(
      sex == "female" ~ 1.5,
      sex == "male" ~ 3.5
    ),
    y_pos = 2400,  # Move to top of graph
    transcript_id = qry_transcript_id  # Add this for proper faceting
  )

p_annotated <- p_main +
  geom_text(data = annotation_data,
            aes(x = x_pos, y = y_pos, label = fc_text),
            color = "black",
            inherit.aes = FALSE,
            fontface = "plain",
            size = 3,
            hjust = 0.5,
            vjust = 0)

# Print the figure
print(p_annotated)

# Save the figure at full width (180mm) for paper
ggsave("./14_Couchpotato_Analysis/Figure_S3_Couchpotato_Expression.pdf", 
       p_annotated, width = 180/25.4, height = 115/25.4)

# Session Info
sessionInfo()