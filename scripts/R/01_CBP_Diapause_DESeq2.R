################################################################################
#
# Title: Differential Gene Expression Analysis of Diapausing and Non-diapausing
#        Colorado Potato Beetle (CPB) Males and Females
#
# Author: Alex Torson
# Institution: USDA-ARS
# Email: Alex.Torson@usda.gov
#
# Description: This script performs differential expression analysis using DESeq2
#              to compare gene expression patterns between diapausing and
#              non-diapausing Colorado potato beetles, with separate analyses
#              for males and females. The analysis includes interaction effects
#              to determine if diapause responses differ between sexes.
#
# Input files:
#   - transcript_count_matrix.csv: Transcript count matrix from StringTie
#   - cpb_sample_information.csv: Sample metadata
#   - cpb_diapause_annotation.csv: Gene annotations
#
# Output files:
#   - Various CSV files with differential expression results
#   - MA plots for each comparison
#   - Normalized count matrix
#
################################################################################

# Package Loading ---------------------------------------------------------

library(DESeq2)
library(tidyverse)
library(scales)

# Plotting Theme and Color Palette Setup ----------------------------------

# Significance colors for MA plots
sig_colors <- c("TRUE" = "purple", "FALSE" = "black")

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

# Data Loading and Preprocessing ------------------------------------------

# Set working directory
setwd(
  "/Users/alextorson/Library/CloudStorage/OneDrive-USDA/Torson_Lab/CPB_Diapause_RNA-seq/"
)

# Load CPB transcriptome annotations
cpb_annot <- read.csv('./00_Annotation/cpb_diapause_annotation.csv',
                      stringsAsFactors = TRUE)

# Load gene count data produced by StringTie/prepDE.py script
countData <- as.matrix(
  read.csv(
    "./01_DESeq2/transcript_count_matrix.csv",
    row.names = "transcript_id",
    check.names = FALSE
  )
)

# Load sample information
colData <- read.csv("./01_DESeq2/cpb_sample_information.csv", row.names = 1)

# Check all sample IDs in colData are also in CountData and match their orders
all(rownames(colData) %in% colnames(countData)) # Must be TRUE to move on

countData <- countData[, rownames(colData)]

all(rownames(colData) == colnames(countData)) # Must be TRUE to move on

# Change to DESeq2 output directory
setwd("./01_DESeq2/")

# DESeq2 Object Creation and Analysis -------------------------------------

# Note: The following section creates the DESeq2 object and runs the analysis.
# This code is commented out because the analysis has already been run and
# the results are saved. Uncomment to re-run the full analysis.

# Create DESeq2 object with interaction design
# The design tests for:
# - Main effect of sex
# - Main effect of diapause status
# - Interaction between sex and diapause status
# dds <- DESeqDataSetFromMatrix(
#   countData = countData,
#   colData = colData,
#   design = ~ sex + status + sex:status
# )
#
# # Set reference levels for comparisons
# # Non-diapause as reference for status
# # Male as reference for sex
# dds$status <- relevel(dds$status, "NonDiapause")
# dds$sex <- relevel(dds$sex, "male")
#
# # Pre-filtering: remove genes with very low counts
# # Genes with total count < 50 across all samples are removed
# # This improves computational efficiency and reduces multiple testing burden
# keep <- rowSums(counts(dds)) >= 50
# dds <- dds[keep,]
#
# # Run DESeq2 analysis
# dds <- DESeq(dds)
#
# # Export normalized counts
# dds_norm_out <- counts(dds, normalized = TRUE)
# write.csv(dds_norm_out, file = 'cpb_diapause_normalized_counts.csv')
#
# # Save DESeq2 object for future use
# save(dds, file = "cpb_diapause_deseq2_object")

# Load pre-computed DESeq2 object
load("cpb_diapause_deseq2_object")

# Display available result coefficients
resultsNames(dds)

# Differential Expression Analysis: Males ---------------------------------

# Extract results for diapause effect in males (main effect)
male_res <- results(dds, name = "status_Diapause_vs_NonDiapause", alpha = 0.05)

# Apply log-fold change shrinkage using adaptive shrinkage (ashr)
# This reduces noise in low-count genes and improves effect size estimates
male_resShrink <- lfcShrink(dds,
                            coef = "status_Diapause_vs_NonDiapause",
                            res = male_res,
                            type = "ashr")

# Summary of results
summary(male_resShrink)

# Convert to data frame for downstream analysis
male_resShrink_df <- as.data.frame(male_resShrink) %>%
  tibble::rownames_to_column("qry_transcript_id")

# Export results
write.csv(male_resShrink_df, file = 'Diapause_vs_NonDiapause_Males_All_Transcripts.csv', row.names = FALSE)

# Export only significantly differentially expressed genes
male_resShrink_df %>%
  subset(padj <= 0.05) %>%
  write.csv(file = 'Diapause_vs_NonDiapause_Males_DE_Transcripts.csv', row.names = FALSE)

# Add gene annotations
male_resShrink_df_annot <- left_join(male_resShrink_df, cpb_annot, by = "qry_transcript_id")

# Export annotated DE genes
male_resShrink_df_annot %>%
  subset(padj <= 0.05) %>%
  write.csv(file = 'Diapause_vs_NonDiapause_Males_DE_Transcripts_Annotated.csv', row.names = FALSE)

# Create MA plot for males
# Replace NA p-values with 0.99 for visualization (standard DESeq2 approach)
male_plot_data <- male_resShrink_df
male_plot_data[is.na(male_plot_data)] <- 0.99

(male_MA_plot <- ggplot(male_plot_data, aes(x = baseMean, y = log2FoldChange)) +
    theme_cpb +
    geom_point(aes(colour = padj < 0.05), size = 1, alpha = 0.6) +
    scale_colour_manual(name = 'padj < 0.05', values = setNames(c(sig_colors["TRUE"], sig_colors["FALSE"]), c(TRUE, FALSE))) +
    scale_x_log10(
      breaks = trans_breaks("log10", function(x)
        10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    ggtitle("Diapause vs. Non-diapause: Males") +
    ylab(expression("Fold change (log"[2] * ")")) +
    xlab('Mean of normalized counts') +
    theme(
      legend.position.inside = c(0.085, 0.16),
      legend.box.background = element_rect(colour = "black", linewidth = 0.5)
    ) +
    guides(colour = guide_legend(override.aes = list(size = 3))))

# Save MA plot
ggsave(
  plot = male_MA_plot,
  filename = "Diapause_vs_NonDiapause_Males_MA_Plot.png",
  width = 7,
  height = 4
)

# Differential Expression Analysis: Females -------------------------------

# Extract results for diapause effect in females
# This combines the main effect + interaction term to get the total effect in females
female_res <- results(dds, list(
  c(
    "status_Diapause_vs_NonDiapause",
    "sexfemale.statusDiapause"
  )
), alpha = 0.05)

# Apply log-fold change shrinkage
female_resShrink <- lfcShrink(
  dds,
  coef = c(
    "status_Diapause_vs_NonDiapause",
    "sexfemale.statusDiapause"
  ),
  res = female_res,
  type = "ashr"
)

# Summary of results
summary(female_resShrink)

# Convert to data frame
female_resShrink_df <- as.data.frame(female_resShrink) %>%
  tibble::rownames_to_column("qry_transcript_id")

# Export results
write.csv(female_resShrink_df, file = 'Diapause_vs_NonDiapause_Females_All_Transcripts.csv', row.names = FALSE)

# Export DE genes
female_resShrink_df %>%
  subset(padj <= 0.05) %>%
  write.csv(file = 'Diapause_vs_NonDiapause_Females_DE_Transcripts.csv', row.names = FALSE)

# Add annotations and export
female_resShrink_df_annot <- left_join(female_resShrink_df, cpb_annot, by = "qry_transcript_id")

female_resShrink_df_annot %>%
  subset(padj <= 0.05) %>%
  write.csv(file = 'Diapause_vs_NonDiapause_Females_DE_Transcripts_Annotated.csv', row.names = FALSE)

# Create MA plot for females
female_plot_data <- female_resShrink_df
female_plot_data[is.na(female_plot_data)] <- 0.99

(female_MA_plot <- ggplot(female_plot_data, aes(x = baseMean, y = log2FoldChange)) +
    theme_cpb +
    geom_point(aes(colour = padj < 0.05), size = 1, alpha = 0.6) +
    scale_colour_manual(name = 'padj < 0.05', values = setNames(c(sig_colors["TRUE"], sig_colors["FALSE"]), c(TRUE, FALSE))) +
    scale_x_log10(
      breaks = trans_breaks("log10", function(x)
        10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    ggtitle("Diapause vs. Non-diapause: Females") +
    ylab(expression("Fold change (log"[2] * ")")) +
    xlab('Mean of normalized counts') +
    theme(
      legend.position.inside = c(0.085, 0.16),
      legend.box.background = element_rect(colour = "black", linewidth = 0.5)
    ) +
    guides(colour = guide_legend(override.aes = list(size = 3))))

# Save MA plot
ggsave(
  plot = female_MA_plot,
  filename = "Diapause_vs_NonDiapause_Females_MA_Plot.png",
  width = 7,
  height = 4
)

# Interaction Analysis  ---------------------------------------------------

# Test for genes where the diapause response differs between sexes
int_res <- results(dds, name = "sexfemale.statusDiapause", alpha = 0.05)

# Apply shrinkage
int_resShrink <- lfcShrink(dds,
                           coef = "sexfemale.statusDiapause",
                           res = int_res,
                           type = "ashr")

# Summary
summary(int_resShrink)

# Convert and export
int_resShrink_df <- as.data.frame(int_resShrink) %>%
  tibble::rownames_to_column("qry_transcript_id")

write.csv(int_resShrink_df, file = 'Int_All_Transcripts.csv', row.names = FALSE)

int_resShrink_df %>%
  subset(padj < 0.05) %>%
  write.csv(file = 'Int_DE_Transcripts.csv', row.names = FALSE)

# Add annotations
int_resShrink_df_annot <- left_join(int_resShrink_df, cpb_annot, by = "qry_transcript_id")

int_resShrink_df_annot %>%
  subset(padj < 0.05) %>%
  write.csv(file = 'Int_DE_Transcripts_Annotated.csv', row.names = FALSE)

# Create MA plot for interaction effects
int_plot_data <- int_resShrink_df
int_plot_data[is.na(int_plot_data)] <- 0.99

(interaction_MA_plot <- ggplot(int_plot_data, aes(x = baseMean, y = log2FoldChange)) +
    theme_cpb +
    geom_point(aes(colour = padj < 0.05), size = 1, alpha = 0.6) +
    scale_colour_manual(name = 'padj < 0.05', values = setNames(c(sig_colors["TRUE"], sig_colors["FALSE"]), c(TRUE, FALSE))) +
    scale_x_log10(
      breaks = trans_breaks("log10", function(x)
        10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    ggtitle("Sex × Diapause Interaction Effects") +
    ylab(expression("Fold change (log"[2] * ")")) +
    xlab('Mean of normalized counts') +
    theme(
      legend.position.inside = c(0.085, 0.16),
      legend.box.background = element_rect(colour = "black", linewidth = 0.5)
    ) +
    guides(colour = guide_legend(override.aes = list(size = 3))))

ggsave(
  plot = interaction_MA_plot,
  filename = "Int_MA_Plot.png",
  width = 7,
  height = 4
)

# Session Info ------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "./01_DESeq2/01_session_info.txt")
