################################################################################
#
# Title: GO Enrichment Summary Plots
#        Colorado Potato Beetle (CPB) Diapause RNA-seq
#
# Author: Alex Torson
# Institution: USDA-ARS
# Email: Alex.Torson@usda.gov
#
# Description: Creates bubble plots summarizing GO biological process enrichments
#              across different differential expression comparisons including
#              sex, diapause status, and interaction effects.
#
# Dependencies: Completed GO enrichment analyses from cpb_diapause_topGO.R
#
# Input files:
#   - ./GO_Analyses/diapause_f_vs_m_up_BP_Enrichment.csv
#   - ./GO_Analyses/diapause_f_vs_m_down_BP_Enrichment.csv
#   - ./GO_Analyses/d_vs_nd_female_up_BP_Enrichment.csv
#   - ./GO_Analyses/d_vs_nd_female_down_BP_Enrichment.csv
#   - ./GO_Analyses/d_vs_nd_male_up_BP_Enrichment.csv
#   - ./GO_Analyses/d_vs_nd_male_down_BP_Enrichment.csv
#   - ./GO_Analyses/int_up_BP_Enrichment.csv
#   - ./GO_Analyses/int_down_BP_Enrichment.csv
#
# Output files:
#   - ./GO_Analyses/diapause_f_vs_m_Enrichment_Plot.png
#   - ./GO_Analyses/d_vs_nd_female_Enrichment_Plot.png
#   - ./GO_Analyses/d_vs_nd_male_Enrichment_Plot.png
#   - ./GO_Analyses/int_Enrichment_Plot.png
#
################################################################################

# Load required packages
library(tidyverse)

# Standardized color palettes
diapause_colors <- c("Diapause" = "#0072B2", "NonDiapause" = "#E69F00")
sex_colors <- c("Male" = "#00AFBB", "Female" = "#009E73")

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

# Set working directory
setwd('/Users/alextorson/Library/CloudStorage/OneDrive-USDA/Torson_Lab/CPB_Diapause_RNA-seq/')

# Diapausing Females vs. Non-diapausing Females
d_vs_nd_female_up <- read.csv(file = "./03_GO_Analyses/d_vs_nd_female_up_BP_Enrichment.csv") %>%
  mutate(enrichment = "Diapause",
         ontology = "BP")

d_vs_nd_female_down <- read.csv(file = "./03_GO_Analyses/d_vs_nd_female_down_BP_Enrichment.csv") %>%
  mutate(enrichment = "NonDiapause",
         ontology = "BP")

d_vs_nd_female_down$p_value <- as.numeric(d_vs_nd_female_down$p_value)

d_vs_nd_female_go_df <- bind_rows(d_vs_nd_female_up, d_vs_nd_female_down)

ggplot(d_vs_nd_female_go_df, aes(y = Term, x = enrichment)) +
  theme_cpb +
  geom_point(aes(size = -log(p_value), color = enrichment), alpha = 0.75) +
  scale_y_discrete(
    labels = function(x) {
      is_long <- nchar(x) > 75
      x[is_long] <- paste0(substr(x[is_long], 1, 75), " ...")
      x
    }
  ) +
  scale_color_manual(values = diapause_colors) +
  labs(title = "Females", y = "Biological Process", x = "Diapause Status")

ggsave("./04_GO_Analyses_Plotting/d_vs_nd_female_Enrichmnt_Plot.png", width = 8, height = 20, type = "cairo-png")

# Diapausing Males vs. Non-diapausing Males
d_vs_nd_male_up <- read.csv(file = "./03_GO_Analyses/d_vs_nd_male_up_BP_Enrichment.csv") %>%
  mutate(enrichment = "Diapause",
         ontology = "BP")

d_vs_nd_male_down <- read.csv(file = "./03_GO_Analyses/d_vs_nd_male_down_BP_Enrichment.csv") %>%
  mutate(enrichment = "NonDiapause",
         ontology = "BP")

d_vs_nd_male_go_df <- bind_rows(d_vs_nd_male_up, d_vs_nd_male_down)

ggplot(d_vs_nd_male_go_df, aes(y = Term, x = enrichment)) +
  theme_cpb +
  geom_point(aes(size = -log(p_value), color = enrichment), alpha = 0.75) +
  scale_y_discrete(
    labels = function(x) {
      is_long <- nchar(x) > 75
      x[is_long] <- paste0(substr(x[is_long], 1, 75), " ...")
      x
    }
  ) +
  scale_color_manual(values = diapause_colors) +
  labs(title = "Males", y = "Biological Process", x = "Diapause Status")

ggsave("./04_GO_Analyses_Plotting/d_vs_nd_male_Enrichmnt_Plot.png", width = 8, height = 20, type = "cairo-png")

# How do males and females differ in their diapause response?
int_up <- read.csv(file = "./03_GO_Analyses/int_up_BP_Enrichment.csv") %>%
  mutate(enrichment = "Female",
         ontology = "BP")

int_down <- read.csv(file = "./03_GO_Analyses/int_down_BP_Enrichment.csv") %>%
  mutate(enrichment = "Male",
         ontology = "BP")

int_go_df <- bind_rows(int_up, int_down)

ggplot(int_go_df, aes(y = Term, x = enrichment)) +
  theme_cpb +
  geom_point(aes(size = -log(p_value), color = enrichment), alpha = 0.75) +
  scale_y_discrete(
    labels = function(x) {
      is_long <- nchar(x) > 75
      x[is_long] <- paste0(substr(x[is_long], 1, 75), " ...")
      x
    }
  ) +
  scale_color_manual(values = sex_colors) +
  labs(title = "Diapause response: Females vs. Males", y = "Biological Process", x = "Sex")

ggsave("./04_GO_Analyses_Plotting/int_Enrichmnt_Plot.png", width = 8, height = 20, type = "cairo-png")

# Session Info ------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "./04_GO_Analyses_Plotting/04_session_info.txt")
