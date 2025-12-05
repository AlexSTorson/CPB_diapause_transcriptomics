################################################################################
#
# Title: Differential Expression Overlap Enrichment Summary Plots
#        Colorado Potato Beetle (CPB) Diapause RNA-seq
#
# Author: Alex Torson
# Institution: USDA-ARS
# Email: Alex.Torson@usda.gov
#
# Description: Creates bubble plots and chord plots summarizing KEGG pathway and
#              GO biological process enrichments across different differential
#              expression overlap categories. Uses semantic similarity reduction
#              to minimize GO term redundancy and includes manual curation of
#              pathway descriptions for cleaner visualizations.
#
# Dependencies: Completed overlap enrichment analyses from CPB_Diapause_DE_Overlap_Venn.R
#
# Input files:
#   - 07_DE_Overlaps/KEGG_Enrichments/Lists/*/enrichment_results.csv
#   - 07_DE_Overlaps/GO_Enrichments/*_BP_Encrichment.csv
#   - 08_DE_Overlap_Enrichment_Summary/KEGG_Enrichments/Master_KEGG_Descriptions_SH.csv
#   - 08_DE_Overlap_Enrichment_Summary/GO_Enrichments/DE_Overlap_GO_SH_Terms.csv
#
# Output files:
#   - 08_DE_Overlap_Enrichment_Summary/KEGG_Enrichments/DE_Overlap_KEGG_Enrichment_Plot.pdf
#   - 08_DE_Overlap_Enrichment_Summary/GO_Enrichments/DE_Overlap_GO_Enrichment_Plot.pdf
#   - 08_DE_Overlap_Enrichment_Summary/KEGG_Enrichments/DE_Overlap_KEGG_Enrichment_Chord_Plot.pdf
#   - 08_DE_Overlap_Enrichment_Summary/GO_Enrichments/DE_Overlap_GO_Enrichment_Chord_Plot.pdf
#
################################################################################

# Package Loading ---------------------------------------------------------

library(tidyverse)
library(circlize)
library(rrvgo)
library(org.Dm.eg.db)

# Working Directory and Output Setup --------------------------------------

setwd("/Users/alextorson/Library/CloudStorage/OneDrive-USDA/Torson_Lab/CPB_Diapause_RNA-seq/")

# Create output directories
output_dirs <- c(
  "./08_DE_Overlap_Enrichment_Summary/",
  "./08_DE_Overlap_Enrichment_Summary/KEGG_Enrichments/",
  "./08_DE_Overlap_Enrichment_Summary/GO_Enrichments/"
)

walk(output_dirs, ~ if (!dir.exists(.x)) dir.create(.x, recursive = TRUE))

# Color Palette and Theme Setup -------------------------------------------

# Standardized color palette for comparisons
comparison_colors <- c(
  "Female" = "#009E73",
  "Male" = "#00AFBB", 
  "Interaction" = "#CC79A7",
  "Overlap" = "#000000"
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
    legend.position = "none"
  )

# KEGG Data Loading and Processing -----------------------------------------

# Load KEGG enrichment results from overlap analyses
male_kegg <- read.csv("./07_DE_Overlaps/KEGG_Enrichments/Lists/Male_Unique/Male_Unique_enrichment_results.csv") %>%
  mutate(comp = "Male")

female_kegg <- read.csv("./07_DE_Overlaps/KEGG_Enrichments/Lists/Female_Unique/Female_Unique_enrichment_results.csv") %>%
  mutate(comp = "Female")

overlap_kegg <- read.csv("./07_DE_Overlaps/KEGG_Enrichments/Lists/Shared_No_Interaction/Shared_No_Interaction_enrichment_results.csv") %>%
  mutate(comp = "Overlap")

# Combine KEGG datasets
kegg_df <- bind_rows(male_kegg, female_kegg, overlap_kegg)

# Save combined KEGG results
write.csv(kegg_df, 
          file = "./08_DE_Overlap_Enrichment_Summary/KEGG_Enrichments/DE_Overlap_KEGG_Enrichment_Results.csv",
          row.names = FALSE)

# Load master KEGG pathway descriptions for visualization
master_kegg_descriptions <- read.csv("./08_DE_Overlap_Enrichment_Summary/KEGG_Enrichments/Master_KEGG_Descriptions_SH.csv")

# GO Data Loading and Processing -------------------------------------------

# Load GO enrichment results from overlap analyses
male_go <- read.csv("./07_DE_Overlaps/GO_Enrichments/Male_Unique_BP_Encrichment.csv") %>%
  mutate(comp = "Male")

female_go <- read.csv("./07_DE_Overlaps/GO_Enrichments/Female_Unique_BP_Encrichment.csv") %>%
  mutate(comp = "Female")

overlap_go <- read.csv("./07_DE_Overlaps/GO_Enrichments/Shared_No_Interaction_BP_Encrichment.csv") %>%
  mutate(comp = "Overlap")

interaction_go <- read.csv("./07_DE_Overlaps/GO_Enrichments/Shared_With_Interaction_BP_Encrichment.csv") %>%
  mutate(comp = "Interaction")

# Combine GO datasets
go_df <- bind_rows(male_go, female_go, overlap_go, interaction_go)

# Save combined GO results
write.csv(go_df, 
          file = "./08_DE_Overlap_Enrichment_Summary/GO_Enrichments/DE_Overlap_GO_Enrichment_Results_Full.csv",
          row.names = FALSE)

# Optional GO Level Filtering ----------------------------------------------

# Load GO level information for optional filtering
go_levels <- read.delim(
  file = "./00_Annotation/cpb_diapause_blast2go_go_propagation.txt",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
) %>%
  dplyr::rename(qry_gene_id = Seq.Name, GO_ID = GO, Level = Level) %>%
  filter(Level >= 4 & Level <= 7) %>%
  dplyr::select(GO_ID, Level) %>%
  distinct()

# Option 1: Filter by GO level (4-7) for more specific terms
go_df_filtered <- go_df %>%
  filter(GO.ID %in% go_levels$GO_ID)

# Option 2: Use all GO terms (uncomment to use instead of level filtering)
# go_df_filtered <- go_df

# KEGG Bubble Plot --------------------------------------------------------

# Add short descriptions for cleaner visualization
kegg_df_annotated <- kegg_df %>%
  left_join(master_kegg_descriptions, by = "ID") %>%
  mutate(label_with_id = paste0(Sh_description, " (", ID, ")"))

(kegg_bubble_plot <- ggplot(kegg_df_annotated, aes(y = label_with_id, x = comp)) +
    theme_cpb +
    theme(
      legend.position = "right", 
      legend.background = element_rect(fill = "white", color = "black"),
      legend.margin = margin(5, 5, 5, 5)
    ) +
    geom_point(aes(size = -log(qvalue), color = comp), alpha = 0.75) +
    scale_color_manual(values = comparison_colors, guide = "none") +
    scale_size_continuous(name = "-log(q-value)") +
    labs(y = "KEGG pathway", x = "Comparison"))

# Save KEGG bubble plot
ggsave("./08_DE_Overlap_Enrichment_Summary/KEGG_Enrichments/DE_Overlap_KEGG_Enrichment_Plot.pdf",
       plot = kegg_bubble_plot,
       width = 6, height = 12)

# GO Term Redundancy Reduction --------------------------------------------

# Apply rrvgo semantic similarity reduction for cleaner chord plots
go_terms_for_slimming <- go_df_filtered %>%
  group_by(GO.ID) %>%
  dplyr::summarise(min_pvalue = min(p_value, na.rm = TRUE), .groups = 'drop')

# Create named vector of p-values for rrvgo
go_pvalues <- setNames(go_terms_for_slimming$min_pvalue, go_terms_for_slimming$GO.ID)

# Calculate semantic similarity matrix
simMatrix <- calculateSimMatrix(names(go_pvalues),
                                orgdb = "org.Dm.eg.db",
                                ont = "BP",
                                method = "Rel")

# Reduce GO term redundancy using semantic similarity
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores = go_pvalues,
                                threshold = 0.7,
                                orgdb = "org.Dm.eg.db")

# Filter to representative GO terms only
go_terms_to_keep <- reducedTerms$go
go_df_slimmed <- go_df_filtered %>%
  filter(GO.ID %in% go_terms_to_keep)

# Save slimmed GO results
write.csv(go_df_slimmed, 
          file = "./08_DE_Overlap_Enrichment_Summary/GO_Enrichments/DE_Overlap_GO_Enrichment_Results_Slimmed.csv",
          row.names = FALSE)

# GO Short-hand Terms Setup -----------------------------------------------

# Load manually curated short GO term descriptions
go_names <- read.csv("./08_DE_Overlap_Enrichment_Summary/GO_Enrichments/DE_Overlap_GO_SH_Terms.csv")

# Add short descriptions to slimmed GO data for visualization
go_df_annotated <- go_df_slimmed %>%
  left_join(go_names, by = "Term") %>%
  filter(!is.na(Sh_term)) %>%
  mutate(label_with_id = paste0(Sh_term, " (", GO.ID, ")"))

# GO Bubble Plot with Short Terms -----------------------------------------

(go_bubble_plot <- ggplot(go_df_annotated, aes(y = label_with_id, x = comp)) +
   theme_cpb +
   theme(
     legend.position = "right", 
     legend.background = element_rect(fill = "white", color = "black"),
     legend.margin = margin(5, 5, 5, 5)
   ) +
   geom_point(aes(size = -log(p_value), color = comp), alpha = 0.75) +
   scale_color_manual(values = comparison_colors, guide = "none") +
   scale_size_continuous(name = "-log(p-value)") +
   scale_x_discrete(limits = c("Female", "Male", "Overlap", "Interaction")) +
   labs(y = "GO Biological Process", x = "Comparison"))

# Save GO bubble plot
ggsave("./08_DE_Overlap_Enrichment_Summary/GO_Enrichments/DE_Overlap_GO_Enrichment_Plot.pdf",
       plot = go_bubble_plot,
       width = 8, height = 12)

# Chord Plot Data Preparation ---------------------------------------------

# Filter KEGG pathways to remove vertebrate-specific terms
vertebrate_kegg_keywords <- c(
  "phototransduction", "phototrans",
  "parathyroid", "PTH", 
  "aldosterone", "aldosterone synth",
  "calcium reabsorption", "calcium reabsorp", "reabsorption", "reabsorp",
  "collecting duct", "acid secretion", "acid secret",
  "bile secretion", "bile secret",
  "pancreatic secretion", "pancreatic secret",
  "renin secretion", "renin secret", 
  "adrenergic signaling", "adrenergic sig",
  "cardiac muscle contraction", "cardiac contract", "cardiac contraction",
  "thyroid hormone", "thyroid",
  "endocrine.*calcium",
  "vertebrate", "mammalian", "human"
)

# Filter vertebrate-specific KEGG pathways
vertebrate_kegg_pattern <- paste(vertebrate_kegg_keywords, collapse = "|")
kegg_enrich_filtered <- kegg_df %>%
  dplyr::select(comp, Description) %>%
  filter(!grepl(vertebrate_kegg_pattern, Description, ignore.case = TRUE))

# Add short descriptions for chord plot
kegg_enrich <- kegg_enrich_filtered %>%
  left_join(master_kegg_descriptions, by = "Description") %>%
  dplyr::select(comp, Sh_description) %>%
  filter(!is.na(Sh_description)) %>%
  distinct()

# Filter GO terms to remove vertebrate-specific terms
vertebrate_go_keywords <- c(
  "ocellus", "pigment.*ocellus", "pigmentation",
  "bone", "cartilage", "skeleton",
  "blood", "hemoglobin", "hematopoietic",
  "kidney", "renal", "nephron",
  "liver", "hepatic",
  "lung", "pulmonary",
  "heart", "cardiac",
  "brain", "neural crest", "neural tube",
  "limb", "appendage",
  "mammary", "lactation",
  "feather", "hair", "fur",
  "vertebrate", "chordate"
)

# Filter vertebrate-specific GO terms and prepare for chord plot
vertebrate_go_pattern <- paste(vertebrate_go_keywords, collapse = "|")
go_enrich_filtered <- go_df %>%
  dplyr::select(comp, Term) %>%
  filter(!grepl(vertebrate_go_pattern, Term, ignore.case = TRUE))

# Export GO terms for reference (optional)
go_df_annotated %>%
  dplyr::select(Term) %>%
  distinct() %>%
  write.csv(file = "./08_DE_Overlap_Enrichment_Summary/GO_Enrichments/DE_Overlap_GO_Terms_Level4to7_Filtered.csv",
            row.names = FALSE)

# Load manually curated short GO term descriptions
go_names <- read.csv("./08_DE_Overlap_Enrichment_Summary/GO_Enrichments/DE_Overlap_GO_SH_Terms.csv")

# Add short descriptions for chord plot
go_enrich <- go_enrich_filtered %>%
  left_join(go_names, by = "Term") %>%
  dplyr::select(comp, Sh_term) %>%
  filter(!is.na(Sh_term)) %>%
  distinct()

# Chord Plot Color Setup --------------------------------------------------

# Create color palettes for chord plots
kegg_othercol <- structure(rep("grey50", length(unique(kegg_enrich$Sh_description))))
names(kegg_othercol) <- unique(kegg_enrich$Sh_description)
kegg_grid.col <- c(comparison_colors, kegg_othercol)

go_othercol <- structure(rep("grey60", length(unique(go_enrich$Sh_term))))
names(go_othercol) <- unique(go_enrich$Sh_term)
go_grid.col <- c(comparison_colors, go_othercol)

# KEGG Chord Plot ---------------------------------------------------------

pdf("./08_DE_Overlap_Enrichment_Summary/KEGG_Enrichments/DE_Overlap_KEGG_Enrichment_Chord_Plot.pdf",
    width = 7.09, height = 7.09)

circos.par(start.degree = 270, 
           cell.padding = c(0, 0, 0, 0), 
           canvas.xlim = c(-1.5, 1.5), 
           canvas.ylim = c(-1.5, 1.5))

chordDiagram(kegg_enrich,
             big.gap = 5,
             small.gap = 2.2,
             annotationTrack = c("grid"),
             grid.col = kegg_grid.col,
             preAllocateTracks = list(track.height = 0.15))

circos.track(track.index = 1,
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter,
                           CELL_META$ylim[1],
                           CELL_META$sector.index,
                           facing = "clockwise",
                           niceFacing = TRUE,
                           adj = c(0, 0.5),
                           cex = 0.7)
             },
             bg.border = NA)

circos.clear()
dev.off()

# GO Chord Plot -----------------------------------------------------------

pdf("./08_DE_Overlap_Enrichment_Summary/GO_Enrichments/DE_Overlap_GO_Enrichment_Chord_Plot.pdf",
    width = 7.09, height = 7.09)

circos.par(start.degree = 270, 
           cell.padding = c(0, 0, 0, 0), 
           canvas.xlim = c(-1.5, 1.5), 
           canvas.ylim = c(-1.5, 1.5))

chordDiagram(go_enrich,
             big.gap = 5,
             small.gap = 2.2,
             annotationTrack = c("grid"),
             grid.col = go_grid.col,
             preAllocateTracks = list(track.height = 0.15))

circos.track(track.index = 1,
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter,
                           CELL_META$ylim[1],
                           CELL_META$sector.index,
                           facing = "clockwise",
                           niceFacing = TRUE,
                           adj = c(0, 0.5),
                           cex = 0.7)
             },
             bg.border = NA)

circos.clear()
dev.off()

# Session Info ------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "./08_DE_Overlap_Enrichment_Summary/08_session_info.txt")