################################################################################
#
# Title: Gene Expression Overlap Analysis and Functional Enrichment
#        Colorado Potato Beetle (CPB) Diapause RNA-seq
#
# Author: Alex Torson
# Institution: USDA-ARS
# Email: Alex.Torson@usda.gov
#
# Description: Analyzes overlapping differentially expressed genes between male
#              and female CPB during diapause, creates Euler diagrams, and
#              performs GO and KEGG pathway enrichment analyses for biologically
#              meaningful gene categories. SIMPLIFIED VERSION with corrected
#              KEGG analysis data preparation.
#
# Dependencies: Differential expression results from DESeq2.R
#
# Input files:
#   - 01_DESeq2/Diapause_vs_NonDiapause_Males_DE_Transcripts.csv
#   - 01_DESeq2/Diapause_vs_NonDiapause_Females_DE_Transcripts.csv
#   - 01_DESeq2/Int_DE_Transcripts.csv
#   - 00_Annotation/cpb_diapause_annotation.csv
#   - 00_Annotation/cpb_diapause_go_omicsbox_export_one_seq_per_row.txt
#   - 00_Annotation/cpb_diapause_KAAS.txt
#
# Output files:
#   - 07_DE_Overlaps/CPB_DE_Euler_*.png (Euler diagrams)
#   - 07_DE_Overlaps/Annotated_Overlap_Lists/*.csv (annotated gene lists)
#   - 07_DE_Overlaps/GO_Enrichments/*.csv (GO enrichment results)
#   - 07_DE_Overlaps/KEGG_Enrichments/*.csv (KEGG enrichment results)
#
################################################################################

# Package Loading ---------------------------------------------------------

library(tidyverse)
library(eulerr)
library(pathview)
library(clusterProfiler)
library(topGO)
library(gage)
library(cowplot)

# Working Directory and Output Setup --------------------------------------

setwd("/Users/alextorson/Library/CloudStorage/OneDrive-USDA/Torson_Lab/CPB_Diapause_RNA-seq/")

# Standardized Color Palette and Theme ------------------------------------

diapause_colors <- c("Diapause" = "#0072B2", "NonDiapause" = "#E69F00")
sex_colors <- c("Male" = "#00AFBB", "Female" = "#009E73")

theme_cpb <- theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'plain', size = 10),
    axis.text = element_text(face = "plain", size = 10),
    axis.title = element_text(face = "plain", size = 10),
    axis.title.x = element_text(margin = margin(t = 8, r = 20, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 6, b = 0, l = 0)),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# Data Loading and Setup --------------------------------------------------

if (!dir.exists("./07_DE_Overlaps/")) {
  dir.create("./07_DE_Overlaps/", recursive = TRUE)
}

# Load CPB transcriptome annotations
cpb_annot <- read.csv('./00_Annotation/cpb_diapause_annotation.csv',
                      stringsAsFactors = TRUE)

# Load differential expression results from each comparison
male_de <- read.csv("./01_DESeq2/Diapause_vs_NonDiapause_Males_DE_Transcripts.csv")
female_de <- read.csv("./01_DESeq2/Diapause_vs_NonDiapause_Females_DE_Transcripts.csv")
int_de <- read.csv("./01_DESeq2/Int_DE_Transcripts.csv")

# Prepare gene lists for comparison ---------------------------------------

male_genes <- male_de$qry_transcript_id
female_genes <- female_de$qry_transcript_id
int_genes <- int_de$qry_transcript_id

# Calculate Euler diagram regions ------------------------------------------

# Find all unique sets for Euler regions
all_de <- unique(c(male_genes, female_genes, int_genes))

# Calculate Euler regions
male_female_both <- intersect(male_genes, female_genes)
shared_no_interaction <- setdiff(male_female_both, int_genes)
shared_with_interaction <- intersect(male_female_both, int_genes)
total_shared <- c(shared_no_interaction, shared_with_interaction)
male_only_euler <- setdiff(male_genes, total_shared)
female_only_euler <- setdiff(female_genes, total_shared)
male_female_only_euler <- shared_no_interaction
shared_with_interaction_euler <- shared_with_interaction

# Count genes in each region
regions <- c(
  "Male" = length(male_only_euler),
  "Female" = length(female_only_euler),
  "Male&Female" = length(male_female_only_euler),
  "Male&Female&Int" = length(shared_with_interaction_euler)
)

# Create Euler plot -------------------------------------------------------

euler_plot <- euler(regions)

# Define colors for each region
region_colors <- c(
  "Female" = "#009E73",
  "Male" = "#00AFBB", 
  "Interaction" = "#CC79A7",
  "Overlap" = "#000000"
)

# Generate full Euler plot
pdf("./07_DE_Overlaps/CPB_DE_Euler.pdf", 
    width = 3.54, height = 3.54)

(plot(euler_plot, 
      quantities = list(type = c("counts"),
                        cex = 0.8,
                        col = "black",
                        font = 1),
      labels = list(cex = 1.0, col = "black", font = 2),
      edges = list(col = "black", lwd = 1.5),
      fills = list(fill = region_colors, alpha = 0.9),
      main = ""))

dev.off()

# Create blank version (no text) for Illustrator editing
pdf("./07_DE_Overlaps/CPB_DE_Euler_BLANK_for_Illustrator.pdf", 
    width = 3.54, height = 3.54)

(plot(euler_plot, 
      quantities = FALSE,
      labels = FALSE,
      edges = list(col = "black", lwd = 1.5),
      fills = list(fill = region_colors, alpha = 0.9),
      main = ""))

dev.off()

# Overlap Category Creation Using Euler Regions -------------------------

profile_overlap <- list("Male_Unique" = male_only_euler,
                        "Female_Unique" = female_only_euler,
                        "Shared_No_Interaction" = male_female_only_euler,
                        "Shared_With_Interaction" = shared_with_interaction_euler)

# Annotated Gene Lists Creation -------------------------------------------

# Create directory for annotated lists
dir.create("./07_DE_Overlaps/Annotated_Overlap_Lists/", 
           showWarnings = FALSE, recursive = TRUE)

# Extract fold change data for joining
male_fc <- male_de %>% 
  dplyr::select(qry_transcript_id, male_log2FoldChange = log2FoldChange)

female_fc <- female_de %>% 
  dplyr::select(qry_transcript_id, female_log2FoldChange = log2FoldChange)

int_fc <- int_de %>% 
  dplyr::select(qry_transcript_id, interaction_log2FoldChange = log2FoldChange)

# Create annotated gene lists for each overlap category
annotated_overlaps <- lapply(1:length(profile_overlap), function(i) {
  # Get current overlap group name and transcript IDs
  overlap_name <- names(profile_overlap)[i]
  overlap_df <- data.frame(qry_transcript_id = profile_overlap[[i]], 
                           stringsAsFactors = FALSE)
  
  # Start with annotations
  annotated_df <- overlap_df %>%
    inner_join(cpb_annot, by = "qry_transcript_id")
  
  # Add relevant fold changes based on overlap category
  if (overlap_name %in% c("Male_Unique", "Shared_No_Interaction", "Shared_With_Interaction")) {
    annotated_df <- annotated_df %>% 
      left_join(male_fc, by = "qry_transcript_id")
  }
  
  if (overlap_name %in% c("Female_Unique", "Shared_No_Interaction", "Shared_With_Interaction")) {
    annotated_df <- annotated_df %>% 
      left_join(female_fc, by = "qry_transcript_id")
  }
  
  if (overlap_name == "Shared_With_Interaction") {
    annotated_df <- annotated_df %>% 
      left_join(int_fc, by = "qry_transcript_id")
  }
  
  # Write annotated results to CSV
  write.csv(annotated_df,
            file = paste0("./07_DE_Overlaps/Annotated_Overlap_Lists/",
                          overlap_name, "_annotated_with_FC.csv"),
            row.names = FALSE)
  
  return(annotated_df)
})

# Add names to annotated overlaps list
names(annotated_overlaps) <- names(profile_overlap)

# Save annotated overlaps for future use
save(annotated_overlaps, file = "./07_DE_Overlaps/annotated_overlaps.RData")

# Male vs Female Fold Change Correlation Analysis -------------------------

# Read the data
shared_no_int_data <- read.csv("./07_DE_Overlaps/Annotated_Overlap_Lists/Shared_No_Interaction_annotated_with_FC.csv")
shared_with_int_data <- read.csv("./07_DE_Overlaps/Annotated_Overlap_Lists/Shared_With_Interaction_annotated_with_FC.csv")

# Calculate Spearman correlations
cor_no_int <- cor.test(shared_no_int_data$male_log2FoldChange, 
                       shared_no_int_data$female_log2FoldChange, 
                       method = "spearman")

cor_with_int <- cor.test(shared_with_int_data$male_log2FoldChange, 
                         shared_with_int_data$female_log2FoldChange, 
                         method = "spearman")

# Shared genes correlation plot - No Interaction
(shared_cor_plot <- ggplot(shared_no_int_data, aes(x = male_log2FoldChange, y = female_log2FoldChange)) +
    geom_point(size = 1.2, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "#0072B2") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
    xlim(-30, 30) + ylim(-40, 40) +
    labs(x = expression("Male log"[2] * " Fold Change"),
         y = expression("Female log"[2] * " Fold Change"),
         title = "Male vs. Female overlap (No Int.)") +
    annotate("text", x = 30, y = -40, 
             label = paste0("R = ", round(cor_no_int$estimate, 2), ", p < 2.2e-16"), 
             hjust = 1, size = 3.5) +
    theme_cpb)

# Interaction genes correlation plot - With Interaction  
(int_cor_plot <- ggplot(shared_with_int_data, aes(x = male_log2FoldChange, y = female_log2FoldChange)) +
    geom_point(size = 1.2, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "#0072B2") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
    #xlim(-30, 30) + 
    ylim(-40, 20) +
    labs(x = expression("Male log"[2] * " Fold Change"),
         y = expression("Female log"[2] * " Fold Change"),
         title = "Male vs. Female overlap (With Int.)") +
    annotate("text", x = -3, y = -40, 
             label = paste0("R = ", round(cor_with_int$estimate, 2), ", p = 1.5e-10"), 
             hjust = 0, size = 3.5) +
    theme_cpb)

# Combine plots
(correlation_panel <- plot_grid(shared_cor_plot, int_cor_plot, ncol = 2, labels = c("A.", "B.")))

# Save
ggsave("./07_DE_Overlaps/Male_vs_Female_Correlation_Panel.pdf", 
       correlation_panel, width = 7, height = 3.5)

# GO Enrichment Setup -----------------------------------------------------

# Load transcript-level GO annotations
geneID2GO <- read.delim(
  file = "./00_Annotation/cpb_diapause_go_omicsbox_export_one_seq_per_row.txt",
  sep = "\t", header = TRUE, fill = TRUE, na.strings = "",
  stringsAsFactors = TRUE) %>%
  dplyr::rename(qry_gene_id = Sequence.Name,
                go_id = Annotation.GO.ID)

# Create list of GO terms for each gene
geneID2GOList <- strsplit(as.character(geneID2GO$go_id), split = ",")
names(geneID2GOList) <- geneID2GO$qry_gene_id

# Define gene universe for enrichment analysis
geneUniverse <- as.character(geneID2GO$qry_gene_id)

# GO Enrichment Analysis Function -----------------------------------------

# Create GO enrichment directory
dir.create('./07_DE_Overlaps/GO_Enrichments/', showWarnings = FALSE, recursive = TRUE)
setwd('./07_DE_Overlaps/GO_Enrichments/')

intersection_GO_func <- function(genesOfInterest) {
  # Extract category label and gene list
  label <- names(genesOfInterest)
  genesOfInterest <- unlist(genesOfInterest, use.names = FALSE)
  
  # Create directory for annotated GO results
  dir.create("./Annotated_GO_Results/", showWarnings = FALSE, recursive = TRUE)
  
  # Create gene list for topGO analysis
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse
  
  # Biological Process Enrichment
  myGOdata_BP <- new("topGOdata",
                     description = paste(label, "BP"),
                     ontology = "BP",
                     allGenes = geneList,
                     annot = annFUN.gene2GO,
                     gene2GO = geneID2GOList)
  
  resultFisher_BP <- runTest(myGOdata_BP, algorithm = "parentChild", statistic = "fisher")
  
  allRes_BP <- GenTable(myGOdata_BP,
                        p_value = resultFisher_BP,
                        orderBy = "p_value",
                        ranksOf = "p_value",
                        topNodes = 200)
  
  # Add gene IDs for each enriched term
  allRes_BP$gene_id <- sapply(allRes_BP$GO.ID, function(x) {
    genes <- genesInTerm(myGOdata_BP, x)
    paste(genes[[1]][genes[[1]] %in% genesOfInterest], collapse = ", ")
  })
  
  # Filter by p-value cutoff before saving
  allRes_BP_filtered <- allRes_BP %>%
    filter(as.numeric(p_value) <= 0.01)
  
  write.csv(allRes_BP_filtered, row.names = FALSE,
            file = paste0(label, "_BP_Encrichment.csv"))
  
  # Molecular Function Enrichment  
  myGOdata_MF <- new("topGOdata",
                     description = paste(label, "MF"),
                     ontology = "MF",
                     allGenes = geneList,
                     annot = annFUN.gene2GO,
                     gene2GO = geneID2GOList)
  
  resultFisher_MF <- runTest(myGOdata_MF, algorithm = "parentChild", statistic = "fisher")
  
  allRes_MF <- GenTable(myGOdata_MF,
                        p_value = resultFisher_MF,
                        orderBy = "p_value",
                        ranksOf = "p_value",
                        topNodes = 200)
  
  # Add gene IDs for each enriched term
  allRes_MF$gene_id <- sapply(allRes_MF$GO.ID, function(x) {
    genes <- genesInTerm(myGOdata_MF, x)
    paste(genes[[1]][genes[[1]] %in% genesOfInterest], collapse = ", ")
  })
  
  # Filter by p-value cutoff before saving
  allRes_MF_filtered <- allRes_MF %>%
    filter(as.numeric(p_value) <= 0.01)
  
  write.csv(allRes_MF_filtered, row.names = FALSE,
            file = paste0(label, "_MF_Encrichment.csv"))
  
  # Cellular Component Enrichment
  myGOdata_CC <- new("topGOdata",
                     description = paste(label, "CC"),
                     ontology = "CC",
                     allGenes = geneList,
                     annot = annFUN.gene2GO,
                     gene2GO = geneID2GOList)
  
  resultFisher_CC <- runTest(myGOdata_CC, algorithm = "parentChild", statistic = "fisher")
  
  allRes_CC <- GenTable(myGOdata_CC,
                        p_value = resultFisher_CC,
                        orderBy = "p_value",
                        ranksOf = "p_value",
                        topNodes = 200)
  
  # Add gene IDs for each enriched term
  allRes_CC$gene_id <- sapply(allRes_CC$GO.ID, function(x) {
    genes <- genesInTerm(myGOdata_CC, x)
    paste(genes[[1]][genes[[1]] %in% genesOfInterest], collapse = ", ")
  })
  
  # Filter by p-value cutoff before saving
  allRes_CC_filtered <- allRes_CC %>%
    filter(as.numeric(p_value) <= 0.01)
  
  write.csv(allRes_CC_filtered, row.names = FALSE,
            file = paste0(label, "_CC_Encrichment.csv"))
}

# Run GO enrichment for each overlap category
lapply(1:length(profile_overlap), function(i)
  intersection_GO_func(profile_overlap[i]))

# Return to main directory
setwd("/Users/alextorson/Library/CloudStorage/OneDrive-USDA/Torson_Lab/CPB_Diapause_RNA-seq/")

# KEGG Enrichment Setup --------------------------------------------------

# Load KEGG ortholog mappings from KAAS annotation
cpb_kegg <- read.table(
  "./00_Annotation/cpb_diapause_KAAS.txt",
  header = FALSE,
  colClasses = c("character", "factor"),
  sep = "",
  fill = TRUE,
  col.names = c("qry_transcript_id", "kegg_id"),
  na.strings = ""
) %>%
  drop_na()

# Merge KEGG IDs with gene annotations
cpb_annot <- cpb_annot %>%
  left_join(cpb_kegg, by = "qry_transcript_id")

# Create kegg_genes mapping (gene to KEGG ortholog)
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

# KEGG Enrichment Analysis Function ---------------------------------------

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
  plot_dir_ko <- file.path(original_wd, "07_DE_Overlaps/KEGG_Enrichments/Mapping_Plots_KO_IDs", label)
  plot_dir_enzymes <- file.path(original_wd, "07_DE_Overlaps/KEGG_Enrichments/Mapping_Plots_Enzymes", label)
  list_dir <- file.path(original_wd, "07_DE_Overlaps/KEGG_Enrichments/Lists", label)
  
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
  
  # Reset to enzyme directory and generate enzyme visualizations
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
  
  # Set placeholder values for commented pathview results
  ko_results <- rep(NULL, length(de_pathway_ids))
  enzyme_results <- rep(NULL, length(de_pathway_ids))  
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
  # Extract pathway information for genes in enriched pathways
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

# Prepare DE data for KEGG analysis - SIMPLIFIED AND CORRECTED -----------

# Create DE dataframes with proper fold change data preservation
de_list <- list(
  "Male_Unique" = male_de %>% 
    filter(qry_transcript_id %in% male_only_euler),
  
  "Female_Unique" = female_de %>% 
    filter(qry_transcript_id %in% female_only_euler),
  
  "Shared_No_Interaction" = male_de %>% 
    filter(qry_transcript_id %in% male_female_only_euler),
  
  "Shared_With_Interaction" = male_de %>% 
    filter(qry_transcript_id %in% shared_with_interaction_euler)
)

# Run KEGG Analysis on All Comparisons ------------------------------------

# Apply KEGG analysis to all differential expression datasets
kegg_results <- imap(de_list, perform_kegg_analysis)

# Print final results summary
cat("\n=== Final KEGG Analysis Summary ===\n")
iwalk(kegg_results, ~ {
  cat("\nCategory:", .y, "\n")
  cat("  Pathways found:", .x$n_pathways, "\n")
  cat("  KO plots generated:", .x$ko_plots, "\n")
  cat("  Enzyme plots generated:", .x$enzyme_plots, "\n")
  cat("  Enriched pathways:", .x$n_enriched, "\n")
})


# Create Enrichment Summary Files -----------------------------------------

# Read GO BP enrichment results and add comparison ID
female_unique_go <- read.csv("./07_DE_Overlaps/GO_Enrichments/Female_Unique_BP_Encrichment.csv") %>%
  mutate(comparison = "Female_Unique")

male_unique_go <- read.csv("./07_DE_Overlaps/GO_Enrichments/Male_Unique_BP_Encrichment.csv") %>%
  mutate(comparison = "Male_Unique")

shared_no_int_go <- read.csv("./07_DE_Overlaps/GO_Enrichments/Shared_No_Interaction_BP_Encrichment.csv") %>%
  mutate(comparison = "Shared_No_Interaction")

shared_with_int_go <- read.csv("./07_DE_Overlaps/GO_Enrichments/Shared_With_Interaction_BP_Encrichment.csv") %>%
  mutate(comparison = "Shared_With_Interaction")

# Combine all GO results
combined_go_bp <- rbind(female_unique_go, male_unique_go, shared_no_int_go, shared_with_int_go)

# Save combined GO results
write.csv(combined_go_bp, "./07_DE_Overlaps/GO_Enrichments/DE_Overlap_GO_Enrichment_Results_Full.csv", row.names = FALSE)

# KEGG Pathway Results ====================================================

# Read KEGG enrichment results and add comparison ID
female_unique_kegg <- read.csv("./07_DE_Overlaps/KEGG_Enrichments/Lists/Female_Unique/Female_Unique_enrichment_results.csv") %>%
  mutate(comparison = "Female_Unique")

male_unique_kegg <- read.csv("./07_DE_Overlaps/KEGG_Enrichments/Lists/Male_Unique/Male_Unique_enrichment_results.csv") %>%
  mutate(comparison = "Male_Unique")

shared_no_int_kegg <- read.csv("./07_DE_Overlaps/KEGG_Enrichments/Lists/Shared_No_Interaction/Shared_No_Interaction_enrichment_results.csv") %>%
  mutate(comparison = "Shared_No_Interaction")

shared_with_int_kegg <- read.csv("./07_DE_Overlaps/KEGG_Enrichments/Lists/Shared_With_Interaction/Shared_With_Interaction_enrichment_results.csv") %>%
  mutate(comparison = "Shared_With_Interaction")

# Combine all KEGG results
combined_kegg <- rbind(female_unique_kegg, male_unique_kegg, shared_no_int_kegg, shared_with_int_kegg)

# Save combined KEGG results
write.csv(combined_kegg, "./07_DE_Overlaps/KEGG_Enrichments/DE_Overlap_KEGG_Enrichment_Results_Full.csv", row.names = FALSE)

# Session Info ------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "./07_DE_Overlaps/07_session_info.txt")