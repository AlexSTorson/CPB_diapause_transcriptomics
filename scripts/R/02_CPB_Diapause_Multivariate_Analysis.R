################################################################################
#
# Title: Principal Component Analysis of Gene Expression During Diapause
#        Colorado Potato Beetle (CPB) Diapause RNA-seq
#
# Author: Alex Torson
# Institution: USDA-ARS
# Email: Alex.Torson@usda.gov
#
# Description: Principal component analysis of variance-stabilized gene expression
#              data to examine transcriptomic differences between diapausing and
#              non-diapausing Colorado potato beetles across sexes. Includes
#              k-means clustering of gene loadings for functional interpretation.
#
# Dependencies: DESeq2 object from previous differential expression analysis
#
# Input files:
#   - 01_DESeq2/cpb_diapause_deseq2_object (DESeq2 object)
#   - 00_Annotation/cpb_diapause_annotation.csv
#
# Output files:
#   - 02_Multivariate_Analysis/PCA_*.png (PCA plots)
#   - 02_Multivariate_Analysis/var_kmeans_clusters.csv (gene clusters)
#
################################################################################

# Package Loading ---------------------------------------------------------

library(DESeq2)
library(tidyverse)
library(cowplot)
library(scales)
library(factoextra)

# Working Directory and Output Setup --------------------------------------

setwd(
  "/Users/alextorson/Library/CloudStorage/OneDrive-USDA/Torson_Lab/CPB_Diapause_RNA-seq/"
)

if (!dir.exists("./02_Multivariate_Analysis/")) {
  dir.create("./02_Multivariate_Analysis/", recursive = TRUE)
}

# Standardized Color Palette and Theme ------------------------------------

# Diapause status colors (primary scheme)
diapause_colors <- c("Diapause" = "#0072B2",
                     "NonDiapause" = "#E69F00")

# Sex colors (assigned from variable clusters)
sex_colors <- c("Male" = "#00AFBB", "Female" = "#009E73")

# K-means cluster colors (these end up coding for sex and diapause status)
cluster_colors <- c("#0072B2", "#00AFBB", "#009E73", "#E69F00")

# Standardized theme
theme_cpb <- theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'plain', size = 10),
    axis.text = element_text(face = "plain", size = 10),
    axis.title = element_text(face = "plain", size = 10),
    axis.title.x = element_text(margin = margin(
      t = 8,
      r = 20,
      b = 0,
      l = 0
    )),
    axis.title.y = element_text(margin = margin(
      t = 0,
      r = 6,
      b = 0,
      l = 0
    )),
    panel.border = element_rect(
      fill = NA,
      colour = "black",
      linewidth = 0.5
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# Data Loading and Preparation --------------------------------------------

load("./01_DESeq2/cpb_diapause_deseq2_object")

cpb_annot <- read.csv('./00_Annotation/cpb_diapause_annotation.csv')

# Principal Component Analysis --------------------------------------------

# Variance stabilizing transformation for heteroskedastic count data
vst <- vst(dds, blind = FALSE)

# Select top variable genes for PCA
ntop <- 500
rv <- rowVars(assay(vst))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]

# Perform PCA on samples using top variable genes
pca <- prcomp(t(assay(vst)[select, ]))

# Extract variance explained by each component
pca_importance <- data.frame(summary(pca)$importance)
PC1 <- round(100 * pca_importance[2, 1])
PC2 <- round(100 * pca_importance[2, 2])
PC3 <- round(100 * pca_importance[2, 3])

# Set reference level for consistent plotting
vst$status <- relevel(vst$status, ref = "Diapause")

# PCA Visualization -------------------------------------------------------

# Scree plot
(scree_plot <- fviz_eig(pca) +
   ggtitle("PCA Scree Plot") +
   theme_cpb)

ggsave(
  "PCA_Scree.pdf",
  plot = scree_plot,
  width = 90,
  height = 90,
  units = "mm",
  path = './02_Multivariate_Analysis/'
)

# Sample PCA plot
(
  sample_pca <- fviz_pca_ind(
    pca,
    axes = c(1, 2),
    habillage = vst$sex,
    title = "",
    fill.ind = vst$status,
    addEllipses = FALSE,
    label = "none",
    repel = TRUE,
    mean = FALSE,
    alpha = 0.8,
    geom = "point",
    pointsize = 2.5
  ) +
    theme_cpb +
    xlab(paste0("PC1 (", PC1, "% of variance)")) +
    ylab(paste0("PC2 (", PC2, "% of variance)")) +
    scale_fill_manual(
      name = "Status:",
      labels = c("Diapause", "Non-diapause"),
      values = c(
        "Diapause" = "#0072B2",
        "NonDiapause" = "#E69F00"
      )
    ) +
    scale_color_manual(
      values = c("black", "black"),
      name = "Status:",
      labels = c("Diapause", "Non-diapause")
    ) +
    scale_shape_manual(
      values = c(21, 24),
      name = "Sex:",
      labels = c("Male", "Female")
    ) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.title = element_text(size = 10, face = "plain"),
      legend.text = element_text(size = 10),
      plot.margin = margin(
        t = 5,
        r = 5,
        b = 5,
        l = 5
      )
    ) +
    guides(
      fill = guide_legend(
        title = "Status:",
        override.aes = list(shape = 21, size = 2),
        title.position = "left",
        direction = "vertical",
        order = 1
      ),
      shape = guide_legend(
        title = "Sex:",
        override.aes = list(fill = "white", size = 2),
        title.position = "left",
        direction = "vertical",
        order = 2
      ),
      color = "none"
    )
)

ggsave(
  "PCA_Samples_Plot.pdf",
  plot = sample_pca,
  width = 100,
  height = 100,
  units = "mm",
  path = './02_Multivariate_Analysis/'
)

# Variable Loading Analysis -----------------------------------------------

# Extract PCA variable information (gene loadings)
pca_var <- get_pca_var(pca)

# Identify genes with above-average contribution to PC1 and PC2
PC12_contrib <- facto_summarize(pca,
                                element = "var",
                                result = "contrib",
                                axes = c(1, 2)) %>%
  filter(contrib > 0.2)

PC12_contrib_count <- nrow(PC12_contrib)

# K-means clustering of gene loadings in PC space
set.seed(123)
res.km <- kmeans(pca_var$coord, centers = 4, nstart = 25)
grp <- as.factor(res.km$cluster)

# Variable cluster plot
(
  var_cluster_plot <- fviz_pca_var(
    pca,
    col.var = grp,
    axes = c(1, 2),
    select.var = list(contrib = PC12_contrib_count),
    palette = cluster_colors,
    repel = TRUE,
    title = "",
    label = "none",
    alpha.var = 0.8
  ) +
    theme_cpb +
    xlab(paste0("PC1 (", PC1, "% of variance)")) +
    ylab(paste0("PC2 (", PC2, "% of variance)")) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10, face = "plain"),
      legend.direction = "horizontal",
      plot.margin = margin(
        t = 5,
        r = 5,
        b = 5,
        l = 5
      )
    ) +
    guides(
      color = guide_legend(
        title = "K-means cluster:",
        override.aes = list(size = 2, alpha = 1),
        ncol = 4,
        title.position = "top",
        title.hjust = 0
      )
    )
)

ggsave(
  "PCA_Variables_Clusters.pdf",
  plot = var_cluster_plot,
  width = 90,
  height = 100,
  units = "mm",
  path = './02_Multivariate_Analysis/'
)

# Publication Figure ------------------------------------------------------

(
  pca_figure <- plot_grid(
    sample_pca,
    var_cluster_plot,
    align = "hv",
    axis = "tblr",
    labels = c("A.", "B."),
    label_size = 12,
    label_fontface = "bold",
    hjust = -1,
    vjust = 2.5,
    ncol = 2,
    rel_widths = c(1, 1),
    scale = 0.95
  )
)

ggsave(
  "PCA_Figure_Publication.pdf",
  plot = pca_figure,
  width = 180,
  height = 120,
  units = "mm",
  dpi = 300,
  path = './02_Multivariate_Analysis/'
)

# Export Cluster Results --------------------------------------------------

# Extract gene coordinates in PC space
pca_coord <- as.data.frame(pca_var$coord) %>%
  rownames_to_column(var = "qry_transcript_id")

# Create cluster assignment dataframe
cluster_df <- data.frame(qry_transcript_id = rownames(pca_var$coord),
                         cluster = as.factor(res.km$cluster))

# Get detailed information for top contributing genes
pc12_top_contrib <- facto_summarize(
  pca,
  element = 'var',
  result = c('contrib', 'coord'),
  axes = 1:2,
  select = list(contrib = PC12_contrib_count)
) %>%
  dplyr::rename(qry_transcript_id = name)

# Combine cluster assignments with gene annotations
pc12_clusters <- pc12_top_contrib %>%
  left_join(cluster_df, by = "qry_transcript_id") %>%
  left_join(cpb_annot, by = "qry_transcript_id") %>%
  dplyr::select(cluster, contrib, qry_transcript_id, ref_transcript_name)

write.csv(pc12_clusters, row.names = FALSE, file = './02_Multivariate_Analysis/var_kmeans_clusters.csv')

# Session Info ------------------------------------------------------------

writeLines(capture.output(sessionInfo()),
           "./02_Multivariate_Analysis/02_session_info.txt")
