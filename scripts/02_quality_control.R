# ============================================================
# Project: SPRY2 Reprogramming
# Script: 02_quality_control.R
# Purpose: Data preprocessing and transcriptomic quality control
# ============================================================
# ------------------------------------------------------------
# 1. Package setup
# ------------------------------------------------------------
required_packages <- c(
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "dendextend"
)
missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(dendextend)
# ------------------------------------------------------------
# 2. Project directories
# ------------------------------------------------------------
processed_dir <- file.path(
  "data",
  "processed"
)
main_figure_dir <- file.path(
  "results",
  "figures",
  "main"
)
supplementary_figure_dir <- file.path(
  "results",
  "figures",
  "supplementary"
)
dir.create(
  processed_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  main_figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  supplementary_figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
# ------------------------------------------------------------
# 3. Load cleaned expression matrix
# ------------------------------------------------------------
expression_file <- file.path(
  processed_dir,
  "expression_matrix_clean.rds"
)
if (!file.exists(expression_file)) {
  stop(
    "Cleaned expression matrix not found: ",
    expression_file,
    "\nRun 01_download_data.R first."
  )
}
expr <- readRDS(expression_file)
if (!is.matrix(expr) && !is.data.frame(expr)) {
  stop(
    "The expression object must be a matrix or data frame."
  )
}
expr <- as.matrix(expr)
if (!all(vapply(
  as.data.frame(expr),
  is.numeric,
  logical(1)
))) {
  stop(
    "All expression values must be numeric."
  )
}
message(
  "Input expression matrix: ",
  nrow(expr),
  " genes x ",
  ncol(expr),
  " samples"
)
# ------------------------------------------------------------
# 4. Expression distribution before log2 transformation
# ------------------------------------------------------------
expr_long <- stack(
  as.data.frame(expr)
)
colnames(expr_long) <- c(
  "Expression",
  "Sample"
)
p_raw <- ggplot(
  expr_long,
  aes(
    x = Sample,
    y = Expression
  )
) +
  geom_boxplot(
    fill = "#8ECAE6",
    outlier.size = 0.3
  ) +
  labs(
    title = "Expression distribution before log2 transformation",
    x = NULL,
    y = "RPKM"
  ) +
  theme_classic(
    base_size = 13
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )
ggsave(
  file.path(
    supplementary_figure_dir,
    "Expression_distribution_before_log2.png"
  ),
  p_raw,
  width = 9,
  height = 6,
  dpi = 300
)
# ------------------------------------------------------------
# 5. Log2 transformation
# ------------------------------------------------------------
expr_log2 <- log2(
  expr + 1
)
message(
  "Log2 transformation completed: log2(RPKM + 1)"
)
# ------------------------------------------------------------
# 6. Remove zero-variance genes
# ------------------------------------------------------------
gene_variance <- apply(
  expr_log2,
  1,
  var
)
zero_variance_genes <- sum(
  gene_variance == 0
)
message(
  "Genes with zero variance: ",
  zero_variance_genes
)
expr_log2_filtered <- expr_log2[
  gene_variance > 0,
  ,
  drop = FALSE
]
message(
  "Filtered expression matrix: ",
  nrow(expr_log2_filtered),
  " genes x ",
  ncol(expr_log2_filtered),
  " samples"
)
saveRDS(
  expr_log2_filtered,
  file.path(
    processed_dir,
    "expression_matrix_log2_filtered.rds"
  )
)
# ------------------------------------------------------------
# 7. Expression distribution after log2 transformation
# ------------------------------------------------------------
boxplot_data <- as.data.frame(
  expr_log2_filtered,
  check.names = FALSE
)
boxplot_data$Gene <- rownames(
  boxplot_data
)
boxplot_data <- reshape(
  boxplot_data,
  varying = setdiff(
    colnames(boxplot_data),
    "Gene"
  ),
  v.names = "Expression",
  timevar = "Sample",
  times = setdiff(
    colnames(boxplot_data),
    "Gene"
  ),
  idvar = "Gene",
  direction = "long"
)
boxplot_data$Sample <- factor(
  boxplot_data$Sample,
  levels = c(
    "MEFs",
    "48h",
    "pre-i#1",
    "pre-i#2",
    "ESCs",
    "48h_OSKMFra1",
    "48h_OSKMcJun",
    "48h_OSKMEsrrb",
    "MEFs_Oct4",
    "MEFs_Sox2",
    "MEFs_Klf4"
  )
)
sample_colors <- c(
  "MEFs" = "#F7CAD0",
  "48h" = "#F26CA7",
  "pre-i#1" = "#D63384",
  "pre-i#2" = "#D63384",
  "ESCs" = "#7A003C",
  "48h_OSKMFra1" = "#FFD9E8",
  "48h_OSKMcJun" = "#FFD9E8",
  "48h_OSKMEsrrb" = "#FFD9E8",
  "MEFs_Oct4" = "#F7CAD0",
  "MEFs_Sox2" = "#F7CAD0",
  "MEFs_Klf4" = "#F7CAD0"
)
p_log2 <- ggplot(
  boxplot_data,
  aes(
    x = Sample,
    y = Expression,
    fill = Sample
  )
) +
  geom_boxplot(
    width = 0.7,
    outlier.size = 0.3,
    linewidth = 0.4
  ) +
  scale_fill_manual(
    values = sample_colors
  ) +
  labs(
    title = "Distribution of log2(RPKM + 1) expression values",
    x = NULL,
    y = "log2(RPKM + 1)"
  ) +
  theme_classic(
    base_size = 13
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold"
    ),
    axis.text.y = element_text(
      face = "bold"
    ),
    axis.title.y = element_text(
      face = "bold"
    ),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )
ggsave(
  file.path(
    main_figure_dir,
    "Expression_distribution_log2.png"
  ),
  p_log2,
  width = 9,
  height = 6,
  dpi = 300
)
# ------------------------------------------------------------
# 8. Principal component analysis
# ------------------------------------------------------------
expr_pca <- t(
  expr_log2_filtered
)
pca <- prcomp(
  expr_pca,
  center = TRUE,
  scale. = TRUE
)
variance_explained <- summary(
  pca
)$importance[2, ] * 100
pc1_variance <- round(
  variance_explained[1],
  1
)
pc2_variance <- round(
  variance_explained[2],
  1
)
message(
  "PC1 variance explained: ",
  pc1_variance,
  "%"
)
message(
  "PC2 variance explained: ",
  pc2_variance,
  "%"
)
# ------------------------------------------------------------
# 9. Define sample groups
# ------------------------------------------------------------
sample_groups <- c(
  "MEFs" = "MEFs",
  "48h" = "48 h",
  "pre-i#1" = "pre-iPSC",
  "pre-i#2" = "pre-iPSC",
  "ESCs" = "ESC",
  "48h_OSKMFra1" = "TF overexpression",
  "48h_OSKMcJun" = "TF overexpression",
  "48h_OSKMEsrrb" = "TF overexpression",
  "MEFs_Oct4" = "MEFs",
  "MEFs_Sox2" = "MEFs",
  "MEFs_Klf4" = "MEFs"
)
missing_group_assignments <- setdiff(
  rownames(pca$x),
  names(sample_groups)
)
if (length(missing_group_assignments) > 0) {
  stop(
    "Missing group assignments for: ",
    paste(
      missing_group_assignments,
      collapse = ", "
    )
  )
}
pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Group = unname(
    sample_groups[
      rownames(pca$x)
    ]
  )
)
# ------------------------------------------------------------
# 10. PCA labels
# ------------------------------------------------------------
pca_df$Label <- ""
pca_df$Label[
  pca_df$Sample == "MEFs"
] <- "MEFs"
pca_df$Label[
  pca_df$Sample == "48h"
] <- "48 h"
pca_df$Label[
  pca_df$Sample == "pre-i#1"
] <- "pre-i#1"
pca_df$Label[
  pca_df$Sample == "pre-i#2"
] <- "pre-i#2"
pca_df$Label[
  pca_df$Sample == "ESCs"
] <- "ESCs"
# ------------------------------------------------------------
# 11. PCA visualization
# ------------------------------------------------------------
group_colors <- c(
  "MEFs" = "#F7CAD0",
  "48 h" = "#F26CA7",
  "pre-iPSC" = "#D63384",
  "ESC" = "#7A003C",
  "TF overexpression" = "#FFD9E8"
)
p_pca <- ggplot(
  pca_df,
  aes(
    PC1,
    PC2,
    fill = Group
  )
) +
  geom_point(
    shape = 21,
    size = 5,
    colour = "black",
    stroke = 0.7
  ) +
  geom_text_repel(
    aes(label = Label),
    size = 4.5,
    fontface = "bold",
    colour = "black",
    box.padding = 0.8,
    point.padding = 0.6,
    segment.color = "grey70",
    min.segment.length = 0,
    max.overlaps = Inf
  ) +
  scale_fill_manual(
    values = group_colors
  ) +
  scale_x_continuous(
    expand = expansion(mult = 0.12)
  ) +
  scale_y_continuous(
    expand = expansion(mult = 0.12)
  ) +
  labs(
    title = "Principal Component Analysis",
    x = paste0(
      "PC1 (",
      pc1_variance,
      "% variance explained)"
    ),
    y = paste0(
      "PC2 (",
      pc2_variance,
      "% variance explained)"
    )
  ) +
  coord_cartesian(
    clip = "off"
  ) +
  theme_classic(
    base_size = 15
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    axis.title = element_text(
      face = "bold",
      size = 15
    ),
    axis.text = element_text(
      size = 12,
      colour = "black"
    ),
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(
      size = 12
    )
  )
ggsave(
  file.path(
    main_figure_dir,
    "PCA.png"
  ),
  p_pca,
  width = 9,
  height = 7,
  dpi = 300
)
# ------------------------------------------------------------
# 12. Sample-to-sample Pearson correlation
# ------------------------------------------------------------
cor_matrix <- cor(
  expr_log2_filtered,
  method = "pearson"
)
message(
  "Pearson correlation range: ",
  round(
    min(cor_matrix),
    2
  ),
  "–",
  round(
    max(cor_matrix),
    2
  )
)
annotation_col <- data.frame(
  Group = unname(
    sample_groups[
      colnames(cor_matrix)
    ]
  )
)
rownames(annotation_col) <- colnames(
  cor_matrix
)
annotation_colors <- list(
  Group = group_colors
)
pheatmap(
  cor_matrix,
  color = colorRampPalette(
    c(
      "white",
      "#F7CAD0",
      "#D63384",
      "#7A003C"
    )
  )(100),
  border_color = NA,
  clustering_method = "complete",
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  annotation_colors = annotation_colors,
  display_numbers = TRUE,
  number_format = "%.2f",
  number_color = "black",
  fontsize_number = 10,
  fontsize = 11,
  fontsize_row = 11,
  fontsize_col = 11,
  angle_col = 45,
  main = "Sample-to-sample Pearson correlation",
  filename = file.path(
    main_figure_dir,
    "Sample_correlation.png"
  ),
  width = 9,
  height = 8
)
# ------------------------------------------------------------
# 13. Hierarchical clustering
# ------------------------------------------------------------
sample_dist <- as.dist(
  1 - cor_matrix
)
hc <- hclust(
  sample_dist,
  method = "complete"
)
dend <- as.dendrogram(
  hc
)
publication_labels <- c(
  "Klf4",
  "Oct4",
  "Sox2",
  "MEFs",
  "48 h",
  "Esrrb",
  "Fra1",
  "c-Jun",
  "ESCs",
  "pre-i#1",
  "pre-i#2"
)
if (length(labels(dend)) != length(publication_labels)) {
  stop(
    "The number of dendrogram labels does not match the number of samples."
  )
}
labels(dend) <- publication_labels
# ------------------------------------------------------------
# 14. Hierarchical clustering visualization
# ------------------------------------------------------------
group_colors_dend <- c(
  "#F7CAD0",
  "#F7CAD0",
  "#F7CAD0",
  "#F7CAD0",
  "#F26CA7",
  "#FFD9E8",
  "#FFD9E8",
  "#FFD9E8",
  "#7A003C",
  "#D63384",
  "#D63384"
)
png(
  filename = file.path(
    supplementary_figure_dir,
    "Hierarchical_clustering.png"
  ),
  width = 2400,
  height = 2200,
  res = 300
)
par(
  mar = c(8, 4, 4, 2)
)
plot(
  dend,
  main = "Hierarchical clustering of samples",
  ylab = "Distance",
  xlab = "",
  hang = -1,
  cex = 1,
  leaflab = "perpendicular",
  font = 2,
  axes = TRUE
)
x <- seq_along(
  labels(dend)
)
points(
  x - 0.22,
  rep(
    -0.02,
    length(x)
  ),
  pch = 21,
  bg = group_colors_dend,
  col = "black",
  lwd = 0.4,
  cex = 1.2,
  xpd = TRUE
)
dev.off()
# ------------------------------------------------------------
# 15. Save quality-control objects
# ------------------------------------------------------------
saveRDS(
  cor_matrix,
  file.path(
    processed_dir,
    "sample_correlation_matrix.rds"
  )
)
saveRDS(
  pca,
  file.path(
    processed_dir,
    "PCA_results.rds"
  )
)
saveRDS(
  hc,
  file.path(
    processed_dir,
    "hierarchical_clustering.rds"
  )
)
# ------------------------------------------------------------
# 16. Completion message
# ------------------------------------------------------------
message(
  "Quality control and preprocessing completed successfully."
)
message(
  "Filtered matrix: ",
  nrow(expr_log2_filtered),
  " genes x ",
  ncol(expr_log2_filtered),
  " samples."
)
