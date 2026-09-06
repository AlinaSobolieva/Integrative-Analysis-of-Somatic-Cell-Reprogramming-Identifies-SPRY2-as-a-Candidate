# ============================================================
# Project: SPRY2 Reprogramming
# Script: 03_differential_expression.R
# Purpose: Differential expression analysis across the
#          somatic cell reprogramming trajectory
# ============================================================
# ------------------------------------------------------------
# 1. Package setup
# ------------------------------------------------------------
required_packages <- c(
  "limma",
  "ggplot2",
  "ggrepel",
  "patchwork"
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
  bioc_packages <- intersect(
    missing_packages,
    c("limma")
  )
  cran_packages <- setdiff(
    missing_packages,
    bioc_packages
  )
  if (length(cran_packages) > 0) {
    install.packages(cran_packages)
  }
  if (length(bioc_packages) > 0) {
    if (!requireNamespace(
      "BiocManager",
      quietly = TRUE
    )) {
      install.packages("BiocManager")
    }
    BiocManager::install(
      bioc_packages,
      ask = FALSE,
      update = FALSE
    )
  }
}
library(limma)
library(ggplot2)
library(ggrepel)
library(patchwork)
# ------------------------------------------------------------
# 2. Project directories
# ------------------------------------------------------------
processed_dir <- file.path(
  "data",
  "processed"
)
de_dir <- file.path(
  "results",
  "differential_expression"
)
figure_dir <- file.path(
  "results",
  "figures",
  "main"
)
dir.create(
  de_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
# ------------------------------------------------------------
# 3. Load filtered expression matrix
# ------------------------------------------------------------
expression_file <- file.path(
  processed_dir,
  "expression_matrix_log2_filtered.rds"
)
if (!file.exists(expression_file)) {
  stop(
    "Filtered expression matrix not found: ",
    expression_file,
    "\nRun 02_quality_control.R first."
  )
}
expr_log2_filtered <- readRDS(
  expression_file
)
expr_log2_filtered <- as.matrix(
  expr_log2_filtered
)
message(
  "Input matrix: ",
  nrow(expr_log2_filtered),
  " genes x ",
  ncol(expr_log2_filtered),
  " samples"
)
# ------------------------------------------------------------
# 4. Select samples for the main trajectory
# ------------------------------------------------------------
main_samples <- c(
  "MEFs",
  "48h",
  "pre-i#1",
  "pre-i#2",
  "ESCs"
)
missing_samples <- setdiff(
  main_samples,
  colnames(expr_log2_filtered)
)
if (length(missing_samples) > 0) {
  stop(
    "The following required samples are missing: ",
    paste(
      missing_samples,
      collapse = ", "
    )
  )
}
expr_main <- expr_log2_filtered[
  ,
  main_samples,
  drop = FALSE
]
message(
  "Main trajectory matrix: ",
  nrow(expr_main),
  " genes x ",
  ncol(expr_main),
  " samples"
)
# ------------------------------------------------------------
# 5. Define experimental groups
# ------------------------------------------------------------
# H48 is used as a syntactically valid model name
# for the biological 48 h sample.
group_main <- factor(
  c(
    "MEFs",
    "H48",
    "pre_iPSC",
    "pre_iPSC",
    "ESC"
  ),
  levels = c(
    "MEFs",
    "H48",
    "pre_iPSC",
    "ESC"
  )
)
names(group_main) <- colnames(
  expr_main
)
print(
  data.frame(
    Sample = names(group_main),
    Group = group_main
  )
)
# ------------------------------------------------------------
# 6. Design matrix
# ------------------------------------------------------------
design <- model.matrix(
  ~ 0 + group_main
)
colnames(design) <- levels(
  group_main
)
message("Design matrix:")
print(design)
# ------------------------------------------------------------
# 7. Define biological contrasts
# ------------------------------------------------------------
contrast_matrix <- makeContrasts(
  Early_reprogramming =
    H48 - MEFs,
  Intermediate_reprogramming =
    pre_iPSC - H48,
  Late_reprogramming =
    ESC - pre_iPSC,
  ESC_vs_MEF =
    ESC - MEFs,
  levels = design
)
message("Contrast matrix:")
print(contrast_matrix)
# ------------------------------------------------------------
# 8. Fit linear model
# ------------------------------------------------------------
fit <- lmFit(
  expr_main,
  design
)
fit_contrasts <- contrasts.fit(
  fit,
  contrast_matrix
)
fit_ebayes <- eBayes(
  fit_contrasts
)
message(
  "Linear modeling and empirical Bayes moderation completed."
)
# ------------------------------------------------------------
# 9. Extract differential expression results
# ------------------------------------------------------------
extract_de_results <- function(
  fit_object,
  coefficient
) {
  result <- topTable(
    fit_object,
    coef = coefficient,
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
  )
  result$Gene <- rownames(
    result
  )
  result <- result[
    ,
    c(
      "Gene",
      "logFC",
      "AveExpr",
      "t",
      "P.Value",
      "adj.P.Val",
      "B"
    )
  ]
  # Strict statistical significance
  result$Significance <- "Not significant"
  result$Significance[
    result$adj.P.Val < 0.05 &
      abs(result$logFC) >= 1
  ] <- "Significant DEG"
  # Exploratory candidates
  result$Significance[
    result$P.Value < 0.05 &
      abs(result$logFC) >= 1 &
      result$adj.P.Val >= 0.05
  ] <- "Exploratory candidate"
  # Volcano plot categories
  result$Category <- "Not significant"
  result$Category[
    result$P.Value < 0.05 &
      result$logFC >= 1
  ] <- "Exploratory up"
  result$Category[
    result$P.Value < 0.05 &
      result$logFC <= -1
  ] <- "Exploratory down"
  result
}
de_results <- list(
  Early_reprogramming =
    extract_de_results(
      fit_ebayes,
      "Early_reprogramming"
    ),
  Intermediate_reprogramming =
    extract_de_results(
      fit_ebayes,
      "Intermediate_reprogramming"
    ),
  Late_reprogramming =
    extract_de_results(
      fit_ebayes,
      "Late_reprogramming"
    ),
  ESC_vs_MEF =
    extract_de_results(
      fit_ebayes,
      "ESC_vs_MEF"
    )
)
# ------------------------------------------------------------
# 10. Summarize differential expression results
# ------------------------------------------------------------
de_summary <- do.call(
  rbind,
  lapply(
    names(de_results),
    function(contrast_name) {
      result <- de_results[[contrast_name]]
      data.frame(
        Contrast = contrast_name,
        Significant_DEGs = sum(
          result$Significance ==
            "Significant DEG"
        ),
        Exploratory_candidates = sum(
          result$Significance ==
            "Exploratory candidate"
        ),
        Upregulated_exploratory = sum(
          result$Category ==
            "Exploratory up"
        ),
        Downregulated_exploratory = sum(
          result$Category ==
            "Exploratory down"
        )
      )
    }
  )
)
print(de_summary)
write.csv(
  de_summary,
  file.path(
    de_dir,
    "DE_summary.csv"
  ),
  row.names = FALSE
)
# ------------------------------------------------------------
# 11. Save complete differential expression tables
# ------------------------------------------------------------
for (contrast_name in names(de_results)) {
  write.csv(
    de_results[[contrast_name]],
    file.path(
      de_dir,
      paste0(
        contrast_name,
        "_DE_results.csv"
      )
    ),
    row.names = FALSE
  )
}
saveRDS(
  de_results,
  file.path(
    de_dir,
    "differential_expression_results.rds"
  )
)
# ------------------------------------------------------------
# 12. SPRY2 differential expression summary
# ------------------------------------------------------------
if (!"Spry2" %in% rownames(expr_main)) {
  warning(
    "Spry2 was not found in the expression matrix."
  )
} else {
  spry2_summary <- do.call(
    rbind,
    lapply(
      names(de_results),
      function(contrast_name) {
        result <- de_results[[contrast_name]]
        row <- result[
          result$Gene == "Spry2",
          ,
          drop = FALSE
        ]
        data.frame(
          Contrast = contrast_name,
          logFC = row$logFC,
          P.Value = row$P.Value,
          FDR = row$adj.P.Val
        )
      }
    )
  )
  print(spry2_summary)
  write.csv(
    spry2_summary,
    file.path(
      de_dir,
      "SPRY2_DE_summary.csv"
    ),
    row.names = FALSE
  )
}
# ------------------------------------------------------------
# 13. Volcano plot function
# ------------------------------------------------------------
create_volcano <- function(
  result,
  title,
  labels_to_show
) {
  result$Label <- ifelse(
    result$Gene %in% labels_to_show,
    result$Gene,
    NA
  )
  ggplot(
    result,
    aes(
      x = logFC,
      y = -log10(P.Value)
    )
  ) +
    geom_point(
      aes(
        colour = Category
      ),
      alpha = 0.55,
      size = 1.5
    ) +
    geom_vline(
      xintercept = c(
        -1,
        1
      ),
      linetype = "dashed",
      colour = "grey55",
      linewidth = 0.5
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed",
      colour = "grey55",
      linewidth = 0.5
    ) +
    geom_point(
      data = subset(
        result,
        Gene == "Spry2"
      ),
      colour = "#C2185B",
      fill = "white",
      size = 4.5,
      shape = 21,
      stroke = 1.5
    ) +
    geom_text_repel(
      data = subset(
        result,
        !is.na(Label)
      ),
      aes(
        label = Label
      ),
      size = 3.8,
      fontface = "bold",
      colour = "grey15",
      box.padding = 0.7,
      point.padding = 0.4,
      min.segment.length = 0,
      max.overlaps = Inf,
      seed = 42
    ) +
    geom_text_repel(
      data = subset(
        result,
        Gene == "Spry2"
      ),
      aes(
        label = Gene
      ),
      size = 4.2,
      fontface = "bold",
      colour = "grey15",
      nudge_x = 0.55,
      nudge_y = -0.25,
      box.padding = 0.8,
      point.padding = 0.5,
      min.segment.length = 0,
      segment.color = "#C2185B",
      segment.size = 0.5,
      seed = 42
    ) +
    scale_colour_manual(
      values = c(
        "Exploratory down" = "#6A5ACD",
        "Exploratory up" = "#C2185B",
        "Not significant" = "grey80"
      )
    ) +
    labs(
      title = title,
      subtitle =
        "Exploratory candidates: nominal P < 0.05 and |log₂FC| ≥ 1",
      x = expression(
        log[2] ~ fold ~ change
      ),
      y = expression(
        -log[10](P)
      ),
      colour = NULL
    ) +
    theme_classic(
      base_size = 13
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 16
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        size = 9,
        colour = "grey35"
      ),
      axis.title = element_text(
        face = "bold",
        size = 12
      ),
      axis.text = element_text(
        colour = "black",
        size = 10
      ),
      legend.position = "top",
      legend.text = element_text(
        size = 9
      ),
      plot.margin = margin(
        10,
        15,
        10,
        10
      )
    )
}
# ------------------------------------------------------------
# 14. Define genes for annotation
# ------------------------------------------------------------
genes_to_label <- list(
  Early_reprogramming = c(
    "Pou5f1",
    "Sostdc1",
    "Myh11",
    "Meox1"
  ),
  Intermediate_reprogramming = c(
    "Sfrp2",
    "Fmod",
    "Col5a2",
    "Cxcl12",
    "Pou5f1",
    "Nov"
  ),
  Late_reprogramming = c(
    "Esrrb",
    "Lin28a",
    "Gm13051",
    "Dppa3",
    "Rbmxl2"
  ),
  ESC_vs_MEF = c(
    "Pou5f1",
    "Esrrb",
    "Lin28a",
    "Epcam"
  )
)
# ------------------------------------------------------------
# 15. Generate volcano plots
# ------------------------------------------------------------
volcano_Early <- create_volcano(
  de_results$Early_reprogramming,
  "48 h vs MEFs",
  genes_to_label$Early_reprogramming
)
volcano_Intermediate <- create_volcano(
  de_results$Intermediate_reprogramming,
  "pre-iPSC vs 48 h",
  genes_to_label$Intermediate_reprogramming
)
volcano_Late <- create_volcano(
  de_results$Late_reprogramming,
  "ESC vs pre-iPSC",
  genes_to_label$Late_reprogramming
)
volcano_ESC_vs_MEF <- create_volcano(
  de_results$ESC_vs_MEF,
  "ESC vs MEFs",
  genes_to_label$ESC_vs_MEF
)
# ------------------------------------------------------------
# 16. Save individual volcano plots
# ------------------------------------------------------------
ggsave(
  file.path(
    figure_dir,
    "Volcano_48h_vs_MEF.png"
  ),
  volcano_Early,
  width = 8,
  height = 6,
  dpi = 300
)
ggsave(
  file.path(
    figure_dir,
    "Volcano_pre-iPSC_vs_48h.png"
  ),
  volcano_Intermediate,
  width = 8,
  height = 6,
  dpi = 300
)
ggsave(
  file.path(
    figure_dir,
    "Volcano_ESC_vs_pre-iPSC.png"
  ),
  volcano_Late,
  width = 8,
  height = 6,
  dpi = 300
)
ggsave(
  file.path(
    figure_dir,
    "Volcano_ESC_vs_MEF.png"
  ),
  volcano_ESC_vs_MEF,
  width = 8,
  height = 6,
  dpi = 300
)
# ------------------------------------------------------------
# 17. Create main Figure 2
# ------------------------------------------------------------
figure_2 <- (
  volcano_Early +
  volcano_Intermediate
) /
(
  volcano_Late +
  volcano_ESC_vs_MEF
) +
  plot_annotation(
    title =
      "Differential expression across the reprogramming trajectory"
  )
ggsave(
  file.path(
    figure_dir,
    "Figure_2_DE.png"
  ),
  figure_2,
  width = 14,
  height = 11,
  dpi = 300
)
# ------------------------------------------------------------
# 18. Completion message
# ------------------------------------------------------------
message(
  "Differential expression analysis completed successfully."
)
message(
  "Significant DEGs (FDR < 0.05, |log2FC| >= 1): ",
  sum(
    de_summary$Significant_DEGs
  )
)
message(
  "Exploratory candidates across contrasts: ",
  paste(
    de_summary$Exploratory_candidates,
    collapse = ", "
  )
)
