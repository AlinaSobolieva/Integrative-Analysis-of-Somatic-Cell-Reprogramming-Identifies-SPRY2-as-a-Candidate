# ============================================================
# Project: SPRY2 Reprogramming
# Script: 04_functional_enrichment.R
# Purpose: Functional enrichment analysis across the
#          somatic cell reprogramming trajectory
# ============================================================
# ------------------------------------------------------------
# 1. Package setup
# ------------------------------------------------------------
required_packages <- c(
  "clusterProfiler",
  "org.Mm.eg.db",
  "AnnotationDbi",
  "enrichplot",
  "ReactomePA",
  "msigdbr",
  "ggplot2",
  "patchwork",
  "stringr"
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
    c(
      "clusterProfiler",
      "org.Mm.eg.db",
      "AnnotationDbi",
      "enrichplot",
      "ReactomePA"
    )
  )
  cran_packages <- setdiff(
    missing_packages,
    bioc_packages
  )
  if (length(cran_packages) > 0) {
    install.packages(
      cran_packages
    )
  }
  if (length(bioc_packages) > 0) {
    if (!requireNamespace(
      "BiocManager",
      quietly = TRUE
    )) {
      install.packages(
        "BiocManager"
      )
    }
    BiocManager::install(
      bioc_packages,
      ask = FALSE,
      update = FALSE
    )
  }
}
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ReactomePA)
library(msigdbr)
library(ggplot2)
library(patchwork)
library(stringr)
# ------------------------------------------------------------
# 2. Project directories
# ------------------------------------------------------------
de_dir <- file.path(
  "results",
  "differential_expression"
)
enrichment_dir <- file.path(
  "results",
  "functional_enrichment"
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
  enrichment_dir,
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
# 3. Load differential expression results
# ------------------------------------------------------------
de_file <- file.path(
  de_dir,
  "differential_expression_results.rds"
)
if (!file.exists(de_file)) {
  stop(
    "Differential expression results not found: ",
    de_file,
    "\nRun 03_differential_expression.R first."
  )
}
de_results <- readRDS(
  de_file
)
required_contrasts <- c(
  "Early_reprogramming",
  "Intermediate_reprogramming",
  "Late_reprogramming",
  "ESC_vs_MEF"
)
missing_contrasts <- setdiff(
  required_contrasts,
  names(de_results)
)
if (length(missing_contrasts) > 0) {
  stop(
    "Missing differential expression contrasts: ",
    paste(
      missing_contrasts,
      collapse = ", "
    )
  )
}
# ------------------------------------------------------------
# 4. Define exploratory gene sets
# ------------------------------------------------------------
get_gene_set <- function(
  de_table,
  direction
) {
  if (direction == "up") {
    genes <- de_table$Gene[
      de_table$P.Value < 0.05 &
        de_table$logFC >= 1
    ]
  } else {
    genes <- de_table$Gene[
      de_table$P.Value < 0.05 &
        de_table$logFC <= -1
    ]
  }
  unique(
    genes[
      !is.na(genes) &
        genes != ""
    ]
  )
}
gene_sets <- list(
  Early_up = get_gene_set(
    de_results$Early_reprogramming,
    "up"
  ),
  Early_down = get_gene_set(
    de_results$Early_reprogramming,
    "down"
  ),
  Intermediate_up = get_gene_set(
    de_results$Intermediate_reprogramming,
    "up"
  ),
  Intermediate_down = get_gene_set(
    de_results$Intermediate_reprogramming,
    "down"
  ),
  Late_up = get_gene_set(
    de_results$Late_reprogramming,
    "up"
  ),
  Late_down = get_gene_set(
    de_results$Late_reprogramming,
    "down"
  ),
  ESC_MEF_up = get_gene_set(
    de_results$ESC_vs_MEF,
    "up"
  ),
  ESC_MEF_down = get_gene_set(
    de_results$ESC_vs_MEF,
    "down"
  )
)
message("Exploratory gene-set sizes:")
print(
  sapply(
    gene_sets,
    length
  )
)
# ------------------------------------------------------------
# 5. Define gene universe
# ------------------------------------------------------------
gene_universe <- unique(
  unlist(
    lapply(
      de_results,
      function(x) x$Gene
    )
  )
)
gene_universe <- gene_universe[
  !is.na(gene_universe) &
    gene_universe != ""
]
message(
  "Gene universe: ",
  length(gene_universe),
  " genes"
)
# ------------------------------------------------------------
# 6. Convert mouse SYMBOLs to ENTREZ IDs
# ------------------------------------------------------------
symbol_to_entrez <- function(
  genes
) {
  genes <- unique(
    genes[
      !is.na(genes) &
        genes != ""
    ]
  )
  if (length(genes) == 0) {
    return(
      character(0)
    )
  }
  annotation <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys = genes,
    keytype = "SYMBOL",
    columns = c(
      "SYMBOL",
      "ENTREZID"
    )
  )
  annotation <- annotation[
    !is.na(annotation$ENTREZID),
  ]
  unique(
    annotation$ENTREZID
  )
}
gene_sets_entrez <- lapply(
  gene_sets,
  symbol_to_entrez
)
gene_universe_entrez <- symbol_to_entrez(
  gene_universe
)
message(
  "ENTREZ gene universe: ",
  length(gene_universe_entrez),
  " genes"
)
# ------------------------------------------------------------
# 7. GO Biological Process enrichment
# ------------------------------------------------------------
run_go <- function(
  genes
) {
  if (length(genes) == 0) {
    return(NULL)
  }
  enrichGO(
    gene = genes,
    universe = gene_universe,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
}
go_results <- lapply(
  gene_sets,
  run_go
)
# ------------------------------------------------------------
# 8. Simplify redundant GO terms
# ------------------------------------------------------------
simplify_go <- function(
  enrichment_result
) {
  if (is.null(enrichment_result)) {
    return(NULL)
  }
  result_df <- as.data.frame(
    enrichment_result
  )
  if (nrow(result_df) == 0) {
    return(
      enrichment_result
    )
  }
  simplify(
    enrichment_result,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
  )
}
go_simplified <- lapply(
  go_results,
  simplify_go
)
# ------------------------------------------------------------
# 9. Save GO results
# ------------------------------------------------------------
for (set_name in names(go_simplified)) {
  result <- go_simplified[[set_name]]
  if (!is.null(result)) {
    write.csv(
      as.data.frame(result),
      file.path(
        enrichment_dir,
        paste0(
          "GO_BP_",
          set_name,
          ".csv"
        )
      ),
      row.names = FALSE
    )
    saveRDS(
      result,
      file.path(
        enrichment_dir,
        paste0(
          "GO_BP_",
          set_name,
          ".rds"
        )
      )
    )
  }
}
# ------------------------------------------------------------
# 10. KEGG enrichment
# ------------------------------------------------------------
run_kegg <- function(
  genes
) {
  if (length(genes) == 0) {
    return(NULL)
  }
  enrichKEGG(
    gene = genes,
    universe = gene_universe_entrez,
    organism = "mmu",
    keyType = "ncbi-geneid",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
}
kegg_results <- lapply(
  gene_sets_entrez,
  run_kegg
)
# ------------------------------------------------------------
# 11. Save KEGG results
# ------------------------------------------------------------
for (set_name in names(kegg_results)) {
  result <- kegg_results[[set_name]]
  if (!is.null(result)) {
    write.csv(
      as.data.frame(result),
      file.path(
        enrichment_dir,
        paste0(
          "KEGG_",
          set_name,
          ".csv"
        )
      ),
      row.names = FALSE
    )
    saveRDS(
      result,
      file.path(
        enrichment_dir,
        paste0(
          "KEGG_",
          set_name,
          ".rds"
        )
      )
    )
  }
}
# ------------------------------------------------------------
# 12. Reactome pathway enrichment
# ------------------------------------------------------------
run_reactome <- function(
  genes
) {
  if (length(genes) == 0) {
    return(NULL)
  }
  enrichPathway(
    gene = genes,
    universe = gene_universe_entrez,
    organism = "mouse",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
}
reactome_results <- lapply(
  gene_sets_entrez,
  run_reactome
)
# ------------------------------------------------------------
# 13. Save Reactome results
# ------------------------------------------------------------
for (set_name in names(reactome_results)) {
  result <- reactome_results[[set_name]]
  if (!is.null(result)) {
    write.csv(
      as.data.frame(result),
      file.path(
        enrichment_dir,
        paste0(
          "Reactome_",
          set_name,
          ".csv"
        )
      ),
      row.names = FALSE
    )
    saveRDS(
      result,
      file.path(
        enrichment_dir,
        paste0(
          "Reactome_",
          set_name,
          ".rds"
        )
      )
    )
  }
}
# ------------------------------------------------------------
# 14. Build ranked gene lists for GSEA
# ------------------------------------------------------------
make_ranked_list <- function(
  de_table
) {
  rank <- de_table$t
  names(rank) <- de_table$Gene
  rank <- rank[
    !is.na(rank) &
      !is.na(names(rank)) &
      names(rank) != ""
  ]
  rank <- rank[
    !duplicated(names(rank))
  ]
  sort(
    rank,
    decreasing = TRUE
  )
}
ranked_lists <- list(
  Early =
    make_ranked_list(
      de_results$Early_reprogramming
    ),
  Intermediate =
    make_ranked_list(
      de_results$Intermediate_reprogramming
    ),
  Late =
    make_ranked_list(
      de_results$Late_reprogramming
    ),
  ESC_MEF =
    make_ranked_list(
      de_results$ESC_vs_MEF
    )
)
message("Ranked gene-list sizes:")
print(
  sapply(
    ranked_lists,
    length
  )
)
# ------------------------------------------------------------
# 15. Load MSigDB Hallmark gene sets
# ------------------------------------------------------------
msigdb_mouse <- msigdbr(
  species = "Mus musculus",
  collection = "H"
)
hallmark_term2gene <- msigdb_mouse[
  ,
  c(
    "gs_name",
    "gene_symbol"
  )
]
hallmark_term2gene <- unique(
  hallmark_term2gene
)
message(
  "Hallmark gene sets loaded: ",
  length(
    unique(
      hallmark_term2gene$gs_name
    )
  )
)
# ------------------------------------------------------------
# 16. GSEA
# ------------------------------------------------------------
run_gsea <- function(
  ranked_list
) {
  GSEA(
    geneList = ranked_list,
    TERM2GENE = hallmark_term2gene,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    minGSSize = 10,
    maxGSSize = 500,
    verbose = FALSE
  )
}
gsea_results <- lapply(
  ranked_lists,
  run_gsea
)
# ------------------------------------------------------------
# 17. Save GSEA results
# ------------------------------------------------------------
for (set_name in names(gsea_results)) {
  result <- gsea_results[[set_name]]
  if (!is.null(result)) {
    write.csv(
      as.data.frame(result),
      file.path(
        enrichment_dir,
        paste0(
          "GSEA_Hallmark_",
          set_name,
          ".csv"
        )
      ),
      row.names = FALSE
    )
    saveRDS(
      result,
      file.path(
        enrichment_dir,
        paste0(
          "GSEA_Hallmark_",
          set_name,
          ".rds"
        )
      )
    )
  }
}
# ------------------------------------------------------------
# 18. Enrichment summary
# ------------------------------------------------------------
count_enriched_terms <- function(
  result
) {
  if (is.null(result)) {
    return(0)
  }
  nrow(
    as.data.frame(result)
  )
}
enrichment_summary <- data.frame(
  Gene_set = names(gene_sets),
  Gene_count = sapply(
    gene_sets,
    length
  ),
  GO_BP_terms = sapply(
    go_simplified,
    count_enriched_terms
  ),
  KEGG_terms = sapply(
    kegg_results,
    count_enriched_terms
  ),
  Reactome_terms = sapply(
    reactome_results,
    count_enriched_terms
  )
)
print(
  enrichment_summary
)
write.csv(
  enrichment_summary,
  file.path(
    enrichment_dir,
    "enrichment_summary.csv"
  ),
  row.names = FALSE
)
# ------------------------------------------------------------
# 19. GO plotting function
# ------------------------------------------------------------
create_go_plot <- function(
  enrichment_result,
  title,
  show_category = 15
) {
  if (is.null(enrichment_result)) {
    return(
      ggplot() +
        theme_void() +
        ggtitle(
          paste0(
            title,
            "\nNo significant GO terms"
          )
        ) +
        theme(
          plot.title = element_text(
            hjust = 0.5,
            face = "bold"
          )
        )
    )
  }
  result_df <- as.data.frame(
    enrichment_result
  )
  if (nrow(result_df) == 0) {
    return(
      ggplot() +
        theme_void() +
        ggtitle(
          paste0(
            title,
            "\nNo significant GO terms"
          )
        ) +
        theme(
          plot.title = element_text(
            hjust = 0.5,
            face = "bold"
          )
        )
    )
  }
  dotplot(
    enrichment_result,
    showCategory = show_category,
    font.size = 9
  ) +
    ggtitle(title) +
    labs(
      x = "Gene Ratio",
      colour = "Adjusted P-value",
      size = "Gene Count"
    ) +
    scale_y_discrete(
      labels = function(x) {
        stringr::str_wrap(
          x,
          width = 42
        )
      }
    ) +
    theme_classic(
      base_size = 12
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 14
      ),
      axis.title.y = element_blank(),
      axis.title.x = element_text(
        face = "bold",
        size = 11
      ),
      axis.text.y = element_text(
        size = 8,
        colour = "black"
      ),
      axis.text.x = element_text(
        size = 9,
        colour = "black"
      ),
      legend.title = element_text(
        face = "bold"
      ),
      panel.grid = element_blank(),
      plot.margin = margin(
        8,
        12,
        8,
        10
      )
    )
}
# ------------------------------------------------------------
# 20. Main Figure 3
# ------------------------------------------------------------
p_intermediate_up <- create_go_plot(
  go_simplified$Intermediate_up,
  "A  Intermediate Upregulated",
  show_category = 15
)
p_late_up <- create_go_plot(
  go_simplified$Late_up,
  "B  Late Upregulated",
  show_category = 15
)
p_esc_mef_down <- create_go_plot(
  go_simplified$ESC_MEF_down,
  "C  ESC vs MEF Downregulated",
  show_category = 15
)
figure_3 <- (
  p_intermediate_up |
  p_late_up |
  p_esc_mef_down
) +
  plot_annotation(
    title =
      "Functional remodeling during somatic cell reprogramming"
  )
ggsave(
  file.path(
    main_figure_dir,
    "Figure_3_Enrichment.png"
  ),
  figure_3,
  width = 18,
  height = 7,
  dpi = 300
)
# ------------------------------------------------------------
# 21. Supplementary Figure 2
# ------------------------------------------------------------
p_early_down <- create_go_plot(
  go_simplified$Early_down,
  "A  Early Downregulated",
  show_category = 15
)
p_intermediate_down <- create_go_plot(
  go_simplified$Intermediate_down,
  "B  Intermediate Downregulated",
  show_category = 15
)
p_late_down <- create_go_plot(
  go_simplified$Late_down,
  "C  Late Downregulated",
  show_category = 15
)
p_esc_mef_up <- create_go_plot(
  go_simplified$ESC_MEF_up,
  "D  ESC vs MEF Upregulated",
  show_category = 15
)
supplementary_figure_2 <- (
  p_early_down +
  p_intermediate_down
) /
(
  p_late_down +
  p_esc_mef_up
) +
  plot_annotation(
    title =
      "Additional functional enrichment analyses"
  )
ggsave(
  file.path(
    supplementary_figure_dir,
    "Supplementary_Figure_S2.png"
  ),
  supplementary_figure_2,
  width = 15,
  height = 12,
  dpi = 300
)
# ------------------------------------------------------------
# 22. Save complete enrichment objects
# ------------------------------------------------------------
saveRDS(
  list(
    GO = go_simplified,
    KEGG = kegg_results,
    Reactome = reactome_results,
    GSEA = gsea_results
  ),
  file.path(
    enrichment_dir,
    "functional_enrichment_results.rds"
  )
)
# ------------------------------------------------------------
# 23. Completion message
# ------------------------------------------------------------
message(
  "Functional enrichment analysis completed successfully."
)
message(
  "GO, KEGG, Reactome, and Hallmark GSEA results were saved to: ",
  enrichment_dir
)
message(
  "Main Figure 3 saved to: ",
  file.path(
    main_figure_dir,
    "Figure_3_Enrichment.png"
  )
)
message(
  "Supplementary Figure 2 saved to: ",
  file.path(
    supplementary_figure_dir,
    "Supplementary_Figure_S2.png"
  )
)
