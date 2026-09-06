# ============================================================
# Project: SPRY2 Reprogramming
# Script: 05_network_analysis.R
# Purpose: Select genes for network analysis and prepare
#          the STRING input gene set
# ============================================================
# ------------------------------------------------------------
# 1. Project directories
# ------------------------------------------------------------
de_dir <- file.path(
  "results",
  "differential_expression"
)
network_dir <- file.path(
  "results",
  "network_analysis"
)
dir.create(
  network_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
# ------------------------------------------------------------
# 2. Load differential expression results
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
deg_Early <- de_results$Early_reprogramming
deg_Intermediate <- de_results$Intermediate_reprogramming
deg_Late <- de_results$Late_reprogramming
deg_ESC_vs_MEF <- de_results$ESC_vs_MEF
# ------------------------------------------------------------
# 3. Validate required columns
# ------------------------------------------------------------
required_columns <- c(
  "Gene",
  "logFC",
  "t",
  "P.Value",
  "adj.P.Val"
)
for (column in required_columns) {
  if (!column %in% colnames(deg_ESC_vs_MEF)) {
    stop(
      "Required column missing from ESC_vs_MEF results: ",
      column
    )
  }
}
# ------------------------------------------------------------
# 4. Inspect leading differential expression results
# ------------------------------------------------------------
message("Top genes from ESC vs MEF:")
print(
  head(
    deg_ESC_vs_MEF,
    20
  )
)
message("Top genes from intermediate reprogramming:")
print(
  head(
    deg_Intermediate,
    20
  )
)
# ------------------------------------------------------------
# 5. Inspect SPRY2 across all contrasts
# ------------------------------------------------------------
spry2_contrasts <- list(
  Early = deg_Early,
  Intermediate = deg_Intermediate,
  Late = deg_Late,
  ESC_vs_MEF = deg_ESC_vs_MEF
)
for (contrast_name in names(spry2_contrasts)) {
  result <- spry2_contrasts[[contrast_name]]
  if ("Spry2" %in% result$Gene) {
    message(
      "\nSPRY2 — ",
      contrast_name
    )
    print(
      result[
        result$Gene == "Spry2",
        ,
        drop = FALSE
      ]
    )
  } else {
    warning(
      "Spry2 was not found in contrast: ",
      contrast_name
    )
  }
}
# ============================================================
# STEP 5.1
# Select top genes for ESC vs MEF network
# ============================================================
# ------------------------------------------------------------
# Top 50 genes increased in ESC
# ------------------------------------------------------------
esc_mef_top_up <- deg_ESC_vs_MEF[
  deg_ESC_vs_MEF$logFC > 0,
  ,
  drop = FALSE
]
esc_mef_top_up <- esc_mef_top_up[
  order(
    esc_mef_top_up$t,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]
esc_mef_top_up <- head(
  esc_mef_top_up,
  50
)
# ------------------------------------------------------------
# Top 50 genes decreased in ESC
# ------------------------------------------------------------
esc_mef_top_down <- deg_ESC_vs_MEF[
  deg_ESC_vs_MEF$logFC < 0,
  ,
  drop = FALSE
]
esc_mef_top_down <- esc_mef_top_down[
  order(
    esc_mef_top_down$t,
    decreasing = FALSE
  ),
  ,
  drop = FALSE
]
esc_mef_top_down <- head(
  esc_mef_top_down,
  50
)
# ------------------------------------------------------------
# Combine selected genes
# ------------------------------------------------------------
network_genes_ESC_MEF <- unique(
  c(
    esc_mef_top_up$Gene,
    esc_mef_top_down$Gene
  )
)
message(
  "Number of selected ESC vs MEF network genes: ",
  length(network_genes_ESC_MEF)
)
# ============================================================
# STEP 5.2
# Inspect selected network genes
# ============================================================
message("\nTop 20 upregulated genes:")
print(
  head(
    esc_mef_top_up[
      ,
      c(
        "Gene",
        "logFC",
        "t",
        "P.Value",
        "adj.P.Val"
      )
    ],
    20
  )
)
message("\nTop 20 downregulated genes:")
print(
  head(
    esc_mef_top_down[
      ,
      c(
        "Gene",
        "logFC",
        "t",
        "P.Value",
        "adj.P.Val"
      )
    ],
    20
  )
)
# ============================================================
# STEP 5.3
# Add SPRY2 as candidate gene of interest
# ============================================================
network_genes_with_SPRY2 <- unique(
  c(
    network_genes_ESC_MEF,
    "Spry2"
  )
)
message(
  "\nNetwork gene set including SPRY2: ",
  length(
    network_genes_with_SPRY2
  ),
  " genes"
)
if (!"Spry2" %in% network_genes_ESC_MEF) {
  message(
    "Spry2 was not among the top 100 ESC vs MEF genes ",
    "and was therefore added separately as the candidate gene."
  )
} else {
  message(
    "Spry2 was already present among the selected ",
    "ESC vs MEF genes."
  )
}
# ============================================================
# STEP 5.4
# SPRY2 expression changes across all stages
# ============================================================
get_spry2_row <- function(
  de_table,
  contrast_name
) {
  spry2_row <- de_table[
    de_table$Gene == "Spry2",
    ,
    drop = FALSE
  ]
  if (nrow(spry2_row) == 0) {
    return(
      data.frame(
        Contrast = contrast_name,
        logFC = NA_real_,
        t = NA_real_,
        P.Value = NA_real_,
        adj.P.Val = NA_real_
      )
    )
  }
  data.frame(
    Contrast = contrast_name,
    logFC = spry2_row$logFC[1],
    t = spry2_row$t[1],
    P.Value = spry2_row$P.Value[1],
    adj.P.Val = spry2_row$adj.P.Val[1]
  )
}
spry2_summary <- rbind(
  get_spry2_row(
    deg_Early,
    "Early_reprogramming"
  ),
  get_spry2_row(
    deg_Intermediate,
    "Intermediate_reprogramming"
  ),
  get_spry2_row(
    deg_Late,
    "Late_reprogramming"
  ),
  get_spry2_row(
    deg_ESC_vs_MEF,
    "ESC_vs_MEF"
  )
)
rownames(
  spry2_summary
) <- NULL
message("\nSPRY2 summary:")
print(
  spry2_summary
)
# ------------------------------------------------------------
# Save SPRY2 summary
# ------------------------------------------------------------
write.csv(
  spry2_summary,
  file.path(
    network_dir,
    "SPRY2_differential_expression_summary.csv"
  ),
  row.names = FALSE
)
# ============================================================
# STEP 5.5
# Prepare gene list for STRING
# ============================================================
string_genes <- network_genes_with_SPRY2
string_genes <- string_genes[
  !is.na(string_genes) &
    string_genes != ""
]
string_genes <- unique(
  string_genes
)
message(
  "\nFinal STRING input gene set: ",
  length(string_genes),
  " genes"
)
message(
  "First selected genes:"
)
print(
  head(
    string_genes,
    20
  )
)
# ============================================================
# STEP 5.6
# Save STRING input list
# ============================================================
string_input_file <- file.path(
  network_dir,
  "STRING_ESC_vs_MEF_genes.txt"
)
write.table(
  string_genes,
  file = string_input_file,
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
message(
  "\nSTRING input list saved to: ",
  string_input_file
)
# ------------------------------------------------------------
# Save complete network-analysis objects
# ------------------------------------------------------------
saveRDS(
  list(
    top_up = esc_mef_top_up,
    top_down = esc_mef_top_down,
    network_genes = network_genes_ESC_MEF,
    network_genes_with_SPRY2 = network_genes_with_SPRY2,
    SPRY2_summary = spry2_summary
  ),
  file.path(
    network_dir,
    "network_analysis_results.rds"
  )
)
# ------------------------------------------------------------
# Completion message
# ------------------------------------------------------------
message(
  "\nNetwork gene selection and STRING input preparation ",
  "completed successfully."
)
