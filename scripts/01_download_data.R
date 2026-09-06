# ============================================================
# Project: SPRY2 Reprogramming
# Script: 01_download_data.R
# Purpose: Download, inspect, and prepare the GEO transcriptomic data
# ============================================================
# ------------------------------------------------------------
# 1. Package setup
# ------------------------------------------------------------
required_packages <- c(
  "GEOquery",
  "Biobase",
  "readxl"
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
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install(
    missing_packages,
    ask = FALSE,
    update = FALSE
  )
}
library(GEOquery)
library(Biobase)
library(readxl)
# ------------------------------------------------------------
# 2. Project directories
# ------------------------------------------------------------
raw_dir <- file.path("data", "raw")
processed_dir <- file.path("data", "processed")
geo_dir <- file.path(raw_dir, "GSE90894")
dir.create(
  geo_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  processed_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
# ------------------------------------------------------------
# 3. Retrieve GEO sample metadata
# ------------------------------------------------------------
message("Retrieving GEO metadata for GSE90894...")
gse_rna <- getGEO(
  "GSE90894",
  GSEMatrix = TRUE,
  destdir = raw_dir
)
eset_rna <- gse_rna[[1]]
sample_metadata <- pData(eset_rna)
message(
  "Number of RNA-seq samples: ",
  nrow(sample_metadata)
)
metadata_columns <- intersect(
  c("title", "source_name_ch1"),
  colnames(sample_metadata)
)
if (length(metadata_columns) > 0) {
  print(
    sample_metadata[
      ,
      metadata_columns,
      drop = FALSE
    ]
  )
}
write.csv(
  sample_metadata,
  file.path(
    processed_dir,
    "sample_metadata.csv"
  ),
  row.names = TRUE
)
# ------------------------------------------------------------
# 4. Download supplementary expression data
# ------------------------------------------------------------
expression_file <- file.path(
  geo_dir,
  "GSE90894_RPKM_mRNAseq_table.xlsx"
)
if (!file.exists(expression_file)) {
  message(
    "Expression matrix not found locally."
  )
  message(
    "Downloading supplementary files for GSE90894..."
  )
  getGEOSuppFiles(
    "GSE90894",
    baseDir = raw_dir
  )
}
# ------------------------------------------------------------
# 5. Locate the primary expression matrix
# ------------------------------------------------------------
expression_candidates <- list.files(
  path = geo_dir,
  pattern = "GSE90894_RPKM_mRNAseq_table\\.xlsx$",
  recursive = TRUE,
  full.names = TRUE
)
if (length(expression_candidates) == 0) {
  stop(
    paste0(
      "The primary expression matrix could not be found.\n",
      "Expected file: ",
      expression_file
    )
  )
}
expression_file <- expression_candidates[1]
message(
  "Primary expression matrix found:\n",
  expression_file
)
# ------------------------------------------------------------
# 6. Import expression matrix
# ------------------------------------------------------------
expression_matrix <- read_excel(
  expression_file,
  skip = 3
)
message(
  "Imported expression matrix: ",
  nrow(expression_matrix),
  " rows x ",
  ncol(expression_matrix),
  " columns"
)
# ------------------------------------------------------------
# 7. Validate gene identifiers
# ------------------------------------------------------------
if (!"name" %in% colnames(expression_matrix)) {
  stop(
    "The expression matrix does not contain the expected 'name' column."
  )
}
gene_ids <- expression_matrix$name
if (anyNA(gene_ids) || any(gene_ids == "")) {
  stop(
    "Missing or empty gene identifiers were detected."
  )
}
# ------------------------------------------------------------
# 8. Construct expression matrix
# ------------------------------------------------------------
expression_matrix <- as.data.frame(
  expression_matrix,
  check.names = FALSE
)
expression_matrix$name <- NULL
rownames(expression_matrix) <- make.unique(gene_ids)
# ------------------------------------------------------------
# 9. Report duplicated gene identifiers
# ------------------------------------------------------------
duplicate_count <- sum(duplicated(gene_ids))
message(
  "Duplicated gene identifiers detected: ",
  duplicate_count
)
if (duplicate_count > 0) {
  message(
    "Duplicated identifiers were resolved using make.unique()."
  )
}
if (anyDuplicated(rownames(expression_matrix)) > 0) {
  stop(
    "Duplicate row names remain after identifier processing."
  )
}
# ------------------------------------------------------------
# 10. Validate expression values
# ------------------------------------------------------------
if (!all(vapply(
  expression_matrix,
  is.numeric,
  logical(1)
))) {
  stop(
    "One or more expression columns are not numeric."
  )
}
if (anyNA(expression_matrix)) {
  warning(
    "Missing values were detected in the expression matrix."
  )
}
message(
  "Final expression matrix: ",
  nrow(expression_matrix),
  " genes x ",
  ncol(expression_matrix),
  " samples"
)
# ------------------------------------------------------------
# 11. Save cleaned expression matrix
# ------------------------------------------------------------
saveRDS(
  expression_matrix,
  file.path(
    processed_dir,
    "expression_matrix_clean.rds"
  )
)
write.csv(
  expression_matrix,
  file.path(
    processed_dir,
    "expression_matrix_clean.csv"
  ),
  row.names = TRUE
)
message(
  "Cleaned expression matrix saved to data/processed/"
)
# ------------------------------------------------------------
# 12. Session information
# ------------------------------------------------------------
message(
  "Data import and preparation completed successfully."
)
print(sessionInfo())
