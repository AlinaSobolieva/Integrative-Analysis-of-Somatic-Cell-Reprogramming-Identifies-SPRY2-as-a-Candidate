# ============================================================
# Project: SPRY2 Reprogramming
# Script: 06_network_analysis.R
# Purpose: Locate STRING interaction file and network analysis
# ============================================================
# ============================================================
# 8. Locate STRING interaction file
# ============================================================
# Expected project location
string_file_project <- file.path(
  "results",
  "network_analysis",
  "string_interactions_short.tsv"
)
# Search for any STRING TSV file within the project
string_candidates <- list.files(
  path = ".",
  pattern = "^string_interactions.*\.tsv$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)
# ------------------------------------------------------------
# Select STRING interaction file
# ------------------------------------------------------------
if (file.exists(string_file_project)) {
  string_file <- string_file_project
} else if (length(string_candidates) == 1) {
  string_file <- string_candidates[1]
} else if (length(string_candidates) > 1) {
  message("Multiple STRING interaction files found:")
  print(string_candidates)
  string_file <- string_candidates[1]
} else {
  stop(
    paste0(
      "STRING interaction file was not found in the current project.\n\n",
      "Expected location:\n",
      "  ", string_file_project, "\n\n",
      "Files matching 'string_interactions*.tsv' were not found.\n\n",
      "Current working directory:\n",
      "  ", getwd(), "\n\n",
      "Please place the STRING TSV file in:\n",
      "  results/network_analysis/\n",
      "and rename it to:\n",
      "  string_interactions_short.tsv"
    )
  )
}
message(
  "STRING interaction file found:\n  ",
  normalizePath(
    string_file,
    winslash = "/",
    mustWork = FALSE
  )
)
