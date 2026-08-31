# ==============================================================================
# Stage 1 -> Stage 2 automatic hand-off: the candidate-drug list.
#
# Both CKM_Drug (Stage 1) pipelines produce a per-drug, per-outcome Cox
# screen (core/R/cox_screen.R) and flag some (ATC code, outcome) pairs as
# `candidate == TRUE`. The side-effect model (Stage 2) restricts its drug
# predictors to that candidate set before Boruta ever runs (a clinically/
# statistically motivated pool, not "every drug ever prescribed").
#
# Historically this hand-off was a MANUAL file copy: Korea's
# candidate_drugs.csv and Japan's cox_results.csv were each hand-copied from
# the sibling CKM_Drug project's results/<run_id>/ output into the
# SideEffect_Model project's own folder. This file replaces that with an
# automatic loader for BOTH countries: read the CKM_DRUG database's
# cox_results table (always current, written by Stage 1's run_pipeline.R)
# when reachable, else fall back to Stage 1's own most recent
# results/<run_id>/cox_results.rds file checkpoint (see
# core/R/db_common.R's read_table_or_checkpoint()/find_latest_checkpoint()).
#
# Column names are already standardized (atc_code/outcome/candidate) because
# core/R/cox_screen.R emits them that way for both countries -- unlike
# before the 2026-08-31 switch (see each country's stage1 CLAUDE-equivalent
# history), there is no more who_atc_code/ATC_CODE naming split to bridge
# here, and no more mixed-ATC-level allow-list to match by LIKE-prefix:
# both countries' Cox screen already operates at a single, exact ATC-code
# granularity (Drug_master.who_atc_code for Japan, NHIS_ATC_MAPPED.xlsx's
# ATC_CODE for Korea), so this file uses a plain exact match throughout.
# ==============================================================================

#' Load the candidate-drug table for one country's Stage 1 pipeline: DB
#' table `table_name` in CKM_DRUG (via `con`, if reachable) or else the
#' latest results/<run_id>/cox_results.rds under `stage1_results_root`.
load_candidate_drugs <- function(con, stage1_results_root, table_name = "cox_results", log = message) {
  df <- read_table_or_checkpoint(
    con, table_name, stage1_results_root,
    checkpoint_name = "cox_results", log = log
  )

  required_cols <- c("atc_code", "outcome", "candidate")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Candidate-drug source '", table_name, "' is missing expected column(s): ",
         paste(missing_cols, collapse = ", "))
  }

  df$atc_code <- toupper(trimws(df$atc_code))
  df$candidate <- as.logical(df$candidate)

  bad <- df$atc_code[!is.na(df$atc_code) & !grepl("^[A-Z0-9]{3,7}$", df$atc_code)]
  if (length(bad) > 0) {
    stop("Candidate-drug source contains malformed ATC code(s): ", paste(unique(bad), collapse = ", "))
  }

  log(sprintf(
    "Loaded %d candidate-drug row(s) (%d flagged candidate == TRUE) from '%s'.",
    nrow(df), sum(df$candidate, na.rm = TRUE), table_name
  ))
  df
}

#' The union of ATC codes flagged candidate == TRUE for ANY outcome -- used
#' to build the SQL-side allow-list filter in each country's Stage 2
#' 02_extract_drug_exposure.R (same role as the old load_atc_allowlist()).
candidate_atc_union <- function(candidate_drugs) {
  codes <- candidate_drugs$atc_code[candidate_drugs$candidate == TRUE]
  unique(codes[!is.na(codes) & nzchar(codes)])
}

#' ATC codes flagged candidate == TRUE for a specific Stage-2 outcome (which
#' is a composite of one or more of Stage 1's four outcomes -- see
#' cox_outcome_map in each country's stage2 config). Used by
#' core/R/boruta_selection.R to restrict each outcome's own Boruta run to
#' its own eligible drug columns.
candidate_atc_for_outcome <- function(candidate_drugs, stage1_outcomes) {
  codes <- candidate_drugs$atc_code[
    candidate_drugs$candidate == TRUE & candidate_drugs$outcome %in% stage1_outcomes
  ]
  unique(codes[!is.na(codes) & nzchar(codes)])
}

#' SQL WHERE-clause fragment: `column` IN (codes). Exact match -- see file
#' header for why prefix/LIKE matching is no longer needed for either
#' country.
atc_allowlist_sql_in <- function(codes, column) {
  sprintf("%s IN (%s)", column, paste(sprintf("'%s'", codes), collapse = ", "))
}
