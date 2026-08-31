# ==============================================================================
# Demo-mode data loading: lets `Rscript countries/<country>/run_pipeline.R
# --demo` run the WHOLE pipeline (Stage 1's hand-off through all of Stage 2)
# with no live database at all, using the synthetic seed data under
# data/synthetic/<country>/ (see data/generate_synthetic_data.R) and the
# schema described in data/metadata.json.
#
# Design: demo mode does NOT run Stage 1's real cohort-build/Cox-screen logic
# (that needs a live NHIS/JMDC connection, which demo mode has none of) --
# it only seeds Stage 1's OUTPUT as a normal file checkpoint
# (results/stage1_drug_screen/<run_id>/*.rds), validated against
# metadata.json's declared schema. Stage 2 then runs completely unmodified:
# it already prefers a live CKM_DRUG connection and falls back to the latest
# file checkpoint via core/R/db_common.R's read_table_or_checkpoint() when
# none is reachable -- in demo mode, connect_db() is still attempted and
# still fails (no real server), so Stage 2 picks up the seeded checkpoint
# through the exact same fallback path a real DB outage would use. Demo mode
# is therefore a Stage-1 substitute only; it required zero changes to Stage
# 2's own code.
# ==============================================================================

library(jsonlite)

#' Reads data/metadata.json from the pipeline root.
load_metadata <- function(pipeline_root) {
  path <- file.path(pipeline_root, "data", "metadata.json")
  if (!file.exists(path)) {
    stop("data/metadata.json not found at '", path, "' -- required for --demo mode.")
  }
  jsonlite::fromJSON(path, simplifyVector = TRUE)
}

#' Checks that every column metadata.json declares for `table_name`/`country`
#' is actually present in `df`. Warns (does not stop) on EXTRA columns --
#' only missing declared columns are an error, since core/R/*.R adapters are
#' free to select a subset. Called before a demo checkpoint is used, so a
#' stale metadata.json or a stale synthetic-data regeneration is caught
#' loudly instead of surfacing as a confusing downstream NA/column-not-found
#' error several pipeline stages later.
validate_against_schema <- function(df, table_name, country, metadata, log = message) {
  table_meta <- metadata$tables[[table_name]][[country]]
  if (is.null(table_meta) || is.null(table_meta$columns)) {
    log(sprintf("No schema declared for table '%s' / country '%s' in metadata.json -- skipping validation.", table_name, country))
    return(invisible(TRUE))
  }
  declared_cols <- setdiff(names(table_meta$columns), "note")
  missing <- setdiff(declared_cols, names(df))
  if (length(missing) > 0) {
    stop(sprintf(
      "Demo data for '%s' (%s) is missing column(s) declared in metadata.json: %s -- regenerate with data/generate_synthetic_data.R or fix metadata.json.",
      table_name, country, paste(missing, collapse = ", ")
    ))
  }
  log(sprintf("Schema check OK: '%s' (%s) has all %d column(s) metadata.json declares.", table_name, country, length(declared_cols)))
  invisible(TRUE)
}

#' Copies the synthetic seed data for `country` into a fresh
#' results/stage1_drug_screen/<run_id>/ checkpoint (same file shape
#' write_checkpoint() produces for a real Stage 1 run), after validating
#' each table against metadata.json's declared schema. Returns the
#' synthetic cox_results data frame (for the caller's own logging/parity
#' with a real Stage 1 run -- Stage 2 re-reads the checkpoint independently
#' via core/R/candidate_drugs.R, it does not receive this return value).
seed_stage1_checkpoint_from_demo <- function(country, pipeline_root, stage1_results_path, log = message) {
  metadata <- load_metadata(pipeline_root)
  demo <- metadata$demo_data[[country]]
  if (is.null(demo)) {
    stop("data/metadata.json has no demo_data entry for country '", country, "'.")
  }

  log(sprintf("Seeding Stage 1 checkpoint from synthetic demo data (metadata.json, seed=%s, n_patients=%s)...",
              metadata$demo_data$seed, metadata$demo_data$n_patients))

  tables <- list(
    cohort_processed = demo$cohort_processed,
    drug_exposure     = demo$drug_exposure,
    cox_results       = demo$cox_results
  )

  out <- list()
  for (nm in names(tables)) {
    rel_path <- tables[[nm]]
    abs_path <- file.path(pipeline_root, rel_path)
    if (!file.exists(abs_path)) {
      stop(sprintf(
        "Demo data file '%s' (declared in metadata.json for '%s'/'%s') not found -- run: Rscript data/generate_synthetic_data.R",
        rel_path, nm, country
      ))
    }
    df <- readRDS(abs_path)
    validate_against_schema(df, nm, country, metadata, log = log)
    write_checkpoint(stage1_results_path, nm, df)
    out[[nm]] <- df
  }

  log(sprintf("Demo checkpoint written to %s (cohort_processed: %d rows, drug_exposure: %d rows, cox_results: %d rows, %d candidate).",
              stage1_results_path, nrow(out$cohort_processed), nrow(out$drug_exposure),
              nrow(out$cox_results), sum(out$cox_results$candidate)))

  out$cox_results
}
