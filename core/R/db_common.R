# Database helpers shared by both countries' Stage 1 (drug screen) and
# Stage 2 (side-effect model) pipelines -- everything here is schema-agnostic
# (works the same regardless of which raw source DB a country reads from).
# Country-specific connection functions (connect_jmdc/connect_nhis/
# connect_db, all pointed at different Server/Database values) live in each
# country's own db_connect.R, not here.
#
# check_connection/write_ckm_table/push_window_table/raw_tbl/union_years_tbl
# are verbatim ports of the identical (or near-identical, see comments)
# functions duplicated across CKM_Drug/Japan and CKM_Drug/Korea/NHIS_SAMPLE's
# own utils/db_utils.R.

library(DBI)

#' Run a trivial query against a freshly-opened connection to confirm it's
#' actually usable. Without this, a dead-on-arrival connection (e.g. its
#' external pointer already invalid) only surfaces on the pipeline's first
#' real query, several steps in, with an error that doesn't say which
#' connection or step caused it.
#'
#' Retries a few times with a short pause before giving up: a Trusted_Connection
#' opened via dbConnect() has occasionally shown "external pointer is not
#' valid" on the very first query when checked immediately (e.g. when run via
#' source() with no gap between dbConnect() and first use), but not when the
#' same calls are run interactively with natural pauses between them --
#' consistent with the Windows-authenticated handshake not being fully
#' settled the instant dbConnect() returns.
check_connection <- function(con, label, retries = 3, wait_seconds = 1) {
  last_error <- NULL
  for (attempt in seq_len(retries)) {
    outcome <- tryCatch({
      dbGetQuery(con, "SELECT 1 AS ok")
      TRUE
    }, error = function(e) e)

    if (isTRUE(outcome)) return(invisible(TRUE))
    last_error <- outcome
    if (attempt < retries) Sys.sleep(wait_seconds)
  }
  stop(sprintf("Connectivity check against %s failed: %s", label, conditionMessage(last_error)), call. = FALSE)
}

#' Write a data frame to the CKM_DRUG database, replacing the table if it
#' already exists. Inf/-Inf/NaN (e.g. a Cox model's confidence interval
#' blowing up on a rare drug-outcome pair with near-zero exposed events)
#' become SQL NULL rather than Inf/NA -- SQL Server's ODBC driver rejects
#' Inf/NaN as an invalid value for a float column ("TDS RPC protocol stream
#' is incorrect"), failing the whole write. CSV/RDS checkpoint outputs (see
#' write_checkpoint() below) keep the true Inf, since those formats have no
#' such restriction.
write_ckm_table <- function(con, table_name, data) {
  data[] <- lapply(data, function(col) {
    if (is.numeric(col)) col[!is.finite(col)] <- NA
    col
  })
  dbWriteTable(con, table_name, data, overwrite = TRUE)
  invisible(table_name)
}

#' Push a small per-member window table (member_id/PERSON_ID + one or more
#' date bounds) to a temp table on `con`, so a large claims table can be
#' joined to it server-side instead of pulling per-subject date windows into
#' R. Used wherever a step needs claims within a per-member date range
#' (lookback period, observation period, ...).
push_window_table <- function(con, window_df, name = "#cohort_windows") {
  dplyr::copy_to(con, window_df, name = name, overwrite = TRUE, temporary = TRUE)
}

#' A single raw table for one (prefix, year), e.g. raw_tbl(con, "NHID_GJ",
#' 2010) -> tbl(con, "NHID_GJ_2010"). Returns NULL (with a warning, not an
#' error) if the table doesn't exist on `con` -- some prefix/year
#' combinations may legitimately not exist yet, and callers union across
#' many of these, so one missing year/type shouldn't abort the whole query.
#' Used by the Korea (NHIS) adapter's per-year table-naming convention; the
#' Japan (JMDC) adapter doesn't need this (JMDC's tables aren't split by
#' year).
raw_tbl <- function(con, prefix, year) {
  name <- sprintf("%s_%d", prefix, year)
  if (!DBI::dbExistsTable(con, name)) {
    warning(sprintf("Table %s not found on this connection -- skipped.", name), call. = FALSE)
    return(NULL)
  }
  dplyr::tbl(con, name)
}

#' UNION ALL a raw table across every year in `years`, for a single prefix.
#' For tables where every row already carries its own PERSON_ID (no KEY_SEQ
#' join needed first) -- see the Korea adapter's diagnosis_utils.R/
#' drug_utils.R for the per-year-join variant used where a join IS needed
#' first.
union_years_tbl <- function(con, prefix, years) {
  parts <- Filter(Negate(is.null), lapply(years, function(y) raw_tbl(con, prefix, y)))
  if (length(parts) == 0) {
    stop(sprintf("No %s_<year> tables found for years %s.", prefix, paste(range(years), collapse = "-")), call. = FALSE)
  }
  Reduce(dplyr::union_all, parts)
}

# ==============================================================================
# File-checkpoint fallback (new in this integrated pipeline).
#
# Both stages' run_pipeline.R ALWAYS write their key outputs to
# results/<run_id>/*.rds (in addition to the CKM_DRUG database, when
# reachable) -- this lets the pipeline be dry-run/smoke-tested without a live
# DB connection, and gives Stage 2 a fallback source of Stage 1's cohort /
# drug-exposure / candidate-drug output when the database isn't reachable.
# The DB stays the primary hand-off path when available (it's always
# current, and multiple downstream consumers can query it); the checkpoint
# is a portable, file-based mirror of the same tables, not a replacement.
# ==============================================================================

#' Write one named table as both .rds (exact round-trip, keeps true Inf/NA)
#' and .csv (human-inspectable) under results_path.
write_checkpoint <- function(results_path, name, data) {
  if (!dir.exists(results_path)) dir.create(results_path, recursive = TRUE)
  saveRDS(data, file.path(results_path, paste0(name, ".rds")))
  write.csv(data, file.path(results_path, paste0(name, ".csv")), row.names = FALSE)
  invisible(file.path(results_path, paste0(name, ".rds")))
}

#' Find the most recent results/<run_id>/ directory (run_id is a
#' YYYYMMDD_HHMMSS timestamp, so lexicographic order == chronological order)
#' that actually contains `name`.rds, under `results_root`. Returns NULL if
#' none exists, so callers can decide how to fail (e.g. "run Stage 1 first").
find_latest_checkpoint <- function(results_root, name) {
  if (!dir.exists(results_root)) return(NULL)
  run_dirs <- sort(list.dirs(results_root, recursive = FALSE), decreasing = TRUE)
  for (d in run_dirs) {
    candidate <- file.path(d, paste0(name, ".rds"))
    if (file.exists(candidate)) return(candidate)
  }
  NULL
}

#' Read `name` from the CKM_DRUG database if `con` is a live, working
#' connection; otherwise (or on any query error) fall back to the most
#' recent file checkpoint under `results_root`. Always logs which source was
#' actually used, per this repo's no-silent-fallback convention.
read_table_or_checkpoint <- function(con, table_name, results_root, checkpoint_name = table_name, log = message) {
  if (!is.null(con)) {
    out <- tryCatch({
      check_connection(con, table_name)
      DBI::dbReadTable(con, table_name)
    }, error = function(e) e)

    if (!inherits(out, "error")) {
      log(sprintf("Read '%s' from CKM_DRUG database (%d rows).", table_name, nrow(out)))
      return(out)
    }
    log(sprintf("Could not read '%s' from CKM_DRUG database (%s) -- falling back to file checkpoint.",
                table_name, conditionMessage(out)))
  } else {
    log(sprintf("No database connection supplied for '%s' -- using file checkpoint.", table_name))
  }

  path <- find_latest_checkpoint(results_root, checkpoint_name)
  if (is.null(path)) {
    stop(sprintf(
      "Could not read '%s': no database connection available and no results/<run_id>/%s.rds checkpoint found under '%s'. Run Stage 1 first.",
      table_name, checkpoint_name, results_root
    ), call. = FALSE)
  }
  log(sprintf("Read '%s' from file checkpoint: %s", table_name, path))
  readRDS(path)
}
