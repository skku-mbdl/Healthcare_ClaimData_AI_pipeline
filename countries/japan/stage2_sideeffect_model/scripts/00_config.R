# ==============================================================================
# 00_config.R -- side-effect-prediction pipeline (Stage 2): constants + DB
# connection + logging (Japan / JMDC).
#
# Predicts risk of MACE/ASCVD/HF/CV_DEATH from an individual's baseline
# drug-exposure history, in the JMDC CKM stage 2/3 cohort already built by
# this project's own Stage 1 (countries/japan/stage1_drug_screen), which
# writes CKM_DRUG.dbo.cohort_processed and CKM_DRUG.dbo.drug_exposure.
#
# ARCHITECTURE CHANGE from the original CKM_PREVENT/JAPAN/SideEffect_Model:
# that project read its cohort from PREVENT_CKM.dbo.jmdc_cohort_final, a
# table built by jmdc_ckm_cohort.sql -- part of CKM_PREVENT/JAPAN/PREVENT,
# which this integrated pipeline deliberately excludes (see the top-level
# CLAUDE.md). This Stage 2 now reads CKM_DRUG.dbo.cohort_processed instead
# -- Stage 1's OWN output, mirroring how Korea's Stage 2 always worked. This
# is a real methodology change, not just plumbing: jmdc_cohort_final used
# the Japanese Society of Nephrology (2009) eGFR equation for CKM
# eligibility/staging, while cohort_processed (via core/R/clinical_calcs.R)
# uses the 2021 CKD-EPI equation -- a different formula at the eligibility
# step can mean a genuinely different set of patients. This was a deliberate
# choice (confirmed with the project owner), not an oversight.
#
# Drug predictors are WHO ATC codes (Drug_master.who_atc_code, as computed
# by Stage 1's own PDC calculation), restricted to ATC codes flagged
# candidate == TRUE in Stage 1's own Cox screen (core/R/candidate_drugs.R) --
# loaded automatically now (DB or file-checkpoint fallback) instead of a
# manually-copied cox_results.csv. The old event별_ATC정리.xlsx
# clinically-curated allow-list is NOT used any more (it was already
# superseded by Stage 1's own Cox screen output before this integration,
# per the original project's own 2026-08-30 decision) -- see
# reference_data/japan/event별_ATC정리.xlsx if that clinical curation is
# ever needed again as a cross-check.
#
# T0_YEARS = 10 matches the diagnosis-ascertainment window Stage 1's cohort
# already used (10-year outcome window).
#
# Run with working directory = countries/japan/ (countries/japan/
# run_pipeline.R sets this before sourcing anything under
# stage2_sideeffect_model/scripts/). This file is sourced first and must
# not be run standalone before 01_extract_cohort_outcomes.R etc.
# ==============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

library(DBI)
library(odbc)
library(dplyr)
library(lubridate)
library(survival)

if (!requireNamespace("Boruta", quietly = TRUE)) install.packages("Boruta")
library(Boruta)

if (!requireNamespace("kernelshap", quietly = TRUE)) install.packages("kernelshap")

# ------------------------------------------------------------------------------
# Logging.
# ------------------------------------------------------------------------------
log_dir <- file.path("logs", "stage2_sideeffect_model")
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
log_file <- file.path(log_dir, paste0("sideeffect_predict_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))

log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

progress <- function(msg) {
  mem_mb <- sum(gc(full = FALSE)[, 2])
  cat(sprintf("[%s] (mem: %.0f MB) %s\n", format(Sys.time(), "%H:%M:%S"), mem_mb, msg))
  flush(stdout())
}

progress(paste("Side-effect-prediction run started:", format(Sys.time())))
progress(paste("Log file:", log_file))

output_dir <- file.path("results", "stage2_sideeffect_model", "sideeffect_predict")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

figures_dir <- file.path("figures", "stage2_sideeffect_model", "sideeffect_predict")
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

# Stage 1's own results root, for the file-checkpoint fallback (see
# core/R/db_common.R's read_table_or_checkpoint()) when CKM_DRUG isn't
# reachable.
stage1_results_root <- file.path("results", "stage1_drug_screen")

# ------------------------------------------------------------------------------
# DB connection -- CKM_DRUG, the SAME database Stage 1 writes to (NOT
# PREVENT_CKM -- see the architecture-change note above).
# ------------------------------------------------------------------------------
connect_db <- function() {
  dbConnect(
    odbc(),
    Driver = "ODBC Driver 17 for SQL Server",
    Server = "SY_PC",
    Database = "CKM_DRUG",
    Trusted_Connection = "yes"
  )
}

# ------------------------------------------------------------------------------
# Stage-1-outcome mapping for this pipeline's 4 composite outcomes -- see
# core/R/boruta_selection.R's header for the rationale (identical for both
# countries, since both Stage 1 pipelines define exactly stroke/
# heart_failure/chd/cv_death).
# ------------------------------------------------------------------------------
cox_outcome_map <- list(
  mace    = c("chd", "stroke", "heart_failure", "cv_death"),
  ascvd   = c("chd", "stroke"),
  hf      = c("heart_failure"),
  cvdeath = c("cv_death")
)

# ------------------------------------------------------------------------------
# Cohort / windowing constants
# ------------------------------------------------------------------------------

DRUG_LOOKBACK_DAYS <- 365

# ATC codes exposed to fewer than this fraction of the cohort (per
# drug_exposure's `exposed` flag) are dropped before ever reaching R.
MIN_DRUG_PREVALENCE <- 0.01

T0_YEARS <- 10

# Administrative end of follow-up: Stage 1's own protocol-defined
# observation cutoff, already used upstream to cap
# cohort_processed$followup_end.
STUDY_END_DATE <- as.Date("2022-07-31")

BORUTA_SEED <- 20260729
BORUTA_MAX_RUNS <- 250
BORUTA_P_VALUE <- 0.01

BORUTA_NUM_THREADS <- if (packageVersion("Boruta") >= "5.0.0") {
  max(1, parallel::detectCores() - 1)
} else {
  NULL
}

TRAIN_FRACTION <- 0.7
SPLIT_SEED <- 20260730

SHAP_MAX_ROWS <- 200
SHAP_BG_N <- 100
SHAP_SEED <- 20260808

cfg <- list(
  T0_YEARS = T0_YEARS,
  STUDY_END_DATE = STUDY_END_DATE,
  TRAIN_FRACTION = TRAIN_FRACTION,
  SPLIT_SEED = SPLIT_SEED,
  BORUTA_SEED = BORUTA_SEED,
  BORUTA_MAX_RUNS = BORUTA_MAX_RUNS,
  BORUTA_P_VALUE = BORUTA_P_VALUE,
  BORUTA_NUM_THREADS = BORUTA_NUM_THREADS,
  SHAP_MAX_ROWS = SHAP_MAX_ROWS,
  SHAP_BG_N = SHAP_BG_N,
  SHAP_SEED = SHAP_SEED
)
