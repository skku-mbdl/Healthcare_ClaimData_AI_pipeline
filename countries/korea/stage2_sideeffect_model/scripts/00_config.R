# ==============================================================================
# 00_config.R -- side-effect-prediction pipeline (Stage 2): constants + DB
# connection + logging (Korea / NHIS Sample Cohort).
#
# Predicts risk of MACE/ASCVD/HF/CV_DEATH from an individual's baseline
# drug-exposure history, in the NHIS Sample Cohort CKM stage 2/3 cohort
# already built by this project's own Stage 1
# (countries/korea/stage1_drug_screen), which writes CKM_DRUG.dbo.
# nhis_cohort_processed, CKM_DRUG.dbo.nhis_drug_exposure, and CKM_DRUG.dbo.
# nhis_cox_results -- the "nhis_" prefix is kept (same as the original
# CKM_Drug/Korea/NHIS_SAMPLE pipeline) because CKM_DRUG is ONE shared
# database Japan's Stage 1 also writes to, with its own unprefixed table
# names (cohort_processed/drug_exposure/cox_results) -- see
# countries/korea/run_pipeline.R and countries/japan/run_pipeline.R.
#
# Drug predictors are per-ATC-code adherence flags (PDC > 0.8 during the
# 1-year lookback) -- Korea's Stage 1 only ever computes exposure at ATC
# granularity (see stage1_drug_screen/scripts/06_drug_exposure.R), one drug
# meaning one ATC code, per the study protocol's "Exposure" section.
#
# T0_YEARS is kept at 10 for methodological consistency with the Japan
# pipeline, even though NHIS_SAMPLE's actual observation period only runs
# index (2009-2012) through 2013-12-31 (max ~1-4 years actually observed).
# A 10-year CIF here extrapolates the fitted baseline hazard well past any
# observed event -- treat absolute risk numbers from this pipeline as
# considerably less reliable than Japan's, and prefer the relative (hazard
# ratio / confirmed-drug) output over the absolute 10-year risk if the two
# countries' numbers are ever compared directly.
#
# Run with working directory = countries/korea/ (countries/korea/
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
# DB connection -- same CKM_DRUG database Stage 1 writes to.
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
# drug_exposure's `exposed` flag) are dropped before ever reaching R --
# Boruta tractability; a random-forest importance test on a near-constant
# column is noise. Applied SQL-side in 02_extract_drug_exposure.R.
MIN_DRUG_PREVALENCE <- 0.01

T0_YEARS <- 10

# Administrative end of follow-up: Stage 1's protocol-defined observation
# cutoff, already used upstream to cap cohort_processed$followup_end.
STUDY_END_DATE <- as.Date("2013-12-31")

# Boruta settings (fixed seed for reproducibility of the confirmed/rejected
# variable split across runs).
BORUTA_SEED <- 20260818
BORUTA_MAX_RUNS <- 250
BORUTA_P_VALUE <- 0.01

BORUTA_NUM_THREADS <- if (packageVersion("Boruta") >= "5.0.0") {
  max(1, parallel::detectCores() - 1)
} else {
  NULL
}

# Train/test split: Boruta selection and Cox fitting run on dfd_train only;
# validation runs on dfd_test only. Stratified on MACE status.
TRAIN_FRACTION <- 0.7
SPLIT_SEED <- 20260819

# SHAP: same Kernel SHAP approach as Japan's pipeline (predict_cif_10yr() is
# a nonlinear competing-risks combination of two Cox models, so a
# closed-form linear SHAP doesn't apply).
SHAP_MAX_ROWS <- 200
SHAP_BG_N <- 100
SHAP_SEED <- 20260820

# Bundled as a plain list so core/R functions that take `cfg` can read
# whichever fields they need without every one of them threading through
# a dozen separate global variables.
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
