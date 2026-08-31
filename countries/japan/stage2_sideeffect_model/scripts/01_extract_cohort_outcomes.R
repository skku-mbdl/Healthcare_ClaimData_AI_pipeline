# ==============================================================================
# 01_extract_cohort_outcomes.R -- cohort extraction (Japan).
#
# Pulls the JMDC CKM stage 2/3 cohort (CKM_DRUG.dbo.cohort_processed,
# already built by this project's own Stage 1) via the DB when reachable,
# else the latest Stage 1 file checkpoint. See 00_config.R's header for why
# this reads cohort_processed rather than the old PREVENT_CKM.dbo.
# jmdc_cohort_final. All outcome/time/status derivation is now shared logic
# (core/R/extract_cohort_outcomes.R's build_cohort_outcomes()) -- this file
# only owns the raw extraction + column renaming/approximation needed to
# match that shared function's input contract.
#
# death_date approximation: cohort_processed does NOT carry a separate
# death date (unlike Korea's nhis_cohort_processed) -- it only has
# died_in_data (from Enrollment.withdrawal_death_flag) and followup_end
# (already capped at death when died_in_data is TRUE, or at the
# administrative cutoff otherwise -- see Stage 1's 05_outcomes.R). This is
# approximated here as death_date = followup_end when died_in_data is TRUE,
# NA otherwise -- the same month/year-level precision Stage 1's own
# followup_end already has, not a fabricated finer-grained date.
#
# eGFR: cohort_processed's `egfr` column is already the 2021 CKD-EPI
# creatinine equation (core/R/clinical_calcs.R's egfr()) -- the SAME
# equation used as the Cox covariate here, so no separate recomputation is
# needed (unlike the old jmdc_cohort_final-based version, which had to
# compute a second `egfr_ckdepi` column alongside jmdc_cohort_final's own
# Japan-Society-of-Nephrology-equation eGFR).
# ==============================================================================

progress("Connecting to CKM_DRUG database (Server = SY_PC)...")
con <- tryCatch(connect_db(), error = function(e) {
  progress(sprintf("Could not open a CKM_DRUG connection (%s) -- will use Stage 1's file checkpoint instead.", conditionMessage(e)))
  NULL
})

cohort_raw <- read_table_or_checkpoint(con, "cohort_processed", stage1_results_root, log = progress)

if (nrow(cohort_raw) == 0) {
  stop("cohort_processed returned 0 rows -- check that countries/japan/stage1_drug_screen has been run.")
}

required_cols <- c(
  "member_id", "index_date", "age", "sex", "egfr", "systolic_bp",
  "smoking_habit", "drinking_habit", "exercise_habit",
  "hdl_cholesterol", "ldl_cholesterol",
  "followup_end", "died_in_data",
  "stroke_event", "stroke_date",
  "heart_failure_event", "heart_failure_date",
  "chd_event", "chd_date",
  "cv_death_event", "cv_death_date"
)
missing_cols <- setdiff(required_cols, names(cohort_raw))
if (length(missing_cols) > 0) {
  stop("cohort_processed is missing expected column(s): ", paste(missing_cols, collapse = ", "))
}

cohort_raw <- cohort_raw %>%
  dplyr::rename(patient_id = member_id) %>%
  dplyr::mutate(
    followup_end = as.Date(followup_end),
    death_date = if_else(died_in_data, followup_end, as.Date(NA))
  )

progress("Deriving baseline covariates and composite event dates...")
cohort <- build_cohort_outcomes(cohort_raw, cfg, log = progress)
baseline_covariates <- attr(cohort, "baseline_covariates")
