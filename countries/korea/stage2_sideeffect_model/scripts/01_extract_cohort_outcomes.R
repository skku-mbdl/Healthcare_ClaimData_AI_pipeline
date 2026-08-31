# ==============================================================================
# 01_extract_cohort_outcomes.R -- cohort extraction (Korea).
#
# Pulls the NHIS Sample Cohort CKM stage 2/3 cohort (CKM_DRUG.dbo.
# nhis_cohort_processed, already built by this project's own Stage 1 --
# CKM==1, washout-excluded, covariate-complete) via the DB when reachable,
# else the latest Stage 1 file checkpoint. All outcome/time/status
# derivation is now shared logic (core/R/extract_cohort_outcomes.R's
# build_cohort_outcomes()) -- this file only owns the raw extraction +
# column renaming to the standardized patient_id.
#
# nhis_cohort_processed already has a real `death_date` column (from Stage
# 1's all_death_claims()), unlike Japan's cohort_processed (which only has
# died_in_data + followup_end) -- passed through as-is here.
#
# eGFR: nhis_cohort_processed's `egfr` column is already the 2021 CKD-EPI
# creatinine equation computed from serum_creatinine/age/sex (core/R/
# clinical_calcs.R's egfr()), the SAME equation used as the Cox covariate --
# no separate recomputation needed here.
# ==============================================================================

progress("Connecting to CKM_DRUG database (Server = SY_PC)...")
con <- tryCatch(connect_db(), error = function(e) {
  progress(sprintf("Could not open a CKM_DRUG connection (%s) -- will use Stage 1's file checkpoint instead.", conditionMessage(e)))
  NULL
})

cohort_raw <- read_table_or_checkpoint(con, "nhis_cohort_processed", stage1_results_root,
                                        checkpoint_name = "cohort_processed", log = progress)

if (nrow(cohort_raw) == 0) {
  stop("nhis_cohort_processed returned 0 rows -- check that countries/korea/stage1_drug_screen has been run.")
}

required_cols <- c(
  "PERSON_ID", "index_date", "age", "sex", "egfr", "systolic_bp",
  "smoking_habit", "drinking_habit", "exercise_habit",
  "hdl_cholesterol", "ldl_cholesterol",
  "followup_end", "died_in_data", "death_date",
  "stroke_event", "stroke_date",
  "heart_failure_event", "heart_failure_date",
  "chd_event", "chd_date",
  "cv_death_event", "cv_death_date"
)
missing_cols <- setdiff(required_cols, names(cohort_raw))
if (length(missing_cols) > 0) {
  stop("nhis_cohort_processed is missing expected column(s): ", paste(missing_cols, collapse = ", "))
}

cohort_raw <- cohort_raw %>% dplyr::rename(patient_id = PERSON_ID)

progress("Deriving baseline covariates and composite event dates...")
cohort <- build_cohort_outcomes(cohort_raw, cfg, log = progress)
baseline_covariates <- attr(cohort, "baseline_covariates")
