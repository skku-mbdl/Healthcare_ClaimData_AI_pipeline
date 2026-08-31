# ==============================================================================
# generate_synthetic_data.R -- builds the demo-mode synthetic dataset used by
# `Rscript countries/<country>/run_pipeline.R --demo` (see data/metadata.json
# and core/R/demo_data.R).
#
# NOT REAL PATIENT DATA. This is a synthetic stand-in for what Stage 1
# (cohort-build + Cox screen) would have written to the CKM_DRUG database,
# used only so the pipeline can be run and pilot-tested end-to-end without a
# live SQL Server -- neither NHIS/JMDC raw claims nor a CKM_DRUG connection
# are reachable in the environment this pipeline was built in.
#
# The simulation embeds genuine signal (not pure noise): for each country,
# 3 of 15 ATC codes are designated "risky" and are given a real excess hazard
# for CHD/Stroke/Heart-Failure/CV-death event times; the accompanying
# cox_results table marks exactly those 3 as candidate == TRUE. This lets a
# pilot run of Stage 2 (Boruta -> Cox -> validation -> SHAP) exercise real
# variable-selection and model-fitting behavior against data with a known
# answer, rather than pure noise Boruta would confirm nothing on.
#
# Re-run this script (Rscript data/generate_synthetic_data.R) to regenerate;
# SEED below is fixed for reproducibility.
# ==============================================================================

suppressPackageStartupMessages(library(dplyr))

SEED <- 20260901
N_PATIENTS <- 2000
N_ATC_CODES <- 15
N_RISKY <- 3
set.seed(SEED)

data_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
if (length(data_dir) == 0 || !nzchar(data_dir)) data_dir <- "data"  # fallback if sourced interactively

make_atc_codes <- function(n, prefix_pool) {
  sprintf("%s%02d", sample(prefix_pool, n, replace = TRUE), sample(10:99, n))
}

#' Simulate one country's Stage-1-shaped output.
#' `id_col`/`atc_col` match that country's real cohort_processed/drug_exposure
#' column names (see each country's stage1_drug_screen adapter).
#' `has_death_date`: Korea's nhis_cohort_processed has a real death_date
#' column; Japan's cohort_processed does not (see CLAUDE.md's "Known
#' asymmetries").
simulate_country <- function(id_col, atc_col, sex_codes, has_death_date, study_start, study_end) {

  ids <- sprintf("SYN%05d", seq_len(N_PATIENTS))
  index_date <- study_start + sample(0:365, N_PATIENTS, replace = TRUE)
  age <- pmin(pmax(round(rnorm(N_PATIENTS, 62, 10)), 30), 90)
  sex <- sample(sex_codes, N_PATIENTS, replace = TRUE)
  is_male <- sex == sex_codes[1]
  egfr <- pmin(pmax(rnorm(N_PATIENTS, 70, 18), 8), 130)
  systolic_bp <- pmin(pmax(round(rnorm(N_PATIENTS, 135, 18)), 90), 220)
  smoking_habit <- sample(1:3, N_PATIENTS, replace = TRUE, prob = c(0.5, 0.2, 0.3))
  drinking_habit <- sample(1:3, N_PATIENTS, replace = TRUE, prob = c(0.4, 0.3, 0.3))
  exercise_habit <- sample(1:3, N_PATIENTS, replace = TRUE, prob = c(0.4, 0.3, 0.3))
  hdl_cholesterol <- pmin(pmax(rnorm(N_PATIENTS, 48, 12), 15), 100)
  ldl_cholesterol <- pmin(pmax(rnorm(N_PATIENTS, 118, 28), 30), 250)

  followup_years <- pmin(pmax(rexp(N_PATIENTS, 1 / 5), 0.1), 10)
  followup_end <- pmin(index_date + round(followup_years * 365.25), study_end)

  # ---- Drug exposure: N_ATC_CODES codes, N_RISKY of them carry real signal ----
  atc_codes <- unique(make_atc_codes(N_ATC_CODES * 2, c("A10", "C02", "C03", "C07", "C08", "C09", "N05", "M01", "J01")))[seq_len(N_ATC_CODES)]
  risky_atc <- sample(atc_codes, N_RISKY)
  prevalence <- runif(N_ATC_CODES, 0.03, 0.28)
  names(prevalence) <- atc_codes

  exposure_mat <- sapply(atc_codes, function(code) rbinom(N_PATIENTS, 1, prevalence[code]))
  colnames(exposure_mat) <- atc_codes
  risky_exposure_any <- rowSums(exposure_mat[, risky_atc, drop = FALSE]) > 0

  # ---- Simulate 4 outcome event times via exponential hazards with a real
  # baseline-covariate + risky-drug linear predictor, censored at followup_end.
  sim_time_to_event <- function(base_rate, lp) {
    # base_rate is already a per-day rate (e.g. 1/(1.5*365.25)), and rexp()
    # with a per-day rate already returns values in days -- do NOT multiply
    # by 365.25 again here (an earlier version of this function did, which
    # inflated every simulated event time ~365x and left every outcome with
    # almost no events within any realistic follow-up window).
    rate <- base_rate * exp(lp)
    rexp(N_PATIENTS, rate)
  }
  lp_common <- 0.015 * (age - 62) + 0.10 * (!is_male) + 0.01 * (systolic_bp - 135) - 0.01 * (egfr - 70)
  lp_risky <- 0.55 * risky_exposure_any  # true excess hazard for risky drugs

  # Mean event times are deliberately short (not clinically realistic) so
  # that even Korea's short study window (~1-4 years of actual follow-up,
  # see stage1_drug_screen/scripts/utils/config.R's observation_end_cutoff)
  # produces a healthy number of events for EACH of the 4 outcomes -- a
  # first version of this generator used realistic 10-20-year mean times,
  # which left some outcomes (e.g. heart failure) with only 1-2 events in a
  # 1400-row training split, crashing Boruta with "need at least two
  # classes to do classification." This is purely a synthetic-data
  # calibration choice for exercising the pipeline's code paths, not a
  # claim about real-world CKM event rates.
  chd_days   <- sim_time_to_event(1 / (6 * 365.25), lp_common + lp_risky)
  stroke_days <- sim_time_to_event(1 / (8 * 365.25), lp_common + 0.9 * lp_risky)
  hf_days    <- sim_time_to_event(1 / (10 * 365.25), lp_common + 0.7 * lp_risky)
  noncvd_death_days <- sim_time_to_event(1 / (14 * 365.25), 0.02 * (age - 62))
  cv_death_days <- sim_time_to_event(1 / (12 * 365.25), lp_common + 0.6 * lp_risky)

  followup_len <- as.numeric(followup_end - index_date)

  chd_event <- chd_days <= followup_len
  chd_date <- as.Date(ifelse(chd_event, index_date + round(chd_days), NA), origin = "1970-01-01")
  stroke_event <- stroke_days <= followup_len
  stroke_date <- as.Date(ifelse(stroke_event, index_date + round(stroke_days), NA), origin = "1970-01-01")
  heart_failure_event <- hf_days <= followup_len
  heart_failure_date <- as.Date(ifelse(heart_failure_event, index_date + round(hf_days), NA), origin = "1970-01-01")

  # Whichever of {CV death, non-CV death} would happen first, if either happens
  # before administrative censoring -- a patient can only die once.
  death_is_cv <- cv_death_days < noncvd_death_days
  death_days <- pmin(cv_death_days, noncvd_death_days)
  died_in_data <- death_days <= followup_len
  death_date <- as.Date(ifelse(died_in_data, index_date + round(death_days), NA), origin = "1970-01-01")
  cv_death_event <- died_in_data & death_is_cv
  cv_death_date <- as.Date(ifelse(cv_death_event, death_date, NA), origin = "1970-01-01")

  # Death (of either cause) ends observation -- cap followup_end at the death
  # date, same convention as the real Stage 1 pipelines.
  followup_end <- as.Date(ifelse(died_in_data, pmin(death_date, followup_end), followup_end), origin = "1970-01-01")

  cohort <- data.frame(
    id = ids, index_date = index_date, age = age, sex = sex,
    egfr = round(egfr, 1), systolic_bp = systolic_bp,
    smoking_habit = smoking_habit, drinking_habit = drinking_habit, exercise_habit = exercise_habit,
    hdl_cholesterol = round(hdl_cholesterol, 1), ldl_cholesterol = round(ldl_cholesterol, 1),
    followup_end = followup_end, died_in_data = died_in_data,
    stroke_event = stroke_event, stroke_date = stroke_date,
    heart_failure_event = heart_failure_event, heart_failure_date = heart_failure_date,
    chd_event = chd_event, chd_date = chd_date,
    cv_death_event = cv_death_event, cv_death_date = cv_death_date,
    stringsAsFactors = FALSE
  )
  if (has_death_date) cohort$death_date <- death_date
  names(cohort)[names(cohort) == "id"] <- id_col

  # ---- drug_exposure: long format, one row per (patient, exposed ATC code) ----
  exposure_long <- which(exposure_mat == 1, arr.ind = TRUE)
  drug_exposure <- data.frame(
    id = ids[exposure_long[, "row"]],
    atc = atc_codes[exposure_long[, "col"]],
    days_covered = sample(300:365, nrow(exposure_long), replace = TRUE),
    pdc = round(runif(nrow(exposure_long), 0.81, 1.0), 3),
    exposed = TRUE,
    stringsAsFactors = FALSE
  )
  names(drug_exposure)[names(drug_exposure) == "id"] <- id_col
  names(drug_exposure)[names(drug_exposure) == "atc"] <- atc_col

  # ---- cox_results: Stage 1's unified Cox-screen output shape (core/R/
  # cox_screen.R), same column names for both countries. Candidate flags
  # match the true risky_atc set used to simulate outcomes above -- this is
  # what makes the demo dataset useful for testing Stage 2's variable
  # selection, not just its plumbing.
  outcomes <- c("chd", "stroke", "heart_failure", "cv_death")
  cox_results <- expand.grid(atc_code = atc_codes, outcome = outcomes, stringsAsFactors = FALSE) %>%
    mutate(
      is_risky = atc_code %in% risky_atc,
      hr = ifelse(is_risky, round(runif(n(), 1.3, 2.3), 3), round(runif(n(), 0.75, 1.3), 3)),
      ci_lower = round(pmax(hr - runif(n(), 0.15, 0.35), 0.5), 3),
      ci_upper = round(hr + runif(n(), 0.15, 0.5), 3),
      p_value = ifelse(is_risky, round(runif(n(), 0.001, 0.045), 4), round(runif(n(), 0.08, 0.95), 4)),
      n = N_PATIENTS,
      n_event = sample(20:300, n(), replace = TRUE),
      candidate = is_risky & hr > 1.00 & p_value < 0.05
    ) %>%
    select(atc_code, outcome, hr, ci_lower, ci_upper, p_value, n, n_event, candidate)

  list(cohort = cohort, drug_exposure = drug_exposure, cox_results = cox_results, risky_atc = risky_atc)
}

progress <- function(...) cat(sprintf("[generate_synthetic_data] %s\n", sprintf(...)))

study_start <- as.Date("2010-01-01")

progress("Simulating Korea (n=%d patients, %d ATC codes, %d risky)...", N_PATIENTS, N_ATC_CODES, N_RISKY)
korea <- simulate_country(
  id_col = "PERSON_ID", atc_col = "ATC_CODE", sex_codes = c("1", "2"),
  has_death_date = TRUE, study_start = study_start, study_end = as.Date("2013-12-31")
)

progress("Simulating Japan (n=%d patients, %d ATC codes, %d risky)...", N_PATIENTS, N_ATC_CODES, N_RISKY)
japan <- simulate_country(
  id_col = "member_id", atc_col = "who_atc_code", sex_codes = c("Male", "Female"),
  has_death_date = FALSE, study_start = study_start, study_end = as.Date("2022-07-31")
)

write_country <- function(sim, country_dir, log_label) {
  dir.create(country_dir, recursive = TRUE, showWarnings = FALSE)
  for (nm in c("cohort", "drug_exposure", "cox_results")) {
    out_name <- if (nm == "cohort") "cohort_processed" else nm
    saveRDS(sim[[nm]], file.path(country_dir, paste0(out_name, ".rds")))
    write.csv(sim[[nm]], file.path(country_dir, paste0(out_name, ".csv")), row.names = FALSE)
  }
  progress("%s: wrote cohort_processed (%d rows), drug_exposure (%d rows), cox_results (%d rows, %d candidate) -> %s",
            log_label, nrow(sim$cohort), nrow(sim$drug_exposure), nrow(sim$cox_results),
            sum(sim$cox_results$candidate), country_dir)
  progress("%s: true risky ATC codes = %s", log_label, paste(sim$risky_atc, collapse = ", "))
}

write_country(korea, file.path(data_dir, "synthetic", "korea"), "Korea")
write_country(japan, file.path(data_dir, "synthetic", "japan"), "Japan")

progress("Done.")
