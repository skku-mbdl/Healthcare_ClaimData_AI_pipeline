# ==============================================================================
# Stage 2 shared: derive MACE/ASCVD/HF/CV_DEATH competing-risks time/status
# from a country's already-extracted CKM cohort.
#
# This became shareable only because Japan's Stage 2 now reads its cohort
# from CKM_DRUG.dbo.cohort_processed (CKM_Drug/Japan's OWN Stage 1 output)
# instead of PREVENT_CKM.dbo.jmdc_cohort_final (built by the excluded
# CKM_PREVENT/JAPAN/PREVENT project's jmdc_ckm_cohort.sql). cohort_processed
# and Korea's nhis_cohort_processed both define exactly four outcomes (CHD,
# Stroke, Heart Failure, CV death; no separate MI/Fatal_CHD split like
# jmdc_cohort_final had), so the composite-outcome/time/status derivation
# below is now byte-identical logic for both countries -- previously it
# genuinely differed (Japan's had to bridge jmdc_cohort_final's 6-outcome
# shape). Each country's own 01_extract_cohort_outcomes.R still owns the raw
# DB query + column renaming (schema differs there), then hands off a
# standardized data frame to build_cohort_outcomes() below.
#
# Required columns on `cohort_raw`:
#   patient_id, index_date, age, sex, egfr, systolic_bp, smoking_habit,
#   drinking_habit, exercise_habit, hdl_cholesterol, ldl_cholesterol,
#   followup_end, died_in_data (logical), death_date (Date; NA if never
#   died -- for Japan, approximated as followup_end when died_in_data is
#   TRUE, since cohort_processed doesn't carry a separate death date; see
#   Japan's stage2 01_extract_cohort_outcomes.R for that approximation),
#   stroke_event/stroke_date, heart_failure_event/heart_failure_date,
#   chd_event/chd_date, cv_death_event/cv_death_date.
#
# `cfg` must have STUDY_END_DATE and T0_YEARS (see each country's stage2
# config.R).
#
# Composite outcomes (per project decision, same for both countries since
# neither protocol splits CHD further):
#   ASCVD := CHD | Stroke        (CHD's ICD-10 range I20-I25 already covers
#                                  what a JMDC-native extract would split
#                                  into MI + Fatal_CHD)
#   MACE  := ASCVD | HF | CV_DEATH
#
# Competing-risks status coding used throughout this repo: 0 = censored,
# 1 = event of interest, 2 = competing event (non-CV death).
# ==============================================================================

library(dplyr)

# Row-wise min across a handful of Date columns, ignoring NA (vectorized
# pmin(), not a per-row apply()).
.num_pmin_date <- function(...) {
  args <- lapply(list(...), as.numeric)
  out <- do.call(pmin, c(args, list(na.rm = TRUE)))
  as.Date(out, origin = "1970-01-01")
}

make_time_status <- function(event_date, index_date, end_date, noncvd_death_date) {
  event_date <- as.Date(event_date)
  index_date <- as.Date(index_date)
  end_date <- as.Date(end_date)
  noncvd_death_date <- as.Date(noncvd_death_date)

  event_before_end <- !is.na(event_date) & event_date <= end_date
  death_before_end <- !is.na(noncvd_death_date) & noncvd_death_date <= end_date

  status <- ifelse(
    event_before_end & (!death_before_end | event_date <= noncvd_death_date),
    1L,
    ifelse(death_before_end, 2L, 0L)
  )

  final_date <- end_date
  final_date[status == 1L] <- event_date[status == 1L]
  final_date[status == 2L] <- noncvd_death_date[status == 2L]

  time <- as.numeric(final_date - index_date) / 365.25

  neg <- !is.na(time) & time < 0
  if (any(neg)) {
    warning(sum(neg), " observation(s) have a negative survival time ",
            "(event/death/end date before index date); treating them as missing.")
    time[neg] <- NA_real_
    status[neg] <- NA_integer_
  }

  data.frame(time = time, status = status)
}

truncate_at_t0 <- function(time, status, t0) {
  time <- as.numeric(time)
  status <- as.integer(status)
  time_10 <- pmin(time, t0)
  status_10 <- ifelse(time <= t0, status, 0L)
  data.frame(time_10 = time_10, status_10 = as.integer(status_10))
}

#' Returns `cohort_raw` with composite event dates, per-outcome
#' time_<x>/status_<x>/time_<x>_10/status_<x>_10 columns added, plus
#' `baseline_covariates` (the standard 9-covariate list) as an attribute.
build_cohort_outcomes <- function(cohort_raw, cfg, log = message) {

  cohort <- cohort_raw %>%
    mutate(
      index_date  = as.Date(index_date),
      end_date    = cfg$STUDY_END_DATE,
      stroke_date = as.Date(stroke_date),
      chd_date    = as.Date(chd_date),
      hf_date     = as.Date(heart_failure_date),
      cv_date     = as.Date(cv_death_date),
      death_date  = as.Date(death_date),
      # died_in_data is TRUE for ANY death (all-cause); cv_death_event is
      # TRUE only for the circulatory-chapter subset. A non-CV death is
      # therefore "died, but not counted as cv_death".
      nocvd             = if_else(died_in_data & !cv_death_event, 1, 0),
      noncvd_death_date = if_else(nocvd == 1, death_date, as.Date(NA)),
      sex             = factor(sex),
      smoking_habit   = factor(smoking_habit),
      drinking_habit  = factor(drinking_habit),
      exercise_habit  = factor(exercise_habit)
    )

  cohort$ascvd_date <- .num_pmin_date(cohort$chd_date, cohort$stroke_date)
  cohort$mace_date  <- .num_pmin_date(cohort$ascvd_date, cohort$hf_date, cohort$cv_date)

  outcome_event_dates <- list(
    mace    = cohort$mace_date,
    ascvd   = cohort$ascvd_date,
    hf      = cohort$hf_date,
    cvdeath = cohort$cv_date
  )

  for (nm in names(outcome_event_dates)) {
    ts <- make_time_status(outcome_event_dates[[nm]], cohort$index_date, cohort$end_date, cohort$noncvd_death_date)
    ts10 <- truncate_at_t0(ts$time, ts$status, t0 = cfg$T0_YEARS)
    cohort[[paste0("time_", nm)]]      <- ts$time
    cohort[[paste0("status_", nm)]]    <- ts$status
    cohort[[paste0("time_", nm, "_10")]]   <- ts10$time_10
    cohort[[paste0("status_", nm, "_10")]] <- ts10$status_10
  }

  log(sprintf(
    "Outcome events (within %d-year horizon, CKM cohort, n=%d): MACE=%d, ASCVD=%d, HF=%d, CV_DEATH=%d.",
    cfg$T0_YEARS, nrow(cohort),
    sum(cohort$status_mace_10 == 1, na.rm = TRUE),
    sum(cohort$status_ascvd_10 == 1, na.rm = TRUE),
    sum(cohort$status_hf_10 == 1, na.rm = TRUE),
    sum(cohort$status_cvdeath_10 == 1, na.rm = TRUE)
  ))

  attr(cohort, "baseline_covariates") <- c(
    "age", "sex", "egfr", "systolic_bp",
    "smoking_habit", "drinking_habit", "exercise_habit",
    "hdl_cholesterol", "ldl_cholesterol"
  )
  cohort
}
