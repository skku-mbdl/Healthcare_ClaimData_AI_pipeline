# Step 5: identify each member's follow-up end and, per primary outcome
# (CHD, Stroke, Heart Failure, Cardiovascular death), the first post-index
# occurrence within the observation period. Ported unchanged from
# CKM_Drug/Japan/scripts/05_outcomes.R.
#
# Observation period: index date until death, event occurrence, or
# 2022-07-31 (design doc). Enrollment.observation_end / withdrawal_death_flag
# give the last observable month and whether that end was death-driven;
# both are assumed to be YYYYMM integers, same format as birth_yyyymm in
# 01_load_cohort_base.R.
#
# CHD/Stroke/HF use admission claims only (Inpatient + DPC, per
# cfg$admission_claim_types). Cardiovascular death is not code-restricted
# to admission claims -- it is any Diagnosis record with
# outcome_death_flag == 1 under a circulatory-chapter diagnosis (see
# cv_death_claims() in utils/diagnosis_utils.R).

library(dplyr)
library(dbplyr)
library(lubridate)

month_end_date <- function(yyyymm) {
  start <- ymd(sprintf("%d01", yyyymm))
  ceiling_date(start, "month") - days(1)
}

first_event_after_index <- function(con, icd10_codes, claim_types, win_tbl_name) {
  win_tbl <- tbl(con, win_tbl_name)
  diagnosis_claims(con, icd10_codes, claim_types = claim_types) %>%
    inner_join(win_tbl, by = "member_id") %>%
    filter(dx_date > index_date, dx_date <= followup_end) %>%
    group_by(member_id) %>%
    summarise(event_date = min(dx_date, na.rm = TRUE)) %>%
    ungroup() %>%
    collect()
}

first_cv_death_after_index <- function(con, win_tbl_name) {
  win_tbl <- tbl(con, win_tbl_name)
  cv_death_claims(con) %>%
    inner_join(win_tbl, by = "member_id") %>%
    filter(dx_date > index_date, dx_date <= followup_end) %>%
    group_by(member_id) %>%
    summarise(event_date = min(dx_date, na.rm = TRUE)) %>%
    ungroup() %>%
    collect()
}

add_outcomes <- function(con, cohort, cfg, log = message) {

  cohort <- cohort %>%
    mutate(
      raw_observation_end = month_end_date(observation_end),
      followup_end = pmin(raw_observation_end, cfg$observation_end_cutoff, na.rm = TRUE),
      died_in_data = !is.na(withdrawal_death_flag) & withdrawal_death_flag == 1
    )

  windows <- cohort %>% select(member_id, index_date, followup_end)
  push_window_table(con, windows, name = "#outcome_windows")

  outcome_events <- list(
    stroke = first_event_after_index(con, cfg$outcome_icd10$stroke, cfg$admission_claim_types, "#outcome_windows"),
    heart_failure = first_event_after_index(con, cfg$outcome_icd10$heart_failure, cfg$admission_claim_types, "#outcome_windows"),
    chd = first_event_after_index(con, cfg$outcome_icd10$chd, cfg$admission_claim_types, "#outcome_windows"),
    cv_death = first_cv_death_after_index(con, "#outcome_windows")
  )

  for (name in names(outcome_events)) {
    ev <- outcome_events[[name]] %>% rename(!!paste0(name, "_date") := event_date)
    cohort <- cohort %>% left_join(ev, by = "member_id")
    event_col <- paste0(name, "_event")
    date_col <- paste0(name, "_date")
    cohort[[event_col]] <- !is.na(cohort[[date_col]])
    # dplyr::if_else(), not base ifelse() -- base ifelse() silently drops the
    # Date class from its yes/no arguments (returns the raw numeric day-count
    # instead), which then fails with "can only subtract from Date objects"
    # on the next line since cohort$index_date is a proper Date.
    cohort[[paste0(name, "_time")]] <- as.numeric(
      dplyr::if_else(cohort[[event_col]], cohort[[date_col]], cohort$followup_end) - cohort$index_date
    )
  }

  log(sprintf(
    "Outcomes during follow-up: stroke %s, heart failure %s, CHD %s, CV death %s (of %s members).",
    sum(cohort$stroke_event), sum(cohort$heart_failure_event),
    sum(cohort$chd_event), sum(cohort$cv_death_event), nrow(cohort)
  ))

  cohort
}
