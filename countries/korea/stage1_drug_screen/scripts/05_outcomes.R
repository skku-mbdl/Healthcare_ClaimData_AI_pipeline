# Step 5: identify each member's follow-up end and, per primary outcome
# (CHD, Stroke, Heart Failure, Cardiovascular death), the first post-index
# occurrence within the observation period. Ported unchanged from
# CKM_Drug/Korea/NHIS_SAMPLE's 05_outcomes.R.
#
# Observation period: index date until death, event occurrence, or
# 2013-12-31 (protocol). Death is only known at month granularity
# (NHID_JK.DTH_YM); the death date is approximated as the LAST day of that
# month (month_end_date() below), and a member's observation is censored at
# min(that date, cfg$observation_end_cutoff) if they died. Any death
# (whatever the cause) ends observation; only cardiovascular deaths (ICD
# starting "I") also count as the cv_death outcome event.
#
# CHD/Stroke/HF use admission claims only (cfg$admission_form_cd).
# Cardiovascular death is not FORM_CD-restricted -- it's sourced from the
# eligibility table's death fields, not a claim at all (see
# utils/diagnosis_utils.R's cv_death_claims()).

library(dplyr)
library(dbplyr)
library(lubridate)

month_end_date <- function(start) {
  ceiling_date(start, "month") - days(1)
}

first_event_after_index <- function(con, cfg, icd10_codes, form_cd, win_tbl_name) {
  win_tbl <- tbl(con, win_tbl_name)
  diagnosis_claims(con, cfg, icd10_codes, form_cd = form_cd) %>%
    inner_join(win_tbl, by = "PERSON_ID") %>%
    filter(dx_date > index_date, dx_date <= followup_end) %>%
    group_by(PERSON_ID) %>%
    summarise(event_date = min(dx_date, na.rm = TRUE)) %>%
    ungroup() %>%
    collect()
}

first_cv_death_after_index <- function(con, cfg, win_tbl_name) {
  win_tbl <- tbl(con, win_tbl_name)
  cv_death_claims(con, cfg) %>%
    inner_join(win_tbl, by = "PERSON_ID") %>%
    filter(death_month_start > index_date, death_month_start <= followup_end) %>%
    group_by(PERSON_ID) %>%
    summarise(event_month_start = min(death_month_start, na.rm = TRUE)) %>%
    ungroup() %>%
    collect()
}

add_outcomes <- function(con, cohort, cfg, log = message) {

  all_deaths <- all_death_claims(con, cfg) %>%
    collect() %>%
    mutate(death_date = month_end_date(death_month_start)) %>%
    group_by(PERSON_ID) %>%
    summarise(death_date = min(death_date, na.rm = TRUE), .groups = "drop")

  cohort <- cohort %>%
    left_join(all_deaths, by = "PERSON_ID") %>%
    mutate(
      died_in_data = !is.na(death_date),
      followup_end = pmin(
        if_else(died_in_data, death_date, as.Date(NA)),
        cfg$observation_end_cutoff,
        na.rm = TRUE
      )
    )

  windows <- cohort %>% select(PERSON_ID, index_date, followup_end)
  push_window_table(con, windows, name = "#outcome_windows")

  outcome_events <- list(
    stroke = first_event_after_index(con, cfg, cfg$outcome_icd10$stroke, cfg$admission_form_cd, "#outcome_windows"),
    heart_failure = first_event_after_index(con, cfg, cfg$outcome_icd10$heart_failure, cfg$admission_form_cd, "#outcome_windows"),
    chd = first_event_after_index(con, cfg, cfg$outcome_icd10$chd, cfg$admission_form_cd, "#outcome_windows")
  )
  cv_death_raw <- first_cv_death_after_index(con, cfg, "#outcome_windows") %>%
    rename(event_date = event_month_start) %>%
    mutate(event_date = month_end_date(event_date))
  outcome_events$cv_death <- cv_death_raw

  for (name in names(outcome_events)) {
    ev <- outcome_events[[name]] %>% rename(!!paste0(name, "_date") := event_date)
    cohort <- cohort %>% left_join(ev, by = "PERSON_ID")
    event_col <- paste0(name, "_event")
    date_col <- paste0(name, "_date")
    cohort[[event_col]] <- !is.na(cohort[[date_col]])
    # dplyr::if_else(), not base ifelse() -- base ifelse() silently drops
    # the Date class from its yes/no arguments, which then fails with "can
    # only subtract from Date objects" on the next line.
    cohort[[paste0(name, "_time")]] <- as.numeric(
      dplyr::if_else(cohort[[event_col]], cohort[[date_col]], cohort$followup_end) - cohort$index_date
    )
  }

  log(sprintf(
    "Outcomes during follow-up: stroke %s, heart failure %s, CHD %s, CV death %s (of %s members; %s died in data).",
    sum(cohort$stroke_event), sum(cohort$heart_failure_event),
    sum(cohort$chd_event), sum(cohort$cv_death_event), nrow(cohort), sum(cohort$died_in_data)
  ))

  cohort
}
