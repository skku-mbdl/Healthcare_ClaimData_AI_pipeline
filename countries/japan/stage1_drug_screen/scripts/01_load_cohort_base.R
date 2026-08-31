# Step 1: identify each member's index date (first annual health checkup on
# record between 2006 and 2021) and pull the checkup values + basic
# demographics needed to build the CKM cohort flag and covariates. Ported
# unchanged from CKM_Drug/Japan/scripts/01_load_cohort_base.R.
#
# Patient and Enrollment have an identical field list in the JMDC data
# dictionary; Enrollment is used here as the demographic/observation-window
# source.
#
# NOTE: date_of_health_checkup is stored as a Character field in the raw
# schema, format not fully verified against real data. TRY_CONVERT(date, ...)
# is used so SQL Server does the parsing; if `n_unparseable` logged below is
# non-trivial, inspect actual date string values before trusting the cohort.
#
# NOTE: the actual Enrollment column is month_and_year_of_birth_of_membe
# (truncated by exactly one character), not month_and_year_of_birth_of_member
# as spelled out in the JMDC data dictionary -- confirmed against the live
# schema via DBI::dbListFields(). Don't "fix" this to match the dictionary
# spelling; it's a real column-name-length limitation on the actual table.

library(DBI)
library(dplyr)
library(dbplyr)
library(lubridate)

load_cohort_base <- function(con, cfg, log = message) {

  enrollment <- tbl(con, "Enrollment") %>%
    transmute(
      member_id,
      birth_yyyymm = month_and_year_of_birth_of_membe,
      sex = gender_of_member,
      observation_start,
      observation_end,
      withdrawal_death_flag
    )

  checkup <- tbl(con, "Annual_health_checkup") %>%
    mutate(checkup_date = sql("TRY_CONVERT(date, date_of_health_checkup)")) %>%
    filter(
      !is.na(checkup_date),
      checkup_date >= cfg$index_date_min,
      checkup_date <= cfg$index_date_max
    )

  n_total <- tbl(con, "Annual_health_checkup") %>% tally() %>% pull(n)
  n_parsed <- checkup %>% tally() %>% pull(n)
  log(sprintf(
    "Annual_health_checkup: %s rows total, %s with a parseable date in [%s, %s]",
    n_total, n_parsed, cfg$index_date_min, cfg$index_date_max
  ))

  # First checkup per member = index date.
  index_checkup <- checkup %>%
    group_by(member_id) %>%
    window_order(checkup_date) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    select(
      member_id,
      index_date = checkup_date,
      abdominal_circumference,
      hdl_cholesterol,
      ldl_cholesterol,
      triglyceride,
      systolic_bp,
      diastolic_bp,
      fasting_blood_sugar,
      serum_creatinine,
      uric_protein_qualitative,
      smoking_habit,
      drinking_habit,
      exercise_habit
    )

  cohort_base <- enrollment %>%
    inner_join(index_checkup, by = "member_id") %>%
    collect() %>%
    mutate(
      # birth_yyyymm assumed to be a 6-digit YYYYMM integer (e.g. 198501);
      # adjust this parse if the real column turns out to be formatted
      # differently.
      birth_date = ymd(sprintf("%d01", birth_yyyymm)),
      age = as.integer(floor(time_length(interval(birth_date, index_date), "years"))),
      sex = as.character(sex),
      # Lookback period: 1 year before index date, up to (excluding) index
      # date itself. Every later step (comorbidities, medication flags,
      # PDC) filters claims into [lookback_start, index_date).
      lookback_start = index_date - cfg$lookback_days
    )

  log(sprintf("Base cohort: %s members with an index date.", nrow(cohort_base)))
  cohort_base
}
