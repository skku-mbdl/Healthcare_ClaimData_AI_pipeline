# Step 4: claims-based HTN_MED / DM_MED flags (used by MET4 and the CKM
# stage 2/3 rule in 02_build_cohort.R), derived from Drug + Drug_master WHO
# ATC codes during each member's lookback period. Ported unchanged from
# CKM_Drug/Japan/scripts/04_medication_flags.R.
#
# date_of_prescription is a Character field in the raw schema; parsed with
# TRY_CONVERT the same way as date_of_health_checkup in
# 01_load_cohort_base.R.

library(DBI)
library(dplyr)
library(dbplyr)

add_medication_flags <- function(con, cohort_base, cfg, log = message) {

  windows <- cohort_base %>% select(member_id, lookback_start, index_date)
  push_window_table(con, windows, name = "#med_windows")
  win_tbl <- tbl(con, "#med_windows")

  drug_atc <- tbl(con, "Drug") %>%
    inner_join(tbl(con, "Drug_master") %>% select(jmdc_drug_code, who_atc_code),
               by = "jmdc_drug_code") %>%
    mutate(
      rx_date = sql("TRY_CONVERT(date, date_of_prescription)"),
      atc3 = substr(who_atc_code, 1, 3)
    ) %>%
    inner_join(win_tbl, by = "member_id") %>%
    filter(!is.na(rx_date), rx_date >= lookback_start, rx_date < index_date)

  htn_flag <- drug_atc %>%
    filter(atc3 %in% cfg$htn_med_atc_prefix) %>%
    distinct(member_id) %>%
    collect() %>%
    mutate(htn_med = TRUE)

  dm_flag <- drug_atc %>%
    filter(atc3 %in% cfg$dm_med_atc_prefix) %>%
    distinct(member_id) %>%
    collect() %>%
    mutate(dm_med = TRUE)

  out <- cohort_base %>%
    left_join(htn_flag, by = "member_id") %>%
    left_join(dm_flag, by = "member_id") %>%
    mutate(
      htn_med = coalesce(htn_med, FALSE),
      dm_med = coalesce(dm_med, FALSE)
    )

  log(sprintf(
    "Medication flags: %s on antihypertensives, %s on antidiabetics (lookback period).",
    sum(out$htn_med), sum(out$dm_med)
  ))
  out
}
