# Step 2 (run before 02_build_cohort, see run_pipeline.R): claims-based
# HTN_MED / DM_MED flags (used by MET4 and the CKM stage 2/3 rule in
# 02_build_cohort.R), derived from GNL_NM_CD -> ATC code during each
# member's lookback period. Ported unchanged from CKM_Drug/Korea/
# NHIS_SAMPLE's 04_medication_flags.R.
#
# Assumes push_atc_map() has already been called on `con` (done once in
# run_pipeline.R, since 06_drug_exposure.R needs the same "#atc_map" temp
# table and there's no reason to push it twice).

library(DBI)
library(dplyr)
library(dbplyr)

add_medication_flags <- function(con, cohort_base, cfg, log = message) {

  windows <- cohort_base %>% select(PERSON_ID, lookback_start, index_date)
  push_window_table(con, windows, name = "#med_windows")
  win_tbl <- tbl(con, "#med_windows")

  rx <- prescription_claims(con, cfg) %>%
    inner_join(win_tbl, by = "PERSON_ID") %>%
    filter(rx_date >= lookback_start, rx_date < index_date)

  htn_flag <- rx %>%
    filter(atc3 %in% cfg$htn_med_atc_prefix) %>%
    distinct(PERSON_ID) %>%
    collect() %>%
    mutate(htn_med = TRUE)

  dm_flag <- rx %>%
    filter(atc3 %in% cfg$dm_med_atc_prefix) %>%
    distinct(PERSON_ID) %>%
    collect() %>%
    mutate(dm_med = TRUE)

  out <- cohort_base %>%
    left_join(htn_flag, by = "PERSON_ID") %>%
    left_join(dm_flag, by = "PERSON_ID") %>%
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
