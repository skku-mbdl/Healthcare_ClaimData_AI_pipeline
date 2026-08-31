# Step 6: PDC (Proportion of Days Covered) per ATC code during each
# member's lookback period; exposure = PDC > cfg$pdc_threshold. Ported
# unchanged from CKM_Drug/Korea/NHIS_SAMPLE's 06_drug_exposure.R.
#
# Exposure is evaluated per ATC code (GNL_NM_CD -> NHIS_ATC_MAPPED.xlsx),
# one drug at a time -- per the protocol's "Exposure" and "Task Rule"
# sections. rx_date (RECU_FR_DT, see utils/drug_utils.R) is the supply
# start date; MDCN_EXEC_FREQ (총투여일수, total days of medication for that
# claim line) is the days-of-supply length. Overlapping or duplicate
# prescriptions for the same member+ATC are merged via interval union so
# covered days aren't double-counted.
#
# Assumes push_atc_map() has already been called on `con` (see
# run_pipeline.R and the note at the top of 04_medication_flags.R).

library(dplyr)
library(dbplyr)

#' Merge possibly-overlapping/adjacent [start, end] day intervals and
#' return total unique days covered. Intervals are treated as touching
#' (mergeable) when one starts the day right after another ends.
merge_interval_days <- function(starts, ends) {
  ord <- order(starts)
  starts <- starts[ord]
  ends <- ends[ord]
  covered <- 0L
  cur_start <- starts[1]
  cur_end <- ends[1]
  for (i in seq_along(starts)[-1]) {
    if (starts[i] <= cur_end + 1) {
      cur_end <- max(cur_end, ends[i])
    } else {
      covered <- covered + as.integer(cur_end - cur_start) + 1L
      cur_start <- starts[i]
      cur_end <- ends[i]
    }
  }
  covered + as.integer(cur_end - cur_start) + 1L
}

add_drug_exposure <- function(con, cohort, cfg, log = message) {

  windows <- cohort %>% select(PERSON_ID, lookback_start, index_date)
  push_window_table(con, windows, name = "#exposure_windows")
  win_tbl <- tbl(con, "#exposure_windows")

  prescriptions <- prescription_claims(con, cfg) %>%
    inner_join(win_tbl, by = "PERSON_ID") %>%
    # Clip each prescription's supply interval to the lookback window.
    # Written as explicit CASE WHEN rather than pmin()/pmax() because
    # GREATEST/LEAST are only native to SQL Server 2022+; this avoids
    # depending on dbplyr's translation picking a compatible fallback.
    mutate(
      clipped_start = sql(
        "CASE WHEN rx_date > lookback_start THEN rx_date ELSE lookback_start END"
      ),
      clipped_end = sql(
        "CASE WHEN DATEADD(day, MDCN_EXEC_FREQ - 1, rx_date) < DATEADD(day, -1, index_date)
              THEN DATEADD(day, MDCN_EXEC_FREQ - 1, rx_date)
              ELSE DATEADD(day, -1, index_date) END"
      )
    ) %>%
    filter(clipped_start <= clipped_end) %>%
    select(PERSON_ID, ATC_CODE, clipped_start, clipped_end) %>%
    collect()

  log(sprintf(
    "Drug exposure: %s prescription intervals across %s distinct ATC codes in lookback windows.",
    nrow(prescriptions), n_distinct(prescriptions$ATC_CODE)
  ))

  pdc <- prescriptions %>%
    group_by(PERSON_ID, ATC_CODE) %>%
    summarise(
      days_covered = merge_interval_days(clipped_start, clipped_end),
      .groups = "drop"
    ) %>%
    left_join(cohort %>% select(PERSON_ID, index_date, lookback_start), by = "PERSON_ID") %>%
    mutate(
      lookback_len = as.numeric(index_date - lookback_start),
      pdc = days_covered / lookback_len,
      exposed = pdc > cfg$pdc_threshold
    ) %>%
    select(PERSON_ID, ATC_CODE, days_covered, pdc, exposed)

  n_atc <- n_distinct(pdc$ATC_CODE)
  log(sprintf(
    "PDC computed for %s member-ATC pairs across %s ATC codes; %s exposures with PDC > %.2f.",
    nrow(pdc), n_atc, sum(pdc$exposed), cfg$pdc_threshold
  ))

  pdc
}
