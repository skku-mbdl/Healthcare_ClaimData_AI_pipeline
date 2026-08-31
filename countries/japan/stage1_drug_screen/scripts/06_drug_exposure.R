# Step 6: PDC (Proportion of Days Covered) per WHO ATC code during each
# member's lookback period; exposure = PDC > cfg$pdc_threshold. Ported
# unchanged from CKM_Drug/Japan/scripts/06_drug_exposure.R.
#
# Exposure is evaluated per Drug_master.who_atc_code (confirmed with study
# lead -- not per jmdc_drug_code, not aggregated by ingredient/general_name).
# date_of_dispense is used as the supply start date (falls back to
# date_of_prescription when dispense date is missing); administered_days is
# the days-of-supply length. Overlapping or duplicate prescriptions for the
# same member+ATC are merged via interval union so covered days aren't
# double-counted.

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

  windows <- cohort %>% select(member_id, lookback_start, index_date)
  push_window_table(con, windows, name = "#exposure_windows")
  win_tbl <- tbl(con, "#exposure_windows")

  prescriptions <- tbl(con, "Drug") %>%
    inner_join(
      tbl(con, "Drug_master") %>% select(jmdc_drug_code, who_atc_code),
      by = "jmdc_drug_code"
    ) %>%
    mutate(
      dispense_date = sql(
        "TRY_CONVERT(date, COALESCE(date_of_dispense, date_of_prescription))"
      )
    ) %>%
    filter(!is.na(dispense_date), !is.na(administered_days), administered_days > 0) %>%
    inner_join(win_tbl, by = "member_id") %>%
    # Clip each prescription's supply interval to the lookback window.
    # Written as explicit CASE WHEN rather than pmin()/pmax() because
    # GREATEST/LEAST are only native to SQL Server 2022+; this avoids
    # depending on dbplyr's translation picking a compatible fallback.
    mutate(
      clipped_start = sql(
        "CASE WHEN dispense_date > lookback_start THEN dispense_date ELSE lookback_start END"
      ),
      clipped_end = sql(
        "CASE WHEN DATEADD(day, administered_days - 1, dispense_date) < DATEADD(day, -1, index_date)
              THEN DATEADD(day, administered_days - 1, dispense_date)
              ELSE DATEADD(day, -1, index_date) END"
      )
    ) %>%
    filter(clipped_start <= clipped_end) %>%
    select(member_id, who_atc_code, clipped_start, clipped_end) %>%
    collect()

  log(sprintf(
    "Drug exposure: %s prescription intervals across %s distinct ATC codes in lookback windows.",
    nrow(prescriptions), n_distinct(prescriptions$who_atc_code)
  ))

  pdc <- prescriptions %>%
    group_by(member_id, who_atc_code) %>%
    summarise(
      days_covered = merge_interval_days(clipped_start, clipped_end),
      .groups = "drop"
    ) %>%
    left_join(cohort %>% select(member_id, index_date, lookback_start), by = "member_id") %>%
    mutate(
      lookback_len = as.numeric(index_date - lookback_start),
      pdc = days_covered / lookback_len,
      exposed = pdc > cfg$pdc_threshold
    ) %>%
    select(member_id, who_atc_code, days_covered, pdc, exposed)

  n_atc <- n_distinct(pdc$who_atc_code)
  log(sprintf(
    "PDC computed for %s member-ATC pairs across %s ATC codes; %s exposures with PDC > %.2f.",
    nrow(pdc), n_atc, sum(pdc$exposed), cfg$pdc_threshold
  ))

  pdc
}
