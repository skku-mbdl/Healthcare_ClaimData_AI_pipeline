# ==============================================================================
# 02_extract_drug_exposure.R -- baseline drug-exposure extraction (Japan,
# long format).
#
# Pulls, for the same CKM cohort as 01_extract_cohort_outcomes.R, every WHO
# ATC code with PDC > 0.8 during the 1-year lookback before index_date, from
# CKM_DRUG.dbo.drug_exposure (already computed by this project's own Stage
# 1 -- stage1_drug_screen/scripts/06_drug_exposure.R), restricted to the
# candidate-drug allow-list loaded automatically via core/R/
# candidate_drugs.R (ATC codes that already passed Stage 1's single-drug Cox
# screen).
#
# ARCHITECTURE CHANGE from the original CKM_PREVENT/JAPAN/SideEffect_Model:
# that project queried JMDC.dbo.Drug directly and recomputed its own
# DATEDIFF-based look-back exposure, entirely separately from Stage 1's own
# PDC calculation -- a real duplicate-logic / drift risk (two independently
# maintained definitions of "exposed"). Now that Stage 2 reads Stage 1's own
# drug_exposure table directly, exposure is computed in exactly ONE place
# (Stage 1's 06_drug_exposure.R) for both the Cox screen and the side-effect
# model, and this file only filters/aggregates it -- mirroring how Korea's
# Stage 2 always worked.
#
# Runs the ATC-allowlist + MIN_DRUG_PREVALENCE filter SQL-side when a live
# CKM_DRUG connection is available; falls back to the equivalent computation
# in R over Stage 1's file checkpoints when no connection is available.
# ==============================================================================

candidate_drugs <- load_candidate_drugs(con, stage1_results_root, table_name = "cox_results", log = progress)
atc_codes <- candidate_atc_union(candidate_drugs)

progress(sprintf(
  "Extracting baseline PDC>0.8 exposures (ATC allow-list of %d code(s) from Stage 1's Cox screen, >= %.1f%% prevalence cutoff)...",
  length(atc_codes), MIN_DRUG_PREVALENCE * 100
))

if (!is.null(con)) {
  atc_filter_sql <- atc_allowlist_sql_in(atc_codes, "de.who_atc_code")

  drug_exposure_query <- sprintf("
    WITH cohort AS (
      SELECT member_id FROM CKM_DRUG.dbo.cohort_processed
    ),
    cohort_n AS (
      SELECT COUNT(*) AS n_cohort FROM cohort
    ),
    exposures AS (
      SELECT DISTINCT de.member_id, de.who_atc_code
      FROM CKM_DRUG.dbo.drug_exposure de
      INNER JOIN cohort c ON c.member_id = de.member_id
      WHERE de.exposed = 1
        AND %s
    ),
    eligible_drugs AS (
      SELECT who_atc_code
      FROM exposures
      GROUP BY who_atc_code
      HAVING CAST(COUNT(DISTINCT member_id) AS float) / (SELECT n_cohort FROM cohort_n) >= %f
    )
    SELECT e.member_id, e.who_atc_code
    FROM exposures e
    INNER JOIN eligible_drugs ed ON ed.who_atc_code = e.who_atc_code
  ", atc_filter_sql, MIN_DRUG_PREVALENCE)

  drug_exposure_long <- dbGetQuery(con, drug_exposure_query)
} else {
  # File-checkpoint fallback: same filter/prevalence logic, computed in R
  # over Stage 1's full drug_exposure checkpoint instead of SQL-side.
  drug_exposure_raw <- read_table_or_checkpoint(NULL, "drug_exposure", stage1_results_root, log = progress)
  n_cohort <- length(unique(cohort_raw$patient_id))

  exposures <- drug_exposure_raw %>%
    filter(exposed == 1, who_atc_code %in% atc_codes, member_id %in% cohort_raw$patient_id) %>%
    distinct(member_id, who_atc_code)

  eligible_drugs <- exposures %>%
    count(who_atc_code, name = "n_exposed") %>%
    filter(n_exposed / n_cohort >= MIN_DRUG_PREVALENCE) %>%
    pull(who_atc_code)

  drug_exposure_long <- exposures %>% filter(who_atc_code %in% eligible_drugs)
}

progress(sprintf(
  "Pulled %d (patient, drug) exposure rows across %d distinct eligible ATC codes.",
  nrow(drug_exposure_long), length(unique(drug_exposure_long$who_atc_code))
))

if (nrow(drug_exposure_long) == 0) {
  stop("No drug exposures met the ATC allow-list + MIN_DRUG_PREVALENCE cutoff (",
       MIN_DRUG_PREVALENCE, ") -- check MIN_DRUG_PREVALENCE in 00_config.R, that ",
       "Stage 1's cox_results candidate codes actually match drug_exposure.who_atc_code, ",
       "or that drug_exposure has exposed = 1 rows for this cohort.")
}

drug_exposure_long <- drug_exposure_long %>%
  dplyr::rename(patient_id = member_id, atc_code = who_atc_code)

# ATC-code -> readable-name lookup for labeling Boruta/Cox output later. No
# per-drug brand/general name applies at this granularity (predictors are
# ATC codes, not individual drug products), so the ATC code doubles as its
# own label.
drug_name_lookup <- drug_exposure_long %>%
  distinct(atc_code) %>%
  mutate(drug_label = atc_code)
