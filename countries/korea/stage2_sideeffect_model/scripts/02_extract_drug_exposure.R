# ==============================================================================
# 02_extract_drug_exposure.R -- baseline drug-exposure extraction (Korea,
# long format).
#
# Pulls, for the same CKM==1 cohort as 01_extract_cohort_outcomes.R, every
# ATC code with PDC > 0.8 during the 1-year lookback before index_date, from
# CKM_DRUG.dbo.nhis_drug_exposure (already computed by this project's own
# Stage 1), restricted to the candidate-drug allow-list loaded automatically
# via core/R/candidate_drugs.R (ATC codes that already passed Stage 1's
# single-drug Cox screen) -- replacing the old manual candidate_drugs.csv
# copy-paste hand-off. Predictors are per-ATC-code adherence flags, not
# individual drug-product codes, since Korea's Stage 1 only ever computes
# exposure at ATC granularity.
#
# IMPORTANT: the NHIS Sample Cohort's raw prescription data has NO native
# ATC code -- 처방전교부상세내역(60t) only carries GNL_NM_CD. Every ATC_CODE
# this script ever sees was already produced upstream by Stage 1's
# prescription_claims() (stage1_drug_screen/scripts/utils/drug_utils.R),
# which INNER JOINs each GNL_NM_CD-coded prescription against
# NHIS_ATC_MAPPED.xlsx before the PDC/exposure calculation ever runs -- an
# inner join, so an unmapped GNL_NM_CD is silently invisible to this
# pipeline, not merely unclassified.
#
# Runs the ATC-allowlist + MIN_DRUG_PREVALENCE filter SQL-side when a live
# CKM_DRUG connection is available (keeps the wide matrix built in
# core/R/build_analysis_dataset.R tractable for Boruta without pulling a
# patient x hundreds-of-rare-ATC-codes result set over ODBC only to throw
# most of it away in R); falls back to the equivalent computation in R over
# Stage 1's file checkpoints when no connection is available.
# ==============================================================================

candidate_drugs <- load_candidate_drugs(con, stage1_results_root, table_name = "nhis_cox_results", log = progress)
atc_codes <- candidate_atc_union(candidate_drugs)

progress(sprintf(
  "Extracting baseline PDC>0.8 exposures (ATC allow-list of %d code(s) from Stage 1's Cox screen, >= %.1f%% prevalence cutoff)...",
  length(atc_codes), MIN_DRUG_PREVALENCE * 100
))

if (!is.null(con)) {
  atc_filter_sql <- atc_allowlist_sql_in(atc_codes, "de.ATC_CODE")

  drug_exposure_query <- sprintf("
    WITH cohort AS (
      SELECT PERSON_ID FROM CKM_DRUG.dbo.nhis_cohort_processed
    ),
    cohort_n AS (
      SELECT COUNT(*) AS n_cohort FROM cohort
    ),
    exposures AS (
      SELECT DISTINCT de.PERSON_ID, de.ATC_CODE
      FROM CKM_DRUG.dbo.nhis_drug_exposure de
      INNER JOIN cohort c ON c.PERSON_ID = de.PERSON_ID
      WHERE de.exposed = 1
        AND %s
    ),
    eligible_drugs AS (
      SELECT ATC_CODE
      FROM exposures
      GROUP BY ATC_CODE
      HAVING CAST(COUNT(DISTINCT PERSON_ID) AS float) / (SELECT n_cohort FROM cohort_n) >= %f
    )
    SELECT e.PERSON_ID, e.ATC_CODE
    FROM exposures e
    INNER JOIN eligible_drugs ed ON ed.ATC_CODE = e.ATC_CODE
  ", atc_filter_sql, MIN_DRUG_PREVALENCE)

  drug_exposure_long <- dbGetQuery(con, drug_exposure_query)
} else {
  # File-checkpoint fallback: same filter/prevalence logic, computed in R
  # over Stage 1's full nhis_drug_exposure checkpoint instead of SQL-side.
  drug_exposure_raw <- read_table_or_checkpoint(NULL, "nhis_drug_exposure", stage1_results_root,
                                                 checkpoint_name = "drug_exposure", log = progress)
  n_cohort <- length(unique(cohort_raw$patient_id))

  exposures <- drug_exposure_raw %>%
    filter(exposed == 1, ATC_CODE %in% atc_codes, PERSON_ID %in% cohort_raw$patient_id) %>%
    distinct(PERSON_ID, ATC_CODE)

  eligible_drugs <- exposures %>%
    count(ATC_CODE, name = "n_exposed") %>%
    filter(n_exposed / n_cohort >= MIN_DRUG_PREVALENCE) %>%
    pull(ATC_CODE)

  drug_exposure_long <- exposures %>% filter(ATC_CODE %in% eligible_drugs)
}

progress(sprintf(
  "Pulled %d (patient, drug) exposure rows across %d distinct eligible ATC codes.",
  nrow(drug_exposure_long), length(unique(drug_exposure_long$ATC_CODE))
))

if (nrow(drug_exposure_long) == 0) {
  stop("No drug exposures met the ATC allow-list + MIN_DRUG_PREVALENCE cutoff (",
       MIN_DRUG_PREVALENCE, ") -- check MIN_DRUG_PREVALENCE in 00_config.R, that ",
       "Stage 1's cox_results candidate codes actually match nhis_drug_exposure.ATC_CODE, ",
       "or that nhis_drug_exposure has exposed = 1 rows for this cohort.")
}

drug_exposure_long <- drug_exposure_long %>%
  dplyr::rename(patient_id = PERSON_ID, atc_code = ATC_CODE)

# ATC-code -> readable-name lookup for labeling Boruta/Cox output later. No
# per-drug brand/general name exists at this granularity (predictors are
# ATC codes, not individual drug products), so the ATC code doubles as its
# own label.
drug_name_lookup <- drug_exposure_long %>%
  distinct(atc_code) %>%
  mutate(drug_label = atc_code)
