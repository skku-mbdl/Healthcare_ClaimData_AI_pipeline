# Step 4: comorbidity flags (Hypertension, Diabetes, Dyslipidemia, CKD).
# Ported unchanged from CKM_Drug/Korea/NHIS_SAMPLE's 03_covariates.R.
#
# Age, sex, eGFR, systolic BP, smoking/drinking/exercise habit, HDL and LDL
# cholesterol are already present on the cohort table from the index
# checkup (01_load_cohort_base.R) and eGFR calc (02_build_cohort.R); this
# step only adds the diagnosis-based comorbidity flags.
#
# Rule (protocol): a comorbidity is counted if the member has >=2 outpatient
# claims OR >=1 inpatient claim with a matching ICD-10 code, ascertained
# during the lookback period (same window used for medication flags and
# drug exposure). "Claim" is counted as a distinct KEY_SEQ (claim-header
# identifier), matching one 명세서/20t record.

library(dplyr)
library(dbplyr)

comorbidity_flag_members <- function(con, icd10_codes, cfg) {
  win_tbl <- tbl(con, "#comorbidity_windows")

  dx_win <- diagnosis_claims(con, cfg, icd10_codes) %>%
    inner_join(win_tbl, by = "PERSON_ID") %>%
    filter(dx_date >= lookback_start, dx_date < index_date)

  outpatient_n <- dx_win %>%
    filter(FORM_CD %in% cfg$outpatient_form_cd) %>%
    distinct(PERSON_ID, KEY_SEQ) %>%
    count(PERSON_ID, name = "n_outpatient")

  admission_n <- dx_win %>%
    filter(FORM_CD %in% cfg$admission_form_cd) %>%
    distinct(PERSON_ID, KEY_SEQ) %>%
    count(PERSON_ID, name = "n_admission")

  outpatient_n %>%
    full_join(admission_n, by = "PERSON_ID") %>%
    collect() %>%
    mutate(
      n_outpatient = coalesce(n_outpatient, 0L),
      n_admission = coalesce(n_admission, 0L)
    ) %>%
    filter(n_outpatient >= 2 | n_admission >= 1) %>%
    pull(PERSON_ID)
}

add_covariates <- function(con, cohort, cfg, log = message) {

  windows <- cohort %>% select(PERSON_ID, lookback_start, index_date)
  push_window_table(con, windows, name = "#comorbidity_windows")

  comorbidities <- names(cfg$comorbidity_icd10)
  flags <- lapply(comorbidities, function(name) {
    comorbidity_flag_members(con, cfg$comorbidity_icd10[[name]], cfg)
  })
  names(flags) <- comorbidities

  cohort <- cohort %>%
    mutate(
      comorbid_hypertension = PERSON_ID %in% flags$hypertension,
      comorbid_diabetes = PERSON_ID %in% flags$diabetes,
      comorbid_dyslipidemia = PERSON_ID %in% flags$dyslipidemia,
      comorbid_ckd = PERSON_ID %in% flags$ckd
    )

  log(sprintf(
    "Comorbidities: HTN %s, DM %s, dyslipidemia %s, CKD %s (of %s members).",
    sum(cohort$comorbid_hypertension), sum(cohort$comorbid_diabetes),
    sum(cohort$comorbid_dyslipidemia), sum(cohort$comorbid_ckd), nrow(cohort)
  ))

  cohort
}
