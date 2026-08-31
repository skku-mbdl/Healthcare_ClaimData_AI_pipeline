# Shared helpers for querying Diagnosis + Diagnosis_master by ICD-10
# category code, reused by washout (02), comorbidities (03), and outcomes
# (05). Ported unchanged from CKM_Drug/Japan/scripts/utils/diagnosis_utils.R.
# ICD-10 codes in the design doc (e.g. "I21", "I50") are 3-character
# categories, matched here against Diagnosis_master.icd10_level3_code.
#
# date_of_medical_care_start is a Character field in the raw schema; parsed
# with TRY_CONVERT, same caveat noted in 01_load_cohort_base.R for
# date_of_health_checkup.

library(DBI)
library(dplyr)
library(dbplyr)

#' Diagnosis claims matching a set of 3-character ICD-10 category codes,
#' optionally restricted to specific claim types (e.g. admission-only).
diagnosis_claims <- function(con, icd10_level3_codes, claim_types = NULL) {
  dx <- tbl(con, "Diagnosis") %>%
    inner_join(
      tbl(con, "Diagnosis_master") %>%
        select(standard_disease_code, icd10_level3_code),
      by = "standard_disease_code"
    ) %>%
    mutate(
      dx_date = sql("TRY_CONVERT(date, date_of_medical_care_start)"),
      icd3 = substr(icd10_level3_code, 1, 3)
    ) %>%
    filter(icd3 %in% icd10_level3_codes)

  if (!is.null(claim_types)) {
    dx <- dx %>% filter(type_of_claim %in% claim_types)
  }

  dx %>% select(member_id, claim_id, type_of_claim, dx_date)
}

#' Cardiovascular death: a Diagnosis record with outcome_death_flag == 1
#' whose underlying diagnosis falls in the circulatory chapter (ICD-10
#' codes starting with "I"), per the design doc's death definition.
cv_death_claims <- function(con) {
  tbl(con, "Diagnosis") %>%
    filter(outcome_death_flag == 1) %>%
    inner_join(
      tbl(con, "Diagnosis_master") %>%
        select(standard_disease_code, icd10_level1_code),
      by = "standard_disease_code"
    ) %>%
    filter(substr(icd10_level1_code, 1, 1) == "I") %>%
    mutate(dx_date = sql("TRY_CONVERT(date, date_of_medical_care_start)")) %>%
    select(member_id, claim_id, dx_date)
}
