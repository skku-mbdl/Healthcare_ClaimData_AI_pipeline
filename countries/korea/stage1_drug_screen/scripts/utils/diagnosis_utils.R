# Shared helpers for querying 상병내역(40t, diagnosis detail) by ICD-10
# category code, reused by washout/CKM-exclusion (02), comorbidities (03),
# and outcomes (05). Ported unchanged from CKM_Drug/Korea/NHIS_SAMPLE's
# scripts/utils/diagnosis_utils.R. Uses raw_tbl()/union_years_tbl() from
# core/R/db_common.R (must be sourced first).
#
# By instruction, this pipeline only reads T1 (의과/보건기관) claims tables
# -- T2 (치과/한방) and T3 (약국) are never queried, for GY20/GY30/GY40/GY60
# alike. See cfg$tbl in utils/config.R.
#
# Structural note: PERSON_ID and FORM_CD live only on the 명세서(20t)
# claim-header table, not on 40t (which only has KEY_SEQ + SEQ_NO). Every
# diagnosis lookup therefore joins 40t -> 20t via KEY_SEQ. That join is done
# PER YEAR -- i.e. NHID_GY40_T1_2010 joined only to NHID_GY20_T1_2010 --
# before the per-year results are UNION ALL'd together. This avoids relying
# on KEY_SEQ being globally unique across years, which the data dictionary
# doesn't guarantee (it only documents KEY_SEQ as "year + serial number").
#
# RECU_FR_DT is an 8-character char(8) field, assumed YYYYMMDD; parsed with
# TRY_CONVERT(date, ..., 112) (SQL Server style code 112 = ISO yyyymmdd) so
# a wrong assumption fails a parse (NULL) instead of silently misparsing.

library(DBI)
library(dplyr)
library(dbplyr)

#' Diagnosis claims (with PERSON_ID + FORM_CD attached) across every year in
#' cfg$study_years, already joined 40t -> 20t (T1 only).
diagnosis_claims_all_years <- function(con, cfg) {
  parts <- lapply(cfg$study_years, function(y) {
    t20 <- raw_tbl(con, cfg$tbl$claim20, y)
    t40 <- raw_tbl(con, cfg$tbl$claim40, y)
    if (is.null(t20) || is.null(t40)) return(NULL)

    t20 <- t20 %>% select(PERSON_ID, KEY_SEQ, FORM_CD)
    t40 %>%
      inner_join(t20, by = "KEY_SEQ") %>%
      mutate(
        dx_date = sql("TRY_CONVERT(date, RECU_FR_DT, 112)"),
        icd3 = substr(SICK_SYM, 1, 3)
      ) %>%
      select(PERSON_ID, KEY_SEQ, FORM_CD, dx_date, icd3)
  })
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) == 0) {
    stop(sprintf("No usable %s/%s tables found for years %s.",
                  cfg$tbl$claim20, cfg$tbl$claim40, paste(range(cfg$study_years), collapse = "-")), call. = FALSE)
  }
  Reduce(dplyr::union_all, parts)
}

#' Diagnosis claims matching a set of 3-character ICD-10 category codes.
#' Optionally restricted to a set of FORM_CD values (e.g. cfg$admission_form_cd).
diagnosis_claims <- function(con, cfg, icd10_codes, form_cd = NULL) {
  dx <- diagnosis_claims_all_years(con, cfg) %>%
    filter(icd3 %in% icd10_codes)

  if (!is.null(form_cd)) {
    dx <- dx %>% filter(FORM_CD %in% form_cd)
  }

  dx %>% select(PERSON_ID, KEY_SEQ, FORM_CD, dx_date)
}

#' Any death on record: NHID_JK.DTH_YM non-missing, regardless of cause.
#' death_month_start is the 1st of DTH_YM's month (DTH_YM is YYYYMM,
#' month-only granularity); 05_outcomes.R rolls this to month-end in R (via
#' lubridate) as the death date proxy.
all_death_claims <- function(con, cfg) {
  elig <- union_years_tbl(con, cfg$tbl$eligibility, cfg$study_years)

  elig %>%
    filter(!is.na(DTH_YM), DTH_YM != "") %>%
    mutate(death_month_start = sql("TRY_CONVERT(date, DTH_YM + '01', 112)")) %>%
    distinct(PERSON_ID, death_month_start) %>%
    select(PERSON_ID, death_month_start)
}

#' Cardiovascular death: NHID_JK.DTH_YM non-missing, with DTH_CODE1 or
#' DTH_CODE2 starting with "I" (any circulatory-chapter cause of death).
cv_death_claims <- function(con, cfg) {
  all_death_claims(con, cfg) %>%
    inner_join(
      union_years_tbl(con, cfg$tbl$eligibility, cfg$study_years) %>%
        filter((substr(DTH_CODE1, 1, 1) == "I") | (substr(DTH_CODE2, 1, 1) == "I")) %>%
        distinct(PERSON_ID),
      by = "PERSON_ID"
    )
}
