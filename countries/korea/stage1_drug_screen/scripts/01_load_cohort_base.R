# Step 1: identify each member's index year (first health-checkup year on
# record in [checkup_year_min, checkup_year_max]) and pull the checkup
# values + demographics needed to build the CKM cohort flag and covariates.
# Ported unchanged from CKM_Drug/Korea/NHIS_SAMPLE's 01_load_cohort_base.R.
#
# The checkup table (NHID_GJ_<year>) has no day/month, only HCHK_YEAR --
# index date is approximated as cfg$index_date_month_day of that year. Age/
# sex aren't on the checkup table at all; they come from the eligibility
# table (NHID_JK_<year>) for the SAME year as the index checkup (STND_Y ==
# HCHK_YEAR), since AGE_GROUP is "age as of STND_Y" per the data dictionary.
#
# Both NHID_GJ and NHID_JK are unioned across years with union_years_tbl()
# (not the per-(type,year)-join pattern in diagnosis_utils.R/drug_utils.R --
# both tables already carry PERSON_ID directly, so no KEY_SEQ join is
# needed here, and no cross-year collision risk applies to a simple UNION).
#
# age_group_to_age() is Korea-specific (NHID_JK.AGE_GROUP is a categorical
# 5-year band -- the NHIS Sample Cohort's public-use eligibility table
# doesn't expose exact birth year, unlike Japan's Enrollment table), so it
# lives here rather than in core/R/clinical_calcs.R.

library(DBI)
library(dplyr)
library(dbplyr)

#' NHID_JK.AGE_GROUP is a categorical 5-year band (0 = "0세", 1 = "1~4세",
#' 2 = "5~9세", ..., 17 = "80~84세", 18 = "85세+"). This approximates a
#' single representative numeric age per band (the band midpoint; the
#' open-ended top band is approximated as 87), for use in the eGFR formula
#' and as a Cox covariate -- a real resolution limit of the data, not a
#' pipeline bug.
age_group_to_age <- function(age_group) {
  midpoints <- c(0, 2.5, 7, 12, 17, 22, 27, 32, 37, 42,
                 47, 52, 57, 62, 67, 72, 77, 82, 87)
  midpoints[as.integer(age_group) + 1L]
}

load_cohort_base <- function(con, cfg, log = message) {

  checkup <- union_years_tbl(con, cfg$tbl$checkup, cfg$checkup_year_min:cfg$checkup_year_max) %>%
    mutate(hchk_year_int = as.integer(HCHK_YEAR)) %>%
    filter(hchk_year_int >= cfg$checkup_year_min, hchk_year_int <= cfg$checkup_year_max)

  n_checkup <- checkup %>% tally() %>% pull(n)
  log(sprintf(
    "NHID_GJ_%s..%s: %s checkup rows in the index-year window.",
    cfg$checkup_year_min, cfg$checkup_year_max, n_checkup
  ))

  # First checkup year per member = index year.
  index_checkup <- checkup %>%
    group_by(PERSON_ID) %>%
    window_order(hchk_year_int) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    select(
      PERSON_ID,
      index_year = hchk_year_int,
      waist = WAIST,
      hdl_cholesterol = HDL_CHOLE,
      ldl_cholesterol = LDL_CHOLE,
      triglyceride = TRIGLYCERIDE,
      systolic_bp = BP_HIGH,
      diastolic_bp = BP_LWST,
      fasting_blood_sugar = BLDS,
      serum_creatinine = CREATININE,
      uric_protein_qualitative = OLIG_PROTE_CD,
      smoking_habit = SMK_STAT_TYPE_RSPS_CD,
      drinking_habit = DRNK_HABIT_RSPS_CD,
      exercise_habit = MOV30_WEK_FREQ_ID
    )

  eligibility <- union_years_tbl(con, cfg$tbl$eligibility, cfg$checkup_year_min:cfg$checkup_year_max) %>%
    mutate(stnd_y_int = as.integer(STND_Y)) %>%
    select(PERSON_ID, stnd_y_int, sex = SEX, age_group = AGE_GROUP)

  cohort_base <- index_checkup %>%
    inner_join(eligibility, by = c("PERSON_ID" = "PERSON_ID", "index_year" = "stnd_y_int")) %>%
    collect() %>%
    mutate(
      age = age_group_to_age(age_group),
      sex = as.character(sex),
      # OLIG_PROTE_CD is a char(1) column in the raw schema ("1".."6");
      # cast explicitly so the %in% comparison against
      # cfg$ckm_urine_protein_grades (numeric) in 02_build_cohort.R isn't
      # relying on %in%'s implicit character coercion.
      uric_protein_qualitative = as.integer(uric_protein_qualitative),
      index_date = as.Date(sprintf("%d-%s", index_year, cfg$index_date_month_day)),
      # Lookback period: 1 year before index date, up to (excluding) index
      # date itself. Every later step (comorbidities, medication flags, PDC)
      # filters claims into [lookback_start, index_date).
      lookback_start = index_date - cfg$lookback_days
    )

  n_index_years <- index_checkup %>% tally() %>% pull(n)
  n_unmatched_eligibility <- n_index_years - nrow(cohort_base)
  if (n_unmatched_eligibility > 0) {
    log(sprintf(
      "WARNING: %s members had an index checkup but no matching NHID_JK record for that same year -- dropped.",
      n_unmatched_eligibility
    ))
  }

  log(sprintf("Base cohort: %s members with an index year and matched demographics.", nrow(cohort_base)))
  cohort_base
}
