# Step 3: MET1-5 / MetS / CKM stage 2-3 flag, then washout exclusion, then
# the missing-covariate exclusion (protocol's Exclusion Criteria #1 and #2).
# Ported from CKM_Drug/Korea/NHIS_SAMPLE's 02_build_cohort.R -- the only
# change is calling core/R/clinical_calcs.R's egfr(creatinine, age, is_male)
# with an explicit is_male boolean (sex == "1") instead of the country's own
# copy of egfr() that took the raw sex code directly; see core/R/
# clinical_calcs.R's header for why.

library(dplyr)

build_cohort <- function(con, cohort_base, cfg, log = message) {

  cohort <- cohort_base %>%
    mutate(
      is_male = sex == "1",
      met1 = (is_male & waist >= cfg$met_waist_male) |
        (!is_male & waist >= cfg$met_waist_female),
      met2 = (is_male & hdl_cholesterol < cfg$met_hdl_male) |
        (!is_male & hdl_cholesterol < cfg$met_hdl_female),
      met3 = triglyceride >= cfg$met_tg,
      met4 = systolic_bp >= cfg$met_sbp | diastolic_bp >= cfg$met_dbp | htn_med,
      met5 = fasting_blood_sugar >= cfg$met_fbs,
      mets_count = rowSums(cbind(met1, met2, met3, met4, met5), na.rm = TRUE),
      mets = mets_count >= 3,
      egfr = egfr(serum_creatinine, age, is_male),
      ckm = (triglyceride >= cfg$ckm_tg) |
        (systolic_bp >= cfg$ckm_sbp | diastolic_bp >= cfg$ckm_dbp | htn_med) |
        mets |
        (fasting_blood_sugar >= cfg$ckm_fbs | dm_med) |
        (egfr < cfg$ckm_egfr) |
        (egfr >= cfg$ckm_egfr & uric_protein_qualitative %in% cfg$ckm_urine_protein_grades)
    )

  log(sprintf(
    "CKM stage 2/3 before washout: %s / %s members.",
    sum(cohort$ckm, na.rm = TRUE), nrow(cohort)
  ))

  cohort <- cohort %>% filter(ckm)

  # Washout (Exclusion Criterion 1): exclude anyone with a prior diagnosis
  # of CHD/Stroke/HF (any claim type, not admission-restricted) strictly
  # before their index date. A cardiovascular-death record can't logically
  # precede a living member's own index date, so in practice washout is
  # driven by CHD/Stroke/HF history only.
  windows <- cohort %>% select(PERSON_ID, index_date)
  push_window_table(con, windows, name = "#washout_windows")
  win_tbl <- tbl(con, "#washout_windows")

  prior_dx_members <- diagnosis_claims(con, cfg, cfg$washout_icd10) %>%
    inner_join(win_tbl, by = "PERSON_ID") %>%
    filter(dx_date < index_date) %>%
    distinct(PERSON_ID) %>%
    collect() %>%
    pull(PERSON_ID)

  n_before <- nrow(cohort)
  cohort <- cohort %>% filter(!PERSON_ID %in% prior_dx_members)
  log(sprintf(
    "Washout excluded %s members with prior CHD/Stroke/HF; %s remain.",
    n_before - nrow(cohort), nrow(cohort)
  ))

  # Exclusion Criterion 2: drop members missing any of the required
  # covariates (age, sex, eGFR, systolic BP, smoking/drinking/exercise
  # habit, HDL, LDL).
  n_before_missing <- nrow(cohort)
  complete <- complete.cases(cohort[, cfg$covariate_vars])
  cohort <- cohort[complete, ]
  log(sprintf(
    "Missing-covariate exclusion removed %s members; %s remain.",
    n_before_missing - nrow(cohort), nrow(cohort)
  ))

  cohort
}
