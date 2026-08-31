# Step 2: MET1-5 / MetS / CKM stage 2-3 flag, then washout exclusion.
# Ported from CKM_Drug/Japan/scripts/02_build_cohort.R -- the only change is
# calling core/R/clinical_calcs.R's egfr(creatinine, age, is_male) with an
# explicit is_male boolean (sex == "Male") instead of the country's own copy
# of egfr() that took the raw sex string directly; see core/R/
# clinical_calcs.R's header for why.
#
# NOTE (known asymmetry, ported as-is from the original CKM_Drug/Japan
# pipeline): unlike Korea's 02_build_cohort.R, this step does NOT drop
# members with a missing baseline covariate (Korea's protocol has an
# explicit "Exclusion Criterion 2" for this; Japan's design doc/pipeline
# never added the equivalent check). core/R/build_analysis_dataset.R's
# Stage 2 covariate-completeness filter therefore does real work for Japan
# (not just a defensive no-op like it is for Korea).

library(dplyr)

build_cohort <- function(con, cohort_base, cfg, log = message) {

  cohort <- cohort_base %>%
    mutate(
      is_male = sex == "Male",
      met1 = (is_male & abdominal_circumference >= cfg$met_waist_male) |
        (!is_male & abdominal_circumference >= cfg$met_waist_female),
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

  # Washout: exclude anyone with a prior diagnosis of CHD/Stroke/HF (any
  # claim type) strictly before their index date. A cardiovascular-death
  # record can't logically precede a living member's own index date, so in
  # practice washout is driven by CHD/Stroke/HF history only.
  windows <- cohort %>% select(member_id, index_date)
  push_window_table(con, windows, name = "#washout_windows")
  win_tbl <- tbl(con, "#washout_windows")

  prior_dx_members <- diagnosis_claims(con, cfg$washout_icd10) %>%
    inner_join(win_tbl, by = "member_id") %>%
    filter(dx_date < index_date) %>%
    distinct(member_id) %>%
    collect() %>%
    pull(member_id)

  n_before <- nrow(cohort)
  cohort <- cohort %>% filter(!member_id %in% prior_dx_members)
  log(sprintf(
    "Washout excluded %s members with prior CHD/Stroke/HF; %s remain.",
    n_before - nrow(cohort), nrow(cohort)
  ))

  cohort
}
