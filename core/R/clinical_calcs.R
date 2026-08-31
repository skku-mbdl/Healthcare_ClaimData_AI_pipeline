# Sex-specific eGFR formula (2021 CKD-EPI creatinine equation, race-free),
# computed from serum creatinine and age. Used both in the CKM stage 2/3
# flag itself (each country's 02_build_cohort.R) and as a Cox covariate.
#
# Unified signature: takes an `is_male` logical rather than a raw sex code,
# so this one function serves both countries even though their raw sex
# coding differs (Japan: "Male"/"Female" strings; Korea: "1"/"2", the raw
# NHID_JK.SEX coding) -- each country's cohort-build script computes
# is_male itself, from whichever raw representation it has, right before
# calling this.
egfr <- function(creatinine, age, is_male) {
  ifelse(
    is.na(creatinine), NA_real_,
    ifelse(
      is_male,
      ifelse(
        creatinine <= 0.9,
        142 * (creatinine / 0.9)^-0.302 * 0.9938^age,
        142 * (creatinine / 0.9)^-1.2 * 0.9938^age
      ),
      ifelse(
        creatinine <= 0.7,
        142 * (creatinine / 0.7)^-0.241 * 0.9938^age * 1.012,
        142 * (creatinine / 0.7)^-1.2 * 0.9938^age * 1.012
      )
    )
  )
}
