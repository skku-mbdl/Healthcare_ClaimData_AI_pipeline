# Stage 1 shared: single-drug-at-a-time Cox proportional hazards screen,
# one ATC code at a time, against each of the four primary outcomes
# (stroke, heart_failure, chd, cv_death). A drug is flagged as a candidate
# cardiovascular side-effect drug when HR > cfg$candidate_hr_min with
# p < cfg$candidate_p_max.
#
# Unified from CKM_Drug/Japan/scripts/07_cox_analysis.R and CKM_Drug/Korea/
# NHIS_SAMPLE/scripts/07_cox_analysis.R, which were already identical except
# for column names (member_id/who_atc_code vs PERSON_ID/ATC_CODE). Both
# country adapters rename their cohort/pdc tables' id and ATC columns to the
# common `patient_id`/`atc_code` right before calling run_cox_screen() (see
# each country's stage1_drug_screen/scripts/run_pipeline.R) -- this file
# itself never sees a country-specific column name.
#
# Adjustment covariates are exactly cfg$covariate_vars (age, sex, eGFR,
# systolic BP, smoking/drinking/exercise habit, HDL, LDL -- identical list in
# both countries' config.R). Comorbidity flags are kept on the cohort table
# for descriptive/Table 1 purposes but are not forced into this adjustment
# set (same as both original pipelines).

library(dplyr)
library(survival)

prep_model_data <- function(cohort) {
  cohort %>%
    mutate(
      sex = factor(sex),
      smoking_habit = factor(smoking_habit),
      drinking_habit = factor(drinking_habit),
      exercise_habit = factor(exercise_habit)
    )
}

#' Build the Surv() ~ exposed + covariates formula for one outcome. Split
#' out so run_cox_screen() can build each outcome's formula once and reuse
#' it across every drug, instead of re-deriving the same formula string on
#' every one of the (drugs x outcomes) fits.
build_outcome_formula <- function(outcome, cfg) {
  as.formula(sprintf(
    "Surv(%s_time, %s_event) ~ exposed + %s",
    outcome, outcome, paste(cfg$covariate_vars, collapse = " + ")
  ))
}

fit_one_drug_outcome <- function(model_data, outcome, form) {
  fit <- tryCatch(
    coxph(form, data = model_data),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  s <- summary(fit)
  if (!"exposedTRUE" %in% rownames(s$coefficients)) return(NULL)

  row <- s$coefficients["exposedTRUE", , drop = FALSE]
  ci <- s$conf.int["exposedTRUE", , drop = FALSE]

  tibble::tibble(
    outcome = outcome,
    hr = row[, "exp(coef)"],
    ci_lower = ci[, "lower .95"],
    ci_upper = ci[, "upper .95"],
    p_value = row[, "Pr(>|z|)"],
    n = s$n,
    n_event = s$nevent
  )
}

#' cohort must have: patient_id, sex, smoking_habit, drinking_habit,
#' exercise_habit, cfg$covariate_vars, and <outcome>_time/<outcome>_event
#' for outcome in c("stroke","heart_failure","chd","cv_death").
#' pdc must have: patient_id, atc_code, exposed (logical).
run_cox_screen <- function(cohort, pdc, cfg, log = message) {

  model_data <- prep_model_data(cohort)
  outcomes <- c("stroke", "heart_failure", "chd", "cv_death")

  # Built once, reused for every drug -- identical across drugs for a given
  # outcome, so there's no reason to re-derive it (drugs x outcomes) times.
  formulas <- setNames(lapply(outcomes, build_outcome_formula, cfg = cfg), outcomes)

  atc_codes <- pdc %>%
    filter(exposed) %>%
    count(atc_code) %>%
    filter(n >= cfg$min_exposed_n) %>%
    pull(atc_code)

  log(sprintf(
    "Screening %s ATC codes with >= %s exposed members (of %s distinct codes total).",
    length(atc_codes), cfg$min_exposed_n, n_distinct(pdc$atc_code)
  ))

  # Split once into a named list of exposed patient_ids per ATC code, so each
  # loop iteration below is a plain list lookup instead of re-filtering the
  # full pdc table (which was O(rows in pdc) per iteration, x length(atc_codes)).
  exposed_ids_by_atc <- pdc %>%
    filter(exposed) %>%
    { split(.$patient_id, .$atc_code) }

  patient_ids <- model_data$patient_id

  results <- vector("list", length(atc_codes))
  for (i in seq_along(atc_codes)) {
    atc <- atc_codes[i]

    # Base-R column assignment instead of dplyr::mutate(): skips
    # tidyeval/tibble-rebuild overhead on a line that runs once per
    # screened drug. R's copy-on-write means the unchanged covariate
    # columns aren't duplicated here -- only the new `exposed` vector is.
    drug_data <- model_data
    drug_data$exposed <- patient_ids %in% exposed_ids_by_atc[[atc]]

    per_outcome <- lapply(outcomes, function(o) fit_one_drug_outcome(drug_data, o, formulas[[o]]))
    per_outcome <- bind_rows(per_outcome)
    if (nrow(per_outcome) > 0) per_outcome$atc_code <- atc
    results[[i]] <- per_outcome

    if (i %% 50 == 0) log(sprintf("  ...%s / %s ATC codes screened", i, length(atc_codes)))
  }

  all_results <- bind_rows(results)
  if (nrow(all_results) == 0) {
    log("No ATC code met the minimum exposed-member threshold; no models fit.")
    return(all_results)
  }

  all_results <- all_results %>%
    mutate(candidate = hr > cfg$candidate_hr_min & p_value < cfg$candidate_p_max) %>%
    select(atc_code, outcome, hr, ci_lower, ci_upper, p_value, n, n_event, candidate) %>%
    arrange(outcome, p_value)

  log(sprintf(
    "Cox screening complete: %s drug-outcome pairs fit, %s flagged as candidates.",
    nrow(all_results), sum(all_results$candidate)
  ))

  all_results
}
