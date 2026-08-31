# ==============================================================================
# Stage 2 shared: Boruta variable selection over the drug matrix.
#
# For each of the four outcomes (MACE/ASCVD/HF/CV_DEATH), runs Boruta
# against a binary "event by t0=10y" target -- Boruta needs a classification
# target, not a Surv() object, so this is a screening step, not the final
# model. The DRUG columns offered to each outcome's Boruta run are
# restricted to that outcome's own candidate == TRUE rows from Stage 1's Cox
# screen (core/R/candidate_drugs.R), via `cox_outcome_map` (see each
# country's stage2 00_config.R) mapping this pipeline's 4 outcomes onto
# Stage 1's 4 outcomes (identical mapping for both countries, since both
# Stage-1 pipelines define exactly stroke/heart_failure/chd/cv_death):
#   hf      -> heart_failure
#   cvdeath -> cv_death
#   ascvd   -> chd, stroke      (CHD's I20-25 range already covers what a
#                                 finer split would call MI/Fatal_CHD)
#   mace    -> chd, stroke, heart_failure, cv_death
# Baseline covariates ride along on every outcome's run (unrestricted) so
# the random forest can contextualize drug importance against known
# confounders, but only Boruta's verdict on the DRUG columns is used:
# baseline covariates are force-included as adjustment variables in
# core/R/cox_models.R regardless of what Boruta decides about them here.
#
# Runs on dfd_train ONLY -- validation needs a held-out test set this
# selection step never saw, or the reported discrimination/calibration
# would be optimistically biased by the same variable-selection process.
#
# Unified from both countries' 05_boruta_selection.R: Japan's original used
# a LIKE-prefix match (its old allow-list mixed WHO ATC levels 3/4/5);
# Korea's used an exact match (its allow-list was already single-level).
# Now that BOTH countries' candidate source is their own Stage 1 Cox screen
# -- always a single, exact ATC-code granularity on both sides (see
# core/R/candidate_drugs.R's header) -- the prefix-match path is dead code
# for both, so this shared version uses exact match only.
#
# `cfg` must have BORUTA_SEED, BORUTA_MAX_RUNS, BORUTA_P_VALUE,
# BORUTA_NUM_THREADS.
# ==============================================================================

drug_cols_for_outcome <- function(outcome, candidate_drugs, cox_outcome_map, drug_col_lookup) {
  codes <- candidate_atc_for_outcome(candidate_drugs, cox_outcome_map[[outcome]])
  eligible_atc <- toupper(trimws(drug_col_lookup$atc_code)) %in% codes
  drug_col_lookup$col_name[eligible_atc]
}

# Boruta (like randomForest under it) can't handle NA predictors -- drop
# rows with any missing baseline covariate or missing outcome status for
# THIS outcome's run. Logged explicitly rather than silently shrinking n.
complete_rows_for_outcome <- function(data, outcome, predictors, log = message) {
  status_col <- paste0("status_", outcome, "_10")
  cc <- complete.cases(data[, c(predictors, status_col)])
  log(sprintf(
    "  [%s] %d/%d rows complete (%d dropped for missing predictor/outcome data).",
    outcome, sum(cc), nrow(data), sum(!cc)
  ))
  cc
}

#' Returns list(boruta_results = <list of Boruta fits>, confirmed_drugs =
#' <list of confirmed drug col_name character vectors>), one per outcome in
#' `outcome_names`. Also writes results/sideeffect_predict/boruta_<outcome>.rds.
run_boruta_selection <- function(dfd_train, drug_col_lookup, baseline_covariates,
                                  candidate_drugs, cox_outcome_map, outcome_names,
                                  output_dir, cfg, log = message) {

  boruta_results <- vector("list", length(outcome_names))
  names(boruta_results) <- outcome_names

  confirmed_drugs <- vector("list", length(outcome_names))
  names(confirmed_drugs) <- outcome_names

  for (i in seq_along(outcome_names)) {
    outcome <- outcome_names[i]
    outcome_drug_cols <- drug_cols_for_outcome(outcome, candidate_drugs, cox_outcome_map, drug_col_lookup)
    boruta_predictors <- c(baseline_covariates, outcome_drug_cols)

    log(sprintf(
      "Running Boruta for outcome '%s' (%d candidate predictors: %d baseline + %d drug, drug set = Stage 1 candidate==TRUE for %s)...",
      outcome, length(boruta_predictors), length(baseline_covariates), length(outcome_drug_cols),
      paste(cox_outcome_map[[outcome]], collapse = "/")
    ))

    cc <- complete_rows_for_outcome(dfd_train, outcome, boruta_predictors, log = log)
    x <- dfd_train[cc, boruta_predictors]
    y <- factor(dfd_train[cc, paste0("status_", outcome, "_10")] == 1, levels = c(FALSE, TRUE), labels = c("no_event", "event"))

    set.seed(cfg$BORUTA_SEED)
    boruta_args <- list(
      x = x, y = y,
      pValue = cfg$BORUTA_P_VALUE,
      maxRuns = cfg$BORUTA_MAX_RUNS,
      doTrace = 1
    )
    if (!is.null(cfg$BORUTA_NUM_THREADS)) boruta_args$num.threads <- cfg$BORUTA_NUM_THREADS

    fit <- do.call(Boruta::Boruta, boruta_args)

    fit <- Boruta::TentativeRoughFix(fit)
    boruta_results[[outcome]] <- fit

    log(sprintf("  [%s] importance source: %s (num.threads=%s)",
                outcome, fit$impSource, ifelse(is.null(cfg$BORUTA_NUM_THREADS), "NA", cfg$BORUTA_NUM_THREADS)))

    decision <- fit$finalDecision
    n_confirmed <- sum(decision == "Confirmed")
    n_rejected  <- sum(decision == "Rejected")
    log(sprintf("  [%s] Boruta done: %d confirmed, %d rejected (of %d candidates).",
                outcome, n_confirmed, n_rejected, length(decision)))

    confirmed_drugs[[outcome]] <- intersect(names(decision)[decision == "Confirmed"], outcome_drug_cols)
    log(sprintf("  [%s] %d of those confirmed variables are drugs (rest were baseline covariates).",
                outcome, length(confirmed_drugs[[outcome]])))

    saveRDS(fit, file.path(output_dir, sprintf("boruta_%s.rds", outcome)))
  }

  log("Boruta variable selection complete for all four outcomes.")
  list(boruta_results = boruta_results, confirmed_drugs = confirmed_drugs)
}
