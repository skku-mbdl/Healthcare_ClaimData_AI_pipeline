# ==============================================================================
# Stage 2 shared: held-out validation of the four drug-based Cox models.
#
# Runs validate_one() (core/R/validation_helpers.R) on dfd_test ONLY --
# pred_mace/pred_ascvd/pred_hf/pred_cvdeath were computed on the test set in
# core/R/cox_models.R using models fit exclusively on dfd_train, so these
# numbers reflect out-of-sample performance, not re-substitution.
#
# Unified from both countries' 07_validate_models.R (already functionally
# identical). `cfg` must have T0_YEARS.
# ==============================================================================

#' Writes results/sideeffect_predict/sideeffect_predict_result.csv,
#' confirmed_drugs_by_outcome.csv, calibration_<outcome>.csv, and
#' figures/sideeffect_predict/calibration_<outcome>.png. Returns list(
#' validation_results =, summary_table =, confirmed_drug_table =).
validate_models <- function(dfd_test, outcome_names, outcome_display_names,
                             cox_models, confirmed_drugs, drug_col_lookup,
                             output_dir, figures_dir, cfg, log = message) {

  validation_results <- vector("list", length(outcome_names))
  names(validation_results) <- outcome_names

  log("Validating drug-based Cox models on the held-out test set...")
  for (outcome in outcome_names) {
    res <- validate_one(
      data              = dfd_test,
      time_var          = paste0("time_", outcome, "_10"),
      status_var        = paste0("status_", outcome, "_10"),
      pred_var          = paste0("pred_", outcome),
      brier_time_var    = paste0("time_", outcome),
      brier_status_var  = paste0("status_", outcome),
      t0 = cfg$T0_YEARS, g = 10
    )
    validation_results[[outcome]] <- res
    log(sprintf("  [%s] Harrell's C = %.3f (%.3f-%.3f), O/E = %.3f, Brier = %.4f.",
                outcome, res$summary$harrell_c, res$summary$harrell_c_lower,
                res$summary$harrell_c_upper, res$summary$OE_ratio, res$summary$brier_10y))
  }

  summary_table <- do.call(rbind, lapply(outcome_names, function(o) {
    cbind(outcome = outcome_display_names[[o]], validation_results[[o]]$summary)
  }))
  print(summary_table)

  log(sprintf("Drawing calibration plots -> %s/...", figures_dir))
  for (outcome in outcome_names) {
    cal <- save_calibration_plot(
      figures_dir,
      sprintf("calibration_%s.png", outcome),
      plot_calibration_simple(
        data       = dfd_test,
        time_var   = paste0("time_", outcome, "_10"),
        status_var = paste0("status_", outcome, "_10"),
        pred_var   = paste0("pred_", outcome),
        title_txt  = sprintf("Calibration: %s (drug-based model)", outcome_display_names[[outcome]])
      )
    )
    write.csv(cal, file = file.path(output_dir, sprintf("calibration_%s.csv", outcome)), row.names = FALSE)
  }

  log(sprintf("Writing summary output -> %s/...", output_dir))
  write.csv(summary_table, file = file.path(output_dir, "sideeffect_predict_result.csv"), row.names = FALSE)

  # Confirmed-drug list per outcome, with the ATC code and the fitted Cox
  # coefficient/HR from the event-of-interest model -- this is the actual
  # "which drugs predict this side effect" deliverable, separate from the
  # aggregate discrimination/calibration numbers above.
  confirmed_drug_table <- do.call(rbind, lapply(outcome_names, function(o) {
    drugs <- confirmed_drugs[[o]]
    if (length(drugs) == 0) {
      # All FOUR columns must be zero-length, matching the populated branch
      # below column-for-column -- a pre-existing bug (present in the
      # original, pre-merge pipeline too, just never exercised) returned
      # data.frame(outcome = <scalar>, col_name = character(0)) here, which
      # errors ("arguments imply differing number of rows: 1, 0") the
      # moment ANY outcome has zero Boruta-confirmed drugs, since a
      # length-1 column can't be recycled against a length-0 one. Caught by
      # this pipeline's own demo pilot run, where Boruta legitimately
      # confirmed 0 drugs for 2 of Japan's 4 outcomes on synthetic data.
      return(data.frame(outcome = character(0), col_name = character(0), log_hr = numeric(0), hazard_ratio = numeric(0)))
    }
    fit <- cox_models[[o]]$event
    coefs <- coef(fit)
    hrs <- exp(coefs)
    drug_coefs <- coefs[intersect(drugs, names(coefs))]
    drug_hrs <- hrs[intersect(drugs, names(coefs))]
    data.frame(
      outcome  = outcome_display_names[[o]],
      col_name = names(drug_coefs),
      log_hr   = as.numeric(drug_coefs),
      hazard_ratio = as.numeric(drug_hrs)
    )
  }))

  confirmed_drug_table <- merge(
    confirmed_drug_table, drug_col_lookup[, c("col_name", "drug_label", "atc_code")],
    by = "col_name", all.x = TRUE
  )
  confirmed_drug_table <- confirmed_drug_table[order(confirmed_drug_table$outcome, -confirmed_drug_table$hazard_ratio), ]

  write.csv(confirmed_drug_table, file = file.path(output_dir, "confirmed_drugs_by_outcome.csv"), row.names = FALSE)

  log("Side-effect-prediction pipeline complete.")
  log(sprintf("Summary: %s", file.path(output_dir, "sideeffect_predict_result.csv")))
  log(sprintf("Confirmed drugs: %s", file.path(output_dir, "confirmed_drugs_by_outcome.csv")))

  list(validation_results = validation_results, summary_table = summary_table, confirmed_drug_table = confirmed_drug_table)
}
