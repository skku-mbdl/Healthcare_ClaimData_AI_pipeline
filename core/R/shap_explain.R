# ==============================================================================
# Stage 2 shared: SHAP explanation of each outcome's predicted 10-year risk.
#
# Explains predict_cif_10yr() (core/R/cox_models.R) directly -- the actual
# 10-year absolute-risk prediction this pipeline reports -- rather than the
# Cox models' linear predictor/hazard ratio. That prediction is a nonlinear
# combination of two cause-specific Cox models' baseline hazards (competing-
# risks Aalen-Johansen construction), so a closed-form "linear SHAP" on the
# log-hazard scale would NOT decompose the number this project actually
# cares about. This instead treats predict_cif_10yr() as a black box and
# uses model-agnostic Kernel SHAP (the `kernelshap` package).
#
# Runs on dfd_test (held-out, the same population validate_models() reports
# performance on) -- explaining predictions on the rows the models were fit
# on would describe how the model fits its training data, not how it
# behaves on new patients.
#
# Unified from both countries' 08_shap_explain.R (already functionally
# identical). `cfg` must have SHAP_MAX_ROWS, SHAP_BG_N, SHAP_SEED, T0_YEARS.
# ==============================================================================

run_shap_explain <- function(dfd_test, outcome_names, outcome_display_names,
                              cox_models, drug_col_lookup, output_dir, figures_dir,
                              cfg, log = message) {

  shap_results <- vector("list", length(outcome_names))
  names(shap_results) <- outcome_names

  for (outcome in outcome_names) {
    predictors <- cox_models[[outcome]]$predictors
    cc <- complete.cases(dfd_test[, predictors, drop = FALSE])
    X_full <- dfd_test[cc, predictors, drop = FALSE]

    if (nrow(X_full) == 0) {
      log(sprintf("  [%s] no complete-case test rows to explain -- skipping SHAP.", outcome))
      next
    }

    if (nrow(X_full) > cfg$SHAP_MAX_ROWS) {
      set.seed(cfg$SHAP_SEED)
      X <- X_full[sample(nrow(X_full), cfg$SHAP_MAX_ROWS), , drop = FALSE]
      log(sprintf(
        "  [%s] subsampling %d of %d complete-case test rows to explain (SHAP_MAX_ROWS cap).",
        outcome, cfg$SHAP_MAX_ROWS, nrow(X_full)
      ))
    } else {
      X <- X_full
    }

    # kernelshap::kernelshap()'s bg_n argument is ONLY honored when bg_X is
    # left NULL (it then auto-samples bg_n rows from X itself) -- if bg_X is
    # supplied explicitly, bg_n is silently ignored and every row of bg_X is
    # used as-is. So the cap has to be applied here, by hand, before bg_X is
    # ever passed in -- passing the full X_full as bg_X would silently
    # explain against the ENTIRE held-out test set as background instead of
    # SHAP_BG_N rows, which is both far slower than intended and can crash
    # deep inside kernelshap's internals at that scale.
    if (nrow(X_full) > cfg$SHAP_BG_N) {
      set.seed(cfg$SHAP_SEED)
      bg_X <- X_full[sample(nrow(X_full), cfg$SHAP_BG_N), , drop = FALSE]
    } else {
      bg_X <- X_full
    }

    # kernelshap's pred_fun signature is function(object, newdata) --
    # bundling both cause-specific fits + t0 (and the known-correct
    # predictor names) into `object`, since predict_cif_10yr() needs the
    # fits/t0 and kernelshap only threads one "object" through. newdata's
    # column names are defensively re-stamped from object$predictors (same
    # order kernelshap always uses) rather than trusted as-is.
    model_obj <- list(
      fit_event  = cox_models[[outcome]]$event,
      fit_comp   = cox_models[[outcome]]$competing,
      t0         = cfg$T0_YEARS,
      predictors = predictors
    )
    pred_fun <- function(object, newdata) {
      newdata <- as.data.frame(newdata)
      colnames(newdata) <- object$predictors
      predict_cif_10yr(object$fit_event, object$fit_comp, newdata, t0 = object$t0)
    }

    log(sprintf(
      "  [%s] running kernelshap over %d predictors x %d rows (background n=%d)...",
      outcome, length(predictors), nrow(X), nrow(bg_X)
    ))

    t_start <- proc.time()[["elapsed"]]
    ks <- kernelshap::kernelshap(
      object   = model_obj,
      X        = X,
      bg_X     = bg_X,
      pred_fun = pred_fun,
      verbose  = FALSE,
      seed     = cfg$SHAP_SEED
    )
    elapsed <- proc.time()[["elapsed"]] - t_start

    shap_vals <- ks$S

    importance <- data.frame(
      col_name      = colnames(shap_vals),
      mean_abs_shap = colMeans(abs(shap_vals)),
      mean_shap     = colMeans(shap_vals),
      stringsAsFactors = FALSE
    )
    # drug_col_lookup only maps drug_<code> columns to a readable name/ATC
    # code -- baseline covariates (age, sbp, etc.) aren't in it, so they
    # fall back to their own column name via the NA-coalesce below rather
    # than losing a label.
    importance <- merge(
      importance, drug_col_lookup[, c("col_name", "drug_label", "atc_code")],
      by = "col_name", all.x = TRUE
    )
    importance$drug_label <- ifelse(is.na(importance$drug_label), importance$col_name, importance$drug_label)
    importance <- importance[order(-importance$mean_abs_shap), ]

    shap_results[[outcome]] <- list(shap = shap_vals, X = X, baseline = ks$baseline, importance = importance)

    write.csv(importance, file.path(output_dir, sprintf("shap_importance_%s.csv", outcome)), row.names = FALSE)
    # Raw per-instance SHAP matrix + the feature values they were computed
    # on, kept for ad hoc follow-up (e.g. a per-patient waterfall) beyond
    # the global importance table above.
    saveRDS(list(shap = shap_vals, X = X, baseline = ks$baseline),
            file.path(output_dir, sprintf("shap_values_%s.rds", outcome)))

    top_n <- min(20, nrow(importance))
    top <- importance[seq_len(top_n), ]
    save_calibration_plot(
      figures_dir,
      sprintf("shap_importance_%s.png", outcome),
      {
        ord <- order(top$mean_abs_shap)
        par(mar = c(5, 14, 4, 2))
        barplot(
          top$mean_abs_shap[ord],
          names.arg = top$drug_label[ord],
          horiz = TRUE,
          las = 1,
          cex.names = 0.7,
          xlab = "mean |SHAP| (10-year predicted risk)",
          main = sprintf("SHAP importance: %s (drug-based model)", outcome_display_names[[outcome]])
        )
      }
    )

    log(sprintf(
      "  [%s] done in %.0fs. Top SHAP driver: %s (mean |SHAP| = %.4f).",
      outcome, elapsed, top$drug_label[which.max(top$mean_abs_shap)], max(top$mean_abs_shap)
    ))
  }

  log("SHAP explanation complete for all four outcomes.")
  shap_results
}
