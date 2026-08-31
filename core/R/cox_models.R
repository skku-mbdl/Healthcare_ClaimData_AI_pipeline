# ==============================================================================
# Stage 2 shared: cause-specific Cox models + 10-year competing-risks CIF.
#
# For each outcome, fits TWO cause-specific Cox models on dfd_train (Boruta-
# confirmed drugs + baseline covariates as predictors, competing events
# coded as censoring within each model in the usual cause-specific way):
#   - fit_event: Surv(time_10, status_10 == 1) ~ predictors   (event of interest)
#   - fit_comp:  Surv(time_10, status_10 == 2) ~ predictors   (competing: non-CV death)
#
# A single cause-specific hazard ratio does NOT translate directly into an
# absolute risk under competing risks -- predict_cif_10yr() combines both
# models' baseline cumulative hazards into a proper 10-year cumulative-
# incidence prediction per subject, the same Aalen-Johansen-style
# construction aj_risk() (core/R/validation_helpers.R) uses for the OBSERVED
# risk, so predicted vs. observed are on the same scale downstream.
#
# Unified from both countries' 06_fit_cox_models.R. Japan's copy had already
# been fixed (2026-08-31) to compute the zero-reference linear predictor by
# direct column indexing rather than stats::model.matrix() on the fit's
# terms -- model.matrix()/model.frame() fall back to evaluating a term in
# the formula's ORIGINAL environment when a name isn't found in `data`, and
# that fallback is what actually crashed a real 08_shap_explain.R run
# (kernelshap's perturbed `newdata` didn't carry a column under the exact
# name expected, surfacing as an opaque "object '<name>' not found"). Korea's
# copy still had the older model.matrix()-based version. This shared file
# uses Japan's fixed version for both countries -- Korea's Stage 2 therefore
# also gets this latent-crash fix, not just a plumbing merge.
#
# `cfg` must have T0_YEARS.
# ==============================================================================

# Linear predictor at the "all covariates = 0" reference -- matches
# basehaz(..., centered = FALSE)'s reference point. Every predictor here is
# a plain numeric/binary main effect (baseline covariates are numeric, not R
# factors; drug columns are 0/1) -- no interactions, splines, or categorical
# terms anywhere in these formulas -- so the design matrix is just newdata's
# named columns, in coefficient order, with no dummy-encoding needed. Direct
# column indexing (rather than model.matrix()) can't hit the formula-
# environment fallback described above: a genuinely missing column now
# fails loudly with "undefined columns selected" instead of an opaque
# "object not found". Verified equivalent to the old model.matrix()-based
# computation (0 numeric difference) before the switch.
compute_lp_zero_ref <- function(fit, newdata) {
  cf <- stats::coef(fit)
  mm <- as.matrix(newdata[, names(cf), drop = FALSE])
  storage.mode(mm) <- "double"
  as.numeric(mm %*% cf)
}

# 10-year cumulative incidence for the event of interest, combining
# fit_event's and fit_comp's baseline hazards:
#   S_i(u-)      = exp(-(H01(u-)*exp(lp1_i) + H02(u-)*exp(lp2_i)))
#   dCIF1_i(u)   = S_i(u-) * dH01(u) * exp(lp1_i)
#   CIF1_i(t0)   = sum over cause-1 jump times u <= t0 of dCIF1_i(u)
# Vectorized across subjects (outer product over subjects x cause-1 jump
# times) rather than a per-subject loop.
predict_cif_10yr <- function(fit_event, fit_comp, newdata, t0) {
  bh1 <- survival::basehaz(fit_event, centered = FALSE)
  bh2 <- survival::basehaz(fit_comp, centered = FALSE)
  bh1 <- bh1[order(bh1$time), ]
  bh2 <- bh2[order(bh2$time), ]

  u <- bh1$time
  dH01 <- diff(c(0, bh1$hazard))
  cumH01_before <- cumsum(dH01) - dH01

  # Step-function lookup: cause-2's baseline cumulative hazard just before
  # each of cause-1's jump times (right-continuous step function, so the
  # value strictly before u[k] is the last recorded H02 at a time < u[k]).
  idx2 <- findInterval(u - 1e-8, bh2$time)
  H02_before <- ifelse(idx2 >= 1, bh2$hazard[idx2], 0)

  keep <- u <= t0
  if (!any(keep)) return(rep(0, nrow(newdata)))
  cumH01_before <- cumH01_before[keep]
  H02_before <- H02_before[keep]
  dH01 <- dH01[keep]

  lp1 <- compute_lp_zero_ref(fit_event, newdata)
  lp2 <- compute_lp_zero_ref(fit_comp, newdata)
  exp_lp1 <- exp(lp1)
  exp_lp2 <- exp(lp2)

  H1mat <- outer(exp_lp1, cumH01_before)
  H2mat <- outer(exp_lp2, H02_before)
  Smat  <- exp(-(H1mat + H2mat))
  dCIF  <- Smat * outer(exp_lp1, dH01)

  as.numeric(rowSums(dCIF))
}

#' Returns list(cox_models = <list of per-outcome fit lists>, dfd_train =,
#' dfd_test =) -- dfd_train/dfd_test get pred_<outcome> columns added.
#' Also writes results/sideeffect_predict/cox_model_<outcome>.rds.
fit_cox_models <- function(dfd_train, dfd_test, baseline_covariates, confirmed_drugs,
                            outcome_names, output_dir, cfg, log = message) {

  cox_models <- vector("list", length(outcome_names))
  names(cox_models) <- outcome_names

  for (outcome in outcome_names) {
    predictors <- c(baseline_covariates, confirmed_drugs[[outcome]])
    log(sprintf(
      "Fitting cause-specific Cox models for '%s' (%d predictors: %d baseline + %d Boruta-confirmed drug).",
      outcome, length(predictors), length(baseline_covariates), length(confirmed_drugs[[outcome]])
    ))

    time_col <- paste0("time_", outcome, "_10")
    status_col <- paste0("status_", outcome, "_10")
    cc <- complete.cases(dfd_train[, c(predictors, time_col, status_col)])

    model_data <- dfd_train[cc, predictors, drop = FALSE]
    model_data$time_10 <- dfd_train[[time_col]][cc]
    model_data$status_10 <- dfd_train[[status_col]][cc]

    n_events <- sum(model_data$status_10 == 1)
    events_per_predictor <- n_events / length(predictors)
    if (events_per_predictor < 10) {
      warning(sprintf(
        "'%s': only %.1f events per predictor (%d events / %d predictors) -- below the usual >=10 rule of thumb; coefficients may be unstable.",
        outcome, events_per_predictor, n_events, length(predictors)
      ))
    }

    form_event <- as.formula(paste0("Surv(time_10, status_10 == 1) ~ ", paste(predictors, collapse = " + ")))
    form_comp  <- as.formula(paste0("Surv(time_10, status_10 == 2) ~ ", paste(predictors, collapse = " + ")))

    fit_event <- survival::coxph(form_event, data = model_data)
    fit_comp  <- survival::coxph(form_comp, data = model_data)

    cox_models[[outcome]] <- list(
      event = fit_event,
      competing = fit_comp,
      predictors = predictors,
      n_train = nrow(model_data),
      n_events_train = n_events
    )

    # Predicted 10-year absolute risk for BOTH splits: test set is what
    # validation uses for honest performance metrics; train set is kept too
    # for a quick apparent-vs-test-performance comparison if wanted.
    test_cc <- complete.cases(dfd_test[, predictors, drop = FALSE])
    dfd_test[[paste0("pred_", outcome)]] <- NA_real_
    dfd_test[[paste0("pred_", outcome)]][test_cc] <- predict_cif_10yr(
      fit_event, fit_comp, dfd_test[test_cc, predictors, drop = FALSE], t0 = cfg$T0_YEARS
    )

    train_cc <- complete.cases(dfd_train[, predictors, drop = FALSE])
    dfd_train[[paste0("pred_", outcome)]] <- NA_real_
    dfd_train[[paste0("pred_", outcome)]][train_cc] <- predict_cif_10yr(
      fit_event, fit_comp, dfd_train[train_cc, predictors, drop = FALSE], t0 = cfg$T0_YEARS
    )

    log(sprintf("  [%s] Cox models fit on %d train rows (%d events); predicted risk attached to train/test.",
                outcome, nrow(model_data), n_events))

    saveRDS(cox_models[[outcome]], file.path(output_dir, sprintf("cox_model_%s.rds", outcome)))
  }

  log("Cox model fitting complete for all four outcomes.")
  list(cox_models = cox_models, dfd_train = dfd_train, dfd_test = dfd_test)
}
