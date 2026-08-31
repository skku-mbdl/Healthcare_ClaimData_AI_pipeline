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
# Unified from both countries' 06_fit_cox_models.R -- and revised again
# during this pipeline's own pilot test (see CLAUDE.md "Bug fix" section for
# the full history). Japan's original copy computed the zero-reference
# linear predictor by direct column indexing (newdata[, names(coef(fit))]),
# specifically to avoid stats::model.matrix()'s footgun: model.matrix()/
# model.frame() silently fall back to evaluating a term in the formula's
# ORIGINAL environment when a name isn't found in `data`, which is what
# crashed a real 08_shap_explain.R run in the original project (kernelshap's
# perturbed `newdata` didn't carry a column under the exact name expected,
# surfacing as an opaque "object '<name>' not found" instead of a clear
# error). That direct-indexing approach silently assumed every predictor is
# a plain numeric/binary column, never an R factor -- true for Japan's
# ORIGINAL Stage 2 covariates (sex_female/smoker/drinker/exerciser, all
# pre-binarized to 0/1) but NOT true once this pipeline unified cohort
# extraction (core/R/extract_cohort_outcomes.R) onto Korea's original
# convention of keeping sex/smoking_habit/drinking_habit/exercise_habit as
# multi-level factors. coxph() dummy-codes a factor predictor into
# coefficient names like "sexMale"/"smoking_habit2" that don't match any
# column in `newdata` at all -- direct indexing can't work for these,
# confirmed by this pipeline's own pilot run (Japan's Cox-fitting step
# failed immediately with "undefined columns selected" the first time any
# outcome had a Boruta-confirmed predictor set including a factor column).
#
# The fix below goes back to model.matrix() (which DOES dummy-code factors
# correctly) but closes the exact hole that made Japan's original version
# unsafe, rather than reintroducing it: (1) explicitly check that `newdata`
# contains every raw variable name the fit's formula needs, and stop() with
# the missing names if not -- turning the old silent wrong-environment
# fallback into a loud, immediate, actionable error; (2) explicitly
# re-apply the fitted model's own factor levels (fit$xlevels) to newdata's
# factor columns before building the design matrix, so a smaller/perturbed
# `newdata` (a SHAP background sample, a train/test split, a competing-risk
# subset) that happens not to contain every level of a multi-level
# covariate still produces the SAME dummy-coding columns the original fit
# used, rather than silently dropping a column and desyncing the design
# matrix from names(coef(fit)) -- the same class of bug, one level deeper.
#
# `cfg` must have T0_YEARS.
# ==============================================================================

compute_lp_zero_ref <- function(fit, newdata) {
  tt <- stats::delete.response(stats::terms(fit))

  needed_vars <- all.vars(tt)
  missing_vars <- setdiff(needed_vars, names(newdata))
  if (length(missing_vars) > 0) {
    stop(
      "compute_lp_zero_ref(): newdata is missing predictor column(s) the fit needs: ",
      paste(missing_vars, collapse = ", "),
      " -- refusing to silently fall back to model.matrix()'s enclosing-environment lookup."
    )
  }

  for (v in names(fit$xlevels)) {
    if (v %in% names(newdata)) newdata[[v]] <- factor(newdata[[v]], levels = fit$xlevels[[v]])
  }

  mm <- stats::model.matrix(tt, data = newdata)
  mm <- mm[, names(stats::coef(fit)), drop = FALSE]
  as.numeric(mm %*% stats::coef(fit))
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

# Fits both cause-specific models for ONE outcome and returns everything the
# caller needs. Split out of fit_cox_models() (which used to inline this
# body directly in a `for` loop, reusing one `model_data`/`form_event`
# variable across every outcome's iteration) because of a real bug this
# pipeline's own pilot run caught: `survival::basehaz()` (called inside
# predict_cif_10yr()) takes no `newdata` argument -- it always reconstructs
# the model's OWN TRAINING data by re-evaluating the fitted formula in the
# environment where coxph() was originally called. A `for` loop body is one
# shared execution frame across all iterations, so by the time
# core/R/shap_explain.R calls predict_cif_10yr() again later (well after
# fit_cox_models() has returned), that shared frame's `model_data` had
# already been overwritten by the LAST outcome's fit -- silently (if a
# column name happened to coincide) or loudly ("object '<drug column>' not
# found", which is what actually surfaced in the pilot run). Wrapping the
# per-outcome body in its OWN function call gives each outcome's `coxph()`
# fit a genuinely separate, permanent closure environment -- the standard R
# fix for this exact for-loop-closure footgun -- so basehaz() always finds
# the correct training data for whichever outcome's model it's asked about,
# no matter when or in what order that later call happens.
fit_one_outcome_cox <- function(outcome, dfd_train, dfd_test, baseline_covariates,
                                 confirmed_drugs, output_dir, cfg, log) {

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

  model <- list(
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
  pred_test <- rep(NA_real_, nrow(dfd_test))
  pred_test[test_cc] <- predict_cif_10yr(
    fit_event, fit_comp, dfd_test[test_cc, predictors, drop = FALSE], t0 = cfg$T0_YEARS
  )

  train_cc <- complete.cases(dfd_train[, predictors, drop = FALSE])
  pred_train <- rep(NA_real_, nrow(dfd_train))
  pred_train[train_cc] <- predict_cif_10yr(
    fit_event, fit_comp, dfd_train[train_cc, predictors, drop = FALSE], t0 = cfg$T0_YEARS
  )

  log(sprintf("  [%s] Cox models fit on %d train rows (%d events); predicted risk attached to train/test.",
              outcome, nrow(model_data), n_events))

  saveRDS(model, file.path(output_dir, sprintf("cox_model_%s.rds", outcome)))

  list(model = model, pred_train = pred_train, pred_test = pred_test)
}

#' Returns list(cox_models = <list of per-outcome fit lists>, dfd_train =,
#' dfd_test =) -- dfd_train/dfd_test get pred_<outcome> columns added.
#' Also writes results/sideeffect_predict/cox_model_<outcome>.rds.
fit_cox_models <- function(dfd_train, dfd_test, baseline_covariates, confirmed_drugs,
                            outcome_names, output_dir, cfg, log = message) {

  cox_models <- vector("list", length(outcome_names))
  names(cox_models) <- outcome_names

  for (outcome in outcome_names) {
    res <- fit_one_outcome_cox(outcome, dfd_train, dfd_test, baseline_covariates,
                                confirmed_drugs, output_dir, cfg, log)
    cox_models[[outcome]] <- res$model
    dfd_train[[paste0("pred_", outcome)]] <- res$pred_train
    dfd_test[[paste0("pred_", outcome)]] <- res$pred_test
  }

  log("Cox model fitting complete for all four outcomes.")
  list(cox_models = cox_models, dfd_train = dfd_train, dfd_test = dfd_test)
}
