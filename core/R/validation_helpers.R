# ==============================================================================
# Stage 2 shared: external-validation statistics for competing risks.
#
# aj_risk / harrell_c_with_ci / brier_ipcw / risk_group / calibration_table /
# royston / validate_one / plot_calibration_simple were already byte-for-byte
# identical between CKM_PREVENT/JAPAN/SideEffect_Model's and CKM_PREVENT/
# KOREA/SideEffect_Model's own 04_helpers_validation.R (both in turn ported
# verbatim from CKM_PREVENT/JAPAN/PREVENT/prevent_ckm_ver3.R) -- kept as one
# shared copy here instead of two. status coding (0/1/2 = censored/event/
# competing-event) and the t0=10 convention are relied on throughout; any new
# outcome added to either country's pipeline must follow the same coding.
#
# save_calibration_plot() takes `figures_dir` as an explicit argument (the
# two original per-country copies instead read it off a global set in each
# project's own 00_config.R) -- a plain-function version has no implicit
# global to read, since this file is now shared by both countries' pipelines.
# ==============================================================================

# Aalen-Johansen cumulative incidence at fixed time t0.
# status: 0 = censored, 1 = event of interest, 2 = competing event
aj_risk <- function(time, status, t0 = 10) {
  time   <- unname(as.numeric(time))
  status <- unname(as.integer(status))

  keep <- is.finite(time) & is.finite(status) & status %in% c(0L, 1L, 2L)
  time <- time[keep]
  status <- status[keep]

  n <- length(time)
  if (n == 0) return(NA_real_)

  ut  <- sort(unique(time))
  idx <- match(time, ut)
  n_ut <- length(ut)

  n_at_time  <- tabulate(idx, nbins = n_ut)
  d1_at_time <- tabulate(idx[status == 1L], nbins = n_ut)
  d_at_time  <- tabulate(idx[status %in% c(1L, 2L)], nbins = n_ut)

  y_at_time <- n - (cumsum(n_at_time) - n_at_time)

  sel <- ut <= t0 & d_at_time > 0
  if (!any(sel)) return(0)

  y  <- y_at_time[sel]
  d1 <- d1_at_time[sel]
  d  <- d_at_time[sel]

  S_before <- c(1, cumprod(1 - d / y))[seq_along(y)]

  sum(S_before * d1 / y)
}

# Harrell's C
harrell_c_with_ci <- function(time, status, risk) {
  time   <- unname(as.numeric(time))
  status <- unname(as.integer(status))
  risk   <- unname(as.numeric(risk))

  keep <- is.finite(time) & is.finite(status) & is.finite(risk)
  time <- time[keep]
  status <- status[keep]
  risk <- risk[keep]

  if (length(risk) == 0) {
    return(c(c = NA_real_, lower = NA_real_, upper = NA_real_))
  }

  fit <- survival::concordance(
    survival::Surv(time, status == 1) ~ risk,
    reverse = TRUE
  )

  cval <- unname(fit$concordance)
  se   <- sqrt(unname(fit$var))

  lower <- cval - 1.96 * se
  upper <- cval + 1.96 * se

  c(c = cval, lower = lower, upper = upper)
}

# IPCW Brier score at fixed horizon. Binary target: event of interest by t0.
# Competing event before t0 is treated as non-event for this Brier score.
brier_ipcw <- function(time, status, risk, t0 = 10) {
  time   <- unname(as.numeric(time))
  status <- unname(as.integer(status))
  risk   <- unname(as.numeric(risk))

  keep <- is.finite(time) & is.finite(status) & is.finite(risk) &
    status %in% c(0L, 1L, 2L)

  time <- time[keep]
  status <- status[keep]
  risk <- risk[keep]

  if (length(risk) == 0) return(NA_real_)

  kmc <- survival::survfit(survival::Surv(time, status == 0L) ~ 1)

  G <- function(tt) {
    g <- summary(kmc, times = tt, extend = TRUE)$surv
    pmax(g, 1e-6)
  }

  y <- as.numeric(time <= t0 & status == 1L)

  w <- rep(0, length(time))

  i_eventfree_t0 <- which(time >= t0)
  if (length(i_eventfree_t0) > 0) {
    w[i_eventfree_t0] <- 1 / G(t0)
  }

  i_known_before_t0 <- which(time < t0 & status != 0L)
  if (length(i_known_before_t0) > 0) {
    w[i_known_before_t0] <- 1 / G(time[i_known_before_t0])
  }

  mean(w * (y - risk)^2)
}

risk_group <- function(x, g = 10) {
  x <- as.numeric(x)
  q <- quantile(x, probs = seq(0, 1, length.out = g + 1), na.rm = TRUE, type = 2)
  q <- unique(q)
  if (length(q) < 3) return(rep(1L, length(x)))
  cut(x, breaks = q, include.lowest = TRUE, labels = FALSE)
}

calibration_table <- function(time, status, risk, t0 = 10, g = 10) {
  time   <- unname(as.numeric(time))
  status <- unname(as.integer(status))
  risk   <- unname(as.numeric(risk))

  keep <- is.finite(time) & is.finite(status) & is.finite(risk) &
    status %in% c(0L, 1L, 2L)

  time <- time[keep]
  status <- status[keep]
  risk <- risk[keep]

  grp <- risk_group(risk, g = g)

  idx_by_group <- split(seq_along(grp), grp)
  grp_keys <- as.integer(names(idx_by_group))

  out <- Map(function(k, i) {
    data.frame(
      group = k,
      n     = length(i),
      pred  = mean(risk[i]),
      obs   = aj_risk(time[i], status[i], t0 = t0)
    )
  }, grp_keys, idx_by_group)

  do.call(rbind, out)
}

# Royston-Sauerbrei D/R2 measures; degrades to NA if survival::royston() is
# not available in the installed survival build, rather than aborting.
royston <- function(time, status, risk) {
  d <- data.frame(
    time = as.numeric(time),
    event = as.integer(status == 1),
    lp = qlogis(as.numeric(risk))
  )
  d <- d[complete.cases(d) & is.finite(d$lp), ]
  if (nrow(d) == 0 || sum(d$event) == 0) return(NA)
  fit <- survival::coxph(survival::Surv(time, event) ~ lp, data = d)

  if (!exists("royston", where = asNamespace("survival"), inherits = FALSE)) {
    warning("survival::royston() is not available in this survival package version; ",
            "Royston D/R2 will be reported as NA.")
    return(NA)
  }
  survival::royston(fit)
}

validate_one <- function(data,
                          time_var,
                          status_var,
                          pred_var,
                          brier_time_var = NULL,
                          brier_status_var = NULL,
                          t0 = 10,
                          g = 10) {
  time   <- unname(as.numeric(data[[time_var]]))
  status <- unname(as.integer(data[[status_var]]))
  risk   <- unname(as.numeric(data[[pred_var]]))

  keep <- is.finite(time) &
    is.finite(status) &
    is.finite(risk) &
    status %in% c(0L, 1L, 2L)

  time   <- time[keep]
  status <- status[keep]
  risk   <- risk[keep]

  if (length(risk) == 0) {
    stop(paste0("No complete observations for ", pred_var))
  }
  obs <- aj_risk(time, status, t0 = t0)
  exp <- mean(risk)

  cal <- calibration_table(time = time, status = status, risk = risk, t0 = t0, g = g)

  slope <- if (!is.null(cal) && nrow(cal) >= 2) {
    unname(coef(lm(obs ~ pred, data = cal, weights = n))[2])
  } else {
    NA_real_
  }
  hc <- harrell_c_with_ci(time = time, status = status, risk = risk)
  rs <- royston(time = time, status = status, risk = risk)

  royston_D <- NA_real_
  royston_R2 <- NA_real_
  if (!is.null(rs) && length(rs) > 0) {
    if ("D" %in% names(rs)) royston_D <- unname(rs["D"])
    if ("R.D" %in% names(rs)) royston_R2 <- unname(rs["R.D"])
  }

  if (!is.null(brier_time_var) && !is.null(brier_status_var)) {
    brier_time_all   <- unname(as.numeric(data[[brier_time_var]]))
    brier_status_all <- unname(as.integer(data[[brier_status_var]]))
    brier_time   <- brier_time_all[keep]
    brier_status <- brier_status_all[keep]
    brier_value <- brier_ipcw(time = brier_time, status = brier_status, risk = risk, t0 = t0)
  } else {
    brier_value <- brier_ipcw(time = time, status = status, risk = risk, t0 = t0)
  }

  list(
    summary = data.frame(
      N = length(risk),
      observed_10y_risk        = obs,
      mean_predicted_10y_risk  = exp,
      OE_ratio                 = obs / exp,
      calibration_in_the_large = obs - exp,
      calibration_slope        = slope,
      harrell_c                = unname(hc["c"]),
      harrell_c_lower          = unname(hc["lower"]),
      harrell_c_upper          = unname(hc["upper"]),
      brier_10y                = brier_value,
      royston_D                = royston_D,
      royston_R2               = royston_R2
    ),
    calibration = cal
  )
}

# Base-R calibration plot (predicted vs. observed risk per decile, 45deg
# reference line).
plot_calibration_simple <- function(data, time_var, status_var, pred_var,
                                     t0 = 10, g = 10,
                                     title_txt = "Calibration plot",
                                     xlab = "Predicted 10-year risk",
                                     ylab = "Observed 10-year risk") {
  time   <- unname(as.numeric(data[[time_var]]))
  status <- unname(as.integer(data[[status_var]]))
  risk   <- unname(as.numeric(data[[pred_var]]))

  keep <- is.finite(time) & is.finite(status) & is.finite(risk) &
    status %in% c(0L, 1L, 2L)

  time   <- time[keep]
  status <- status[keep]
  risk   <- risk[keep]

  if (length(risk) == 0) stop("No complete observations available.")

  cal <- calibration_table(time, status, risk, t0 = t0, g = g)

  lim_max <- max(c(cal$pred, cal$obs), na.rm = TRUE)
  lim_max <- max(lim_max, 0.01)

  plot(
    cal$pred, cal$obs,
    type = "n",
    xlim = c(0, lim_max),
    ylim = c(0, lim_max),
    xlab = xlab,
    ylab = ylab,
    main = title_txt,
    las  = 1,
    bty  = "l"
  )

  abline(0, 1, lty = 2, lwd = 1.5)

  ord <- order(cal$pred)
  lines(cal$pred[ord], cal$obs[ord], lwd = 1.5)
  points(cal$pred, cal$obs, pch = 16, cex = 1.0)

  invisible(cal)
}

# png() needs an actual device open to produce a file -- without this, in a
# headless Rscript run there's no device backing the plot calls.
save_calibration_plot <- function(figures_dir, filename, plot_expr) {
  png(file.path(figures_dir, filename), width = 800, height = 800, res = 120)
  on.exit(dev.off(), add = TRUE)
  plot_expr
}
