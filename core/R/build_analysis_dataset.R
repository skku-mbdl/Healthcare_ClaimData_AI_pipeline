# ==============================================================================
# Stage 2 shared: pivot drug exposures wide, merge onto cohort, train/test split.
#
# Unified from both countries' 03_build_analysis_dataset.R, which were
# already identical logic modulo column names (member_id/who_atc_code vs
# PERSON_ID/ATC_CODE) -- each country's 02_extract_drug_exposure.R now
# renames its output to the common patient_id/atc_code before calling this.
#
# Turns drug_exposure_long (one row per patient x eligible ATC code) into
# one binary column per ATC code (drug_<atc_code>, 1 = exposed in the
# look-back window, 0 = not), and merges it onto the cohort/outcome frame
# from core/R/extract_cohort_outcomes.R. A patient absent from
# drug_exposure_long entirely (no eligible-drug exposure at all) gets
# all-zero drug columns via the left_join + NA->0 fill below, not dropped.
#
# Patients missing ANY baseline covariate are dropped (defensively -- both
# countries' Stage 1 cohort-build already enforces covariate completeness
# for Korea, and Japan's does not currently enforce it the same way, see
# CLAUDE.md's "known asymmetries" note -- so this filter is load-bearing for
# Japan and a no-op safety check for Korea). Done once here, on the full
# analysis dataset, rather than left to each outcome's own complete.cases()
# check downstream, so every outcome's Boruta/Cox run sees the exact same
# population and reported Ns are directly comparable across outcomes.
#
# `cfg` must have TRAIN_FRACTION and SPLIT_SEED.
# ==============================================================================

library(dplyr)

build_analysis_dataset <- function(cohort, drug_exposure_long, drug_name_lookup,
                                    baseline_covariates, output_dir, cfg, log = message) {

  drug_codes <- sort(unique(drug_exposure_long$atc_code))
  log(sprintf("Pivoting %d-drug exposure matrix to wide binary format...", length(drug_codes)))

  # xtabs() is a vectorized C-level cross-tabulation -- avoids a tidyr
  # dependency and an R-level loop over patients for what is otherwise just
  # a patient x drug indicator matrix.
  wide_tab <- xtabs(~ patient_id + atc_code, data = drug_exposure_long)
  drug_wide <- as.data.frame.matrix(wide_tab)
  drug_wide[] <- lapply(drug_wide, function(x) as.integer(x > 0))
  colnames(drug_wide) <- make.names(paste0("drug_", colnames(drug_wide)), unique = TRUE)
  drug_wide$patient_id <- rownames(wide_tab)
  rownames(drug_wide) <- NULL

  rm(wide_tab)
  gc(full = TRUE)

  drug_cols <- setdiff(colnames(drug_wide), "patient_id")

  # drug_name_lookup keys were also derived from atc_code, in the same sort
  # order the SQL query returned them, so re-derive the column-name <->
  # ATC-code mapping the same way colnames(drug_wide) was built (make.names
  # on the sorted unique codes) rather than re-deriving it from column-name
  # text later, which would be fragile if a code ever collides after
  # sanitization.
  drug_col_lookup <- data.frame(
    atc_code = drug_codes,
    col_name = make.names(paste0("drug_", drug_codes), unique = TRUE),
    stringsAsFactors = FALSE
  ) %>%
    left_join(drug_name_lookup, by = "atc_code")

  dfd <- cohort %>%
    mutate(patient_id = as.character(patient_id)) %>%
    left_join(drug_wide %>% mutate(patient_id = as.character(patient_id)), by = "patient_id")

  dfd[drug_cols] <- lapply(dfd[drug_cols], function(x) {
    x[is.na(x)] <- 0L
    as.integer(x)
  })

  covariate_complete <- complete.cases(dfd[, baseline_covariates])
  n_dropped <- sum(!covariate_complete)
  if (n_dropped > 0) {
    missing_counts <- sapply(dfd[, baseline_covariates], function(x) sum(is.na(x)))
    log(sprintf(
      "Dropping %d/%d patients with a missing value in >=1 baseline covariate (per-covariate NA count: %s).",
      n_dropped, nrow(dfd),
      paste(sprintf("%s=%d", baseline_covariates, missing_counts), collapse = ", ")
    ))
    dfd <- dfd[covariate_complete, ]
  }

  log(sprintf(
    "Analysis dataset ready: %d patients x %d drug predictor columns (+ %d baseline covariates + outcomes).",
    nrow(dfd), length(drug_cols), length(baseline_covariates)
  ))

  # Train/test split: Boruta selection and Cox fitting run on dfd_train
  # only; validation runs on dfd_test only, so "which drugs matter" and "how
  # good is the model" are never answered on the same rows the model was
  # built from. Stratified on MACE status (the broadest composite outcome)
  # so the rarer ASCVD/HF/CV_DEATH events are still reasonably represented
  # in both splits.
  set.seed(cfg$SPLIT_SEED)
  event_idx    <- which(dfd$status_mace_10 == 1)
  nonevent_idx <- which(dfd$status_mace_10 != 1)
  train_ids <- c(
    sample(event_idx, size = round(cfg$TRAIN_FRACTION * length(event_idx))),
    sample(nonevent_idx, size = round(cfg$TRAIN_FRACTION * length(nonevent_idx)))
  )
  dfd$split <- "test"
  dfd$split[train_ids] <- "train"

  dfd_train <- dfd[dfd$split == "train", ]
  dfd_test  <- dfd[dfd$split == "test", ]

  # na.rm=TRUE: rows with NA status_mace_10 (the negative-survival-time edge
  # case in make_time_status()) are excluded from both event_idx/nonevent_idx
  # above, so sample() never assigns them to train -- they all land in
  # dfd_test by dfd$split's "test" default. Without na.rm here that NA
  # silently poisons the whole sum.
  log(sprintf(
    "Train/test split: %d train (%d MACE events) / %d test (%d MACE events).",
    nrow(dfd_train), sum(dfd_train$status_mace_10 == 1, na.rm = TRUE),
    nrow(dfd_test), sum(dfd_test$status_mace_10 == 1, na.rm = TRUE)
  ))

  log("Saving checkpoint (results/sideeffect_predict/analysis_dataset.rds)...")
  saveRDS(
    list(data = dfd, drug_cols = drug_cols, drug_col_lookup = drug_col_lookup, baseline_covariates = baseline_covariates),
    file = file.path(output_dir, "analysis_dataset.rds")
  )
  log("Checkpoint saved.")

  list(
    dfd = dfd, dfd_train = dfd_train, dfd_test = dfd_test,
    drug_cols = drug_cols, drug_col_lookup = drug_col_lookup
  )
}
