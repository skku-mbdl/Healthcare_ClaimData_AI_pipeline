# ==============================================================================
# countries/korea/run_pipeline.R -- Korea's FULL integrated pipeline.
#
#   Rscript countries/korea/run_pipeline.R
#
# Runs Stage 1 (drug screen: cohort-build -> Cox screen over the NHIS
# Sample Cohort) and, immediately after, Stage 2 (side-effect model:
# Boruta selection -> cause-specific Cox models -> validation -> SHAP) in
# one automatic run -- no manual copying of candidate_drugs.csv between
# projects, no separate invocation of a second script. This is the merge of
# what used to be two independent projects (CKM_Drug/Korea/NHIS_SAMPLE and
# CKM_PREVENT/KOREA/SideEffect_Model).
#
# Hand-off: Stage 1 writes its cohort/drug-exposure/Cox-screen output to
# BOTH the CKM_DRUG SQL Server database (nhis_cohort_processed/
# nhis_drug_exposure/nhis_cox_results -- same table names the original
# pipeline used, since this is a drop-in replacement against the same
# database Japan's Stage 1 also writes to) AND local file checkpoints
# under results/stage1_drug_screen/<run_id>/*.rds. Stage 2 reads from the
# database when reachable, else automatically falls back to the latest
# file checkpoint (see core/R/db_common.R's read_table_or_checkpoint()) --
# this whole pipeline can run and be smoke-tested with no live database at
# all, it just can't write shared processed tables for other consumers in
# that mode.
#
# Run scripts/00_setup.R-equivalent package installs are NOT handled here --
# see countries/korea/stage1_drug_screen/scripts/utils/config.R's package
# list and each stage2 00_config.R's library()/install.packages() calls;
# install those once per machine before the first run (an ODBC driver for
# SQL Server -- "ODBC Driver 17 for SQL Server" -- must also already be
# installed at the OS level).
# ==============================================================================

# ---- Bootstrap: locate this file, regardless of the caller's working
# directory (Rscript's own working directory is irrelevant to --file=). ----
.args <- commandArgs(trailingOnly = FALSE)
.file_flag <- sub("^--file=", "", .args[grep("^--file=", .args)])
if (length(.file_flag) != 1) {
  stop("Run this file with: Rscript countries/korea/run_pipeline.R (from any working directory).")
}
this_file <- normalizePath(.file_flag)
country_dir <- dirname(this_file)                 # .../countries/korea
pipeline_root <- dirname(dirname(country_dir))     # Healthcare_Data_AI_pipeline root

old_wd <- getwd()
setwd(country_dir)

core_dir <- file.path(pipeline_root, "core", "R")
source(file.path(core_dir, "logging_utils.R"))
source(file.path(core_dir, "db_common.R"))
source(file.path(core_dir, "clinical_calcs.R"))
source(file.path(core_dir, "cox_screen.R"))
source(file.path(core_dir, "candidate_drugs.R"))
source(file.path(core_dir, "extract_cohort_outcomes.R"))
source(file.path(core_dir, "build_analysis_dataset.R"))
source(file.path(core_dir, "validation_helpers.R"))
source(file.path(core_dir, "boruta_selection.R"))
source(file.path(core_dir, "cox_models.R"))
source(file.path(core_dir, "validate_models.R"))
source(file.path(core_dir, "shap_explain.R"))
source(file.path(core_dir, "export_boruta_results.R"))

reference_data_dir <- file.path(pipeline_root, "reference_data", "korea")

overall_result <- tryCatch({

  # ============================================================================
  # STAGE 1: drug screen (CKM_Drug/Korea/NHIS_SAMPLE equivalent)
  # ============================================================================
  stage1_dir <- file.path(country_dir, "stage1_drug_screen", "scripts")

  source(file.path(stage1_dir, "utils", "config.R"))        # -> cfg (Stage 1 study-definition constants)
  source(file.path(stage1_dir, "utils", "db_connect.R"))    # -> connect_nhis(), connect_ckm_drug()
  source(file.path(stage1_dir, "utils", "diagnosis_utils.R"))
  source(file.path(stage1_dir, "utils", "drug_utils.R"))
  source(file.path(stage1_dir, "01_load_cohort_base.R"))
  source(file.path(stage1_dir, "02_build_cohort.R"))
  source(file.path(stage1_dir, "03_covariates.R"))
  source(file.path(stage1_dir, "04_medication_flags.R"))
  source(file.path(stage1_dir, "05_outcomes.R"))
  source(file.path(stage1_dir, "06_drug_exposure.R"))

  run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
  stage1_results_path <- file.path("results", "stage1_drug_screen", run_id)
  dir.create(stage1_results_path, recursive = TRUE)

  logger <- new_logger(run_id, file.path("logs", "stage1_drug_screen"))
  stage1_log <- logger$log
  stage1_log(sprintf("Stage 1 (drug screen) run %s starting. Results -> %s", run_id, stage1_results_path))

  con_nhis <- connect_nhis()
  con_ckm <- tryCatch(connect_ckm_drug(), error = function(e) {
    stage1_log(sprintf("Could not open a CKM_DRUG connection (%s) -- Stage 1 will still write file checkpoints, but not to the database.", conditionMessage(e)))
    NULL
  })

  stage1_cox_results <- tryCatch({

    check_connection(con_nhis, "NHIS")
    if (!is.null(con_ckm)) check_connection(con_ckm, "CKM_DRUG")
    stage1_log("Database connection(s) verified.")

    atc_map <- load_atc_map(file.path(reference_data_dir, cfg$atc_map_file))
    push_atc_map(con_nhis, atc_map)
    stage1_log(sprintf("Loaded ATC map: %s GNL_NM_CD -> ATC_CODE pairs.", nrow(atc_map)))

    cohort_base <- load_cohort_base(con_nhis, cfg, stage1_log)
    cohort_base <- add_medication_flags(con_nhis, cohort_base, cfg, stage1_log)
    cohort <- build_cohort(con_nhis, cohort_base, cfg, stage1_log)
    cohort <- add_covariates(con_nhis, cohort, cfg, stage1_log)
    cohort <- add_outcomes(con_nhis, cohort, cfg, stage1_log)
    pdc <- add_drug_exposure(con_nhis, cohort, cfg, stage1_log)

    # Standardize id/ATC column names for the shared Cox screen (core/R/cox_screen.R).
    cohort_for_cox <- dplyr::rename(cohort, patient_id = PERSON_ID)
    pdc_for_cox <- dplyr::rename(pdc, patient_id = PERSON_ID, atc_code = ATC_CODE)
    cox_results <- run_cox_screen(cohort_for_cox, pdc_for_cox, cfg, stage1_log)

    stage1_log("Writing Stage 1 file checkpoints...")
    write_checkpoint(stage1_results_path, "cohort_processed", cohort)
    write_checkpoint(stage1_results_path, "drug_exposure", pdc)
    write_checkpoint(stage1_results_path, "cox_results", cox_results)

    if (!is.null(con_ckm)) {
      stage1_log("Writing processed tables to CKM_DRUG (nhis_ prefixed)...")
      write_ckm_table(con_ckm, "nhis_cohort_processed", cohort)
      write_ckm_table(con_ckm, "nhis_drug_exposure", pdc)
      write_ckm_table(con_ckm, "nhis_cox_results", cox_results)
    } else {
      stage1_log("CKM_DRUG not reachable -- skipped database write; file checkpoints are this run's only output.")
    }

    stage1_log(sprintf("Stage 1 run %s complete.", run_id))
    cox_results

  }, error = function(e) {
    stage1_log(sprintf("Stage 1 run %s FAILED: %s", run_id, conditionMessage(e)))
    stop(e)
  }, finally = {
    DBI::dbDisconnect(con_nhis)
    if (!is.null(con_ckm)) DBI::dbDisconnect(con_ckm)
  })

  # ============================================================================
  # STAGE 2: side-effect model (CKM_PREVENT/KOREA/SideEffect_Model equivalent)
  # ============================================================================
  stage2_dir <- file.path(country_dir, "stage2_sideeffect_model", "scripts")

  source(file.path(stage2_dir, "00_config.R"))   # -> redefines cfg (Stage 2 shape), progress(), output_dir, figures_dir, stage1_results_root, connect_db(), cox_outcome_map, opens log_con/sink

  con <- tryCatch(connect_db(), error = function(e) {
    progress(sprintf("Could not open a CKM_DRUG connection (%s) -- Stage 2 will use Stage 1's file checkpoint instead.", conditionMessage(e)))
    NULL
  })

  tryCatch({
    withCallingHandlers({

      progress("========== Extract cohort + outcomes ==========")
      source(file.path(stage2_dir, "01_extract_cohort_outcomes.R"))

      progress("========== Extract drug exposure ==========")
      source(file.path(stage2_dir, "02_extract_drug_exposure.R"))

      progress("========== Build analysis dataset ==========")
      bd <- build_analysis_dataset(cohort, drug_exposure_long, drug_name_lookup, baseline_covariates, output_dir, cfg, progress)

      outcome_names <- c("mace", "ascvd", "hf", "cvdeath")
      outcome_display_names <- c(mace = "MACE", ascvd = "ASCVD", hf = "Heart Failure", cvdeath = "CV Death")

      progress("========== Boruta variable selection ==========")
      bor <- run_boruta_selection(bd$dfd_train, bd$drug_col_lookup, baseline_covariates,
                                   candidate_drugs, cox_outcome_map, outcome_names, output_dir, cfg, progress)

      progress("========== Fit Cox models ==========")
      cm <- fit_cox_models(bd$dfd_train, bd$dfd_test, baseline_covariates, bor$confirmed_drugs,
                            outcome_names, output_dir, cfg, progress)

      progress("========== Validate models ==========")
      vm <- validate_models(cm$dfd_test, outcome_names, outcome_display_names, cm$cox_models,
                             bor$confirmed_drugs, bd$drug_col_lookup, output_dir, figures_dir, cfg, progress)

      progress("========== SHAP explanation (XAI) ==========")
      run_shap_explain(cm$dfd_test, outcome_names, outcome_display_names, cm$cox_models,
                        bd$drug_col_lookup, output_dir, figures_dir, cfg, progress)

      progress("========== Export Boruta results ==========")
      export_boruta_results(output_dir, outcome_names, progress)

    }, error = function(e) {
      calls <- sys.calls()
      trace <- vapply(calls, function(cl) paste(deparse(cl, nlines = 1), collapse = " "), character(1))
      progress("Traceback (innermost call last):")
      for (line in trace) progress(sprintf("  %s", line))
    })
  }, error = function(e) {
    progress(sprintf("Stage 2 FAILED: %s", conditionMessage(e)))
    stop(e)
  }, finally = {
    if (!is.null(con) && DBI::dbIsValid(con)) DBI::dbDisconnect(con)
    progress(paste("Side-effect-prediction run finished:", format(Sys.time())))
    sink(type = "message")
    sink()
    close(log_con)
  })

  invisible(TRUE)

}, finally = {
  setwd(old_wd)
})
