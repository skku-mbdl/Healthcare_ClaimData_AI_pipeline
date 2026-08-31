# ==============================================================================
# Stage 2 shared: export each outcome's saved Boruta fit
# (results/sideeffect_predict/boruta_<outcome>.rds, written by
# core/R/boruta_selection.R) into one Excel workbook, one sheet per outcome.
#
# Unified from both countries' export_boruta_results.R (already functionally
# identical). Each country's stage2 run_pipeline.R calls
# export_boruta_results() as the pipeline's last automatic step -- unlike
# the original two projects, this is no longer a separate manual script run.
# ==============================================================================

# Reimplements Boruta::attStats() by hand -- ImpHistory/finalDecision are
# plain base-R objects (matrix/factor), so the saved .rds can be read
# without the Boruta package installed at all; only the *runtime*
# Boruta::Boruta() call itself needs the package, which this export step
# never invokes.
.attr_stats <- function(fit) {
  imp <- fit$ImpHistory
  real_attrs <- names(fit$finalDecision)
  imp <- imp[, real_attrs, drop = FALSE]

  stats <- apply(imp, 2, function(col) {
    finite <- col[is.finite(col)]
    c(
      meanImp   = mean(finite),
      medianImp = median(finite),
      minImp    = min(finite),
      maxImp    = max(finite),
      normHits  = mean(is.finite(col) & col > 0)
    )
  })

  data.frame(
    variable = real_attrs,
    type     = ifelse(startsWith(real_attrs, "drug_"), "drug", "baseline"),
    decision = as.character(fit$finalDecision[real_attrs]),
    t(stats),
    row.names = NULL,
    check.names = FALSE
  )
}

#' Writes results/sideeffect_predict/boruta_results.xlsx from the
#' boruta_<outcome>.rds files already saved by run_boruta_selection().
export_boruta_results <- function(output_dir, outcome_names, log = message) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required. Install it with install.packages('openxlsx').")
  }

  wb <- openxlsx::createWorkbook()
  decision_order <- c("Confirmed", "Tentative", "Rejected")

  for (outcome in outcome_names) {
    rds_path <- file.path(output_dir, sprintf("boruta_%s.rds", outcome))
    if (!file.exists(rds_path)) {
      warning(sprintf("Skipping '%s': %s not found.", outcome, rds_path))
      next
    }

    fit <- readRDS(rds_path)
    df <- .attr_stats(fit)
    df <- df[order(factor(df$decision, levels = decision_order), -df$meanImp), ]

    sheet_name <- toupper(outcome)
    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeData(wb, sheet_name, df, withFilter = TRUE)
    openxlsx::freezePane(wb, sheet_name, firstRow = TRUE)
    openxlsx::setColWidths(wb, sheet_name, cols = 1:ncol(df), widths = "auto")

    confirmed_rows <- which(df$decision == "Confirmed") + 1L  # +1 for header row
    if (length(confirmed_rows) > 0) {
      openxlsx::addStyle(
        wb, sheet_name,
        style = openxlsx::createStyle(fgFill = "#D9EAD3"),
        rows = confirmed_rows, cols = 1:ncol(df), gridExpand = TRUE, stack = TRUE
      )
    }
  }

  out_path <- file.path(output_dir, "boruta_results.xlsx")
  openxlsx::saveWorkbook(wb, out_path, overwrite = TRUE)
  log(sprintf("Saved Boruta results for %d outcomes to %s", length(outcome_names), out_path))
  invisible(out_path)
}
