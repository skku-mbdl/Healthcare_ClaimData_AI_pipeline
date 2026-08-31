# One-time environment setup for the WHOLE integrated pipeline (both
# countries, both stages). Run once per machine:
#   Rscript 00_setup.R
#
# Requires an ODBC driver for SQL Server ("ODBC Driver 17 for SQL Server")
# already installed at the OS level -- the odbc R package talks to it but
# does not install it. An earlier version of one of the projects this
# pipeline merges hit a Windows application-control policy (WDAC/AppLocker/
# Smart App Control) blocking odbc's compiled DLL from loading under a
# per-user/per-project library -- if that recurs, get R's library path
# allow-listed on the machine before re-running this script.
#
# Unlike the original CKM_Drug/CKM_PREVENT projects (each with their own
# renv project + package list), this is ONE shared package list for the
# whole integrated pipeline -- core/R/ and both countries' adapters all run
# in the same R session per country (see countries/<country>/run_pipeline.R).

options(repos = c(CRAN = "https://cloud.r-project.org"))

pkgs <- c(
  "DBI",
  "odbc",
  "dplyr",
  "dbplyr",
  "lubridate",
  "survival",
  "tibble",
  "readxl",     # Korea Stage 1: NHIS_ATC_MAPPED.xlsx
  "Boruta",     # Stage 2: variable selection
  "kernelshap", # Stage 2: SHAP explanation
  "openxlsx"    # Stage 2: Boruta results export
)

to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)

missing_after <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing_after) > 0) {
  stop("Failed to install: ", paste(missing_after, collapse = ", "))
}
cat("All required packages are installed.\n")
