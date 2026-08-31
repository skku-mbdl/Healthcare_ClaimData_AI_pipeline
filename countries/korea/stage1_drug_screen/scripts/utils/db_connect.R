# Korea-specific database connections. Generic helpers (check_connection,
# write_ckm_table, push_window_table, raw_tbl, union_years_tbl) live in
# core/R/db_common.R and are shared with Japan's stage1 adapter -- only the
# actual Server/Database values are country-specific, so this file is
# intentionally tiny.

library(DBI)
library(odbc)

connect_nhis <- function() {
  dbConnect(
    odbc(),
    Driver = "ODBC Driver 17 for SQL Server",
    Server = "SY_PC",
    Database = "NHIS",
    Trusted_Connection = "yes"
  )
}

connect_ckm_drug <- function() {
  dbConnect(
    odbc(),
    Driver = "ODBC Driver 17 for SQL Server",
    Server = "SY_PC",
    Database = "CKM_DRUG",
    Trusted_Connection = "yes"
  )
}
