# Japan-specific database connections. Generic helpers (check_connection,
# write_ckm_table, push_window_table) live in core/R/db_common.R and are
# shared with Korea's stage1 adapter -- only the actual Server/Database
# values are country-specific.

library(DBI)
library(odbc)

connect_jmdc <- function() {
  dbConnect(
    odbc(),
    Driver = "ODBC Driver 17 for SQL Server",
    Server = "SY_PC",
    Database = "JMDC",
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
