# Shared helpers for querying 처방전교부상세내역(60t, prescription detail) and
# mapping GNL_NM_CD (일반명코드) to WHO ATC code, reused by medication flags
# (04) and drug exposure/PDC (06). Ported unchanged from CKM_Drug/Korea/
# NHIS_SAMPLE's scripts/utils/drug_utils.R. Uses raw_tbl() from
# core/R/db_common.R (must be sourced first).
#
# By instruction, this pipeline only reads T1 (의과/보건기관) claims tables
# -- T2 (치과/한방) and T3 (약국) are never queried. See cfg$tbl in
# utils/config.R.
#
# Same structural note as diagnosis_utils.R: PERSON_ID lives only on 20t, so
# 60t is joined to 20t via KEY_SEQ per year before UNION ALL-ing across
# years.
#
# RECU_FR_DT is parsed the same way as diagnosis_utils.R (TRY_CONVERT(date,
# ..., 112), assumed YYYYMMDD).

library(DBI)
library(dplyr)
library(dbplyr)
library(readxl)

#' Read NHIS_ATC_MAPPED.xlsx (GNL_NM_CD -> ATC_CODE) once, as a plain data
#' frame. Called from run_pipeline.R and passed down, rather than re-read
#' per step.
load_atc_map <- function(path) {
  df <- readxl::read_excel(path, col_names = TRUE)
  names(df) <- toupper(trimws(names(df)))
  stopifnot(all(c("GNL_NM_CD", "ATC_CODE") %in% names(df)))
  df %>%
    transmute(GNL_NM_CD = as.character(GNL_NM_CD), ATC_CODE = as.character(ATC_CODE)) %>%
    distinct()
}

#' Push the ATC map to `con` as a temp table so the GNL_NM_CD -> ATC join
#' happens server-side instead of pulling every prescription-detail row
#' into R first.
push_atc_map <- function(con, atc_map, name = "#atc_map") {
  dplyr::copy_to(con, atc_map, name = name, overwrite = TRUE, temporary = TRUE)
}

#' Prescription-detail claims (with PERSON_ID attached) across every year in
#' cfg$study_years, already joined 60t -> 20t (T1 only).
prescription_claims_all_years <- function(con, cfg) {
  parts <- lapply(cfg$study_years, function(y) {
    t20 <- raw_tbl(con, cfg$tbl$claim20, y)
    t60 <- raw_tbl(con, cfg$tbl$claim60, y)
    if (is.null(t20) || is.null(t60)) return(NULL)

    t20 <- t20 %>% select(PERSON_ID, KEY_SEQ)
    t60 %>%
      inner_join(t20, by = "KEY_SEQ") %>%
      mutate(rx_date = sql("TRY_CONVERT(date, RECU_FR_DT, 112)")) %>%
      select(PERSON_ID, KEY_SEQ, GNL_NM_CD, rx_date, MDCN_EXEC_FREQ)
  })
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) == 0) {
    stop(sprintf("No usable %s/%s tables found for years %s.",
                  cfg$tbl$claim20, cfg$tbl$claim60, paste(range(cfg$study_years), collapse = "-")), call. = FALSE)
  }
  Reduce(dplyr::union_all, parts)
}

#' Prescription-detail claims joined to the ATC map already pushed to `con`
#' via push_atc_map(). Returns PERSON_ID, ATC_CODE (+ atc3, the 3-character
#' WHO ATC prefix used for HTN_MED/DM_MED matching), rx_date, and
#' MDCN_EXEC_FREQ (days of medication for that claim line, used as the
#' days-of-supply length starting at rx_date). Rows whose GNL_NM_CD isn't in
#' NHIS_ATC_MAPPED.xlsx are dropped (inner join) -- per the protocol's Task
#' Rule, only drugs described in that mapping are considered at all.
prescription_claims <- function(con, cfg, atc_map_tbl_name = "#atc_map") {
  rx <- prescription_claims_all_years(con, cfg)
  atc_map <- tbl(con, atc_map_tbl_name)

  rx %>%
    inner_join(atc_map, by = "GNL_NM_CD") %>%
    mutate(atc3 = substr(ATC_CODE, 1, 3)) %>%
    filter(!is.na(rx_date), !is.na(MDCN_EXEC_FREQ), MDCN_EXEC_FREQ > 0) %>%
    select(PERSON_ID, ATC_CODE, atc3, rx_date, MDCN_EXEC_FREQ)
}
