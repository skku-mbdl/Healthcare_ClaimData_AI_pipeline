# Central configuration for the CKM drug-candidate pipeline (Korea / NHIS
# Sample Cohort). All study-definition constants live here so downstream
# scripts never hardcode a code list or raw table name twice. Ported
# unchanged from CKM_Drug/Korea/NHIS_SAMPLE/scripts/utils/config.R -- see
# that project's CLAUDE.md for the full methodology these encode and the
# schema-mapping decisions summarized next to each constant below.
#
# results_dir/logs_dir/atc_map_file are NOT set here (unlike the original) --
# path resolution (results/, logs/, reference_data/korea/NHIS_ATC_MAPPED.xlsx)
# is now owned by run_pipeline.R, which locates itself and the shared
# reference_data/ directory explicitly, rather than assuming a fixed
# relative layout from this config file.

cfg <- list(

  # ---- Study window ---------------------------------------------------
  # Index date = first health-checkup YEAR on record in [checkup_year_min,
  # checkup_year_max] -- the checkup table has no day/month, so the index
  # date itself is approximated as index_date_month_day of that year (see
  # 01_load_cohort_base.R). study_years is the full span of raw claims/
  # checkup/eligibility tables this pipeline ever needs to read: one year
  # before checkup_year_min (lookback can reach into the prior year) through
  # observation_end_cutoff's year.
  checkup_year_min = 2009L,
  checkup_year_max = 2012L,
  study_years = 2008:2013,
  index_date_month_day = "07-01",  # mid-year proxy
  lookback_days = 365,
  observation_end_cutoff = as.Date("2013-12-31"),

  # ---- Raw NHIS Sample Cohort table-name prefixes ----------------------
  # One flat file per (table type) x year x institution type was distributed
  # (e.g. NHID_JK_2009, NHID_GY40_T1_2010); this pipeline assumes the live
  # `NHIS` SQL Server database kept that exact per-year, per-type naming
  # (prefix + "_" + year). If the real database consolidated these
  # differently, this is the only place that needs to change (plus
  # raw_tbl()/union_years_tbl() in core/R/db_common.R).
  # 진료DB tables come in 3 institution-type variants (T1 = 의과/보건기관
  # medical + public health, T2 = 치과/한방 dental + oriental medicine,
  # T3 = 약국 pharmacy). By instruction, this pipeline only ever reads T1 --
  # T2 and T3 tables are never queried, for any of GY20/GY30/GY40/GY60.
  tbl = list(
    eligibility = "NHID_JK",       # 자격DB: PERSON_ID, SEX, AGE_GROUP, DTH_YM, DTH_CODE1/2, one row per person per year
    checkup     = "NHID_GJ",       # 건강검진DB: PERSON_ID, HCHK_YEAR, checkup values
    institution = "NHID_YK",       # 요양기관DB (not currently used by this pipeline)
    claim20 = "NHID_GY20_T1",      # 명세서 (claim header) -- T1 only
    claim30 = "NHID_GY30_T1",      # 진료내역 (treatment detail) -- T1 only, not currently read
    claim40 = "NHID_GY40_T1",      # 상병내역 (diagnosis detail) -- T1 only
    claim60 = "NHID_GY60_T1"       # 처방전교부상세내역 (prescription detail) -- T1 only
  ),

  # ---- FORM_CD (서식코드) classification --------------------------------
  admission_form_cd  = c("02", "04", "06", "07", "10", "12"),
  outpatient_form_cd = c("03", "05", "08", "11", "13"),

  # ---- Outcome ICD-10 codes (admission claims only) -------------------
  outcome_icd10 = list(
    stroke = c("I61", "I62", "I63"),
    heart_failure = c("I50"),
    chd = c("I20", "I21", "I22", "I23", "I24", "I25")
    # cardiovascular death is not code-based via SICK_SYM; see
    # 05_outcomes.R -- it is NHID_JK.DTH_YM non-missing with DTH_CODE1 or
    # DTH_CODE2 starting with "I", per the protocol.
  ),

  washout_icd10 = c("I61", "I62", "I63", "I50",
                     "I20", "I21", "I22", "I23", "I24", "I25"),

  # ---- Comorbidity ICD-10 codes ----------------------------------------
  comorbidity_icd10 = list(
    hypertension = c("I10", "I11", "I12", "I13", "I15"),
    diabetes = c("E10", "E11", "E12", "E13", "E14"),
    dyslipidemia = c("E78"),
    ckd = c("N18", "N19")
  ),

  # ---- Medication flags (HTN_MED / DM_MED), claims-based ---------------
  htn_med_atc_prefix = c("C02", "C03", "C07", "C08", "C09"),
  dm_med_atc_prefix  = c("A10"),

  # ---- CKM stage 2/3 cohort thresholds ----------------------------------
  met_waist_male = 90,
  met_waist_female = 80,
  met_hdl_male = 40,
  met_hdl_female = 50,
  met_tg = 150,
  met_sbp = 130,
  met_dbp = 80,
  met_fbs = 100,

  ckm_tg = 135,
  ckm_sbp = 140,
  ckm_dbp = 90,
  ckm_fbs = 126,
  ckm_egfr = 60,
  # OLIG_PROTE_CD has codes 1-6 (1=-, 2=+/-, 3=+1, 4=+2, 5=+3, 6=+4);
  # protocol's "grade in (3,4,5,6)" maps directly onto codes 3-6.
  ckm_urine_protein_grades = c(3, 4, 5, 6),

  # ---- Drug exposure / PDC ----------------------------------------------
  pdc_threshold = 0.8,
  min_exposed_n = 10,

  # ---- Candidate drug screening rule ------------------------------------
  candidate_hr_min = 1.00,
  candidate_p_max = 0.05,

  # ---- Covariates used in the Cox model ----------------------------------
  covariate_vars = c("age", "sex", "egfr", "systolic_bp",
                      "smoking_habit", "drinking_habit", "exercise_habit",
                      "hdl_cholesterol", "ldl_cholesterol"),

  # ---- Reference files (bare filename; resolved against reference_data/korea/) --
  atc_map_file = "NHIS_ATC_MAPPED.xlsx"
)
