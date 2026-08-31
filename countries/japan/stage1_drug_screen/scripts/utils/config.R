# Central configuration for the CKM drug-candidate pipeline (Japan / JMDC).
# All study-definition constants live here so downstream scripts never
# hardcode a code list twice. Ported unchanged from CKM_Drug/Japan/scripts/
# utils/config.R -- see that project's history for the schema-mapping
# decisions clarified with the study lead (not literal from the design doc).
#
# results_dir/logs_dir are NOT set here (unlike the original) -- path
# resolution (results/, logs/) is now owned by countries/japan/run_pipeline.R.

cfg <- list(

  # ---- Study window ---------------------------------------------------
  index_date_min = as.Date("2006-01-01"),
  index_date_max = as.Date("2021-12-31"),
  lookback_days  = 365,
  observation_end_cutoff = as.Date("2022-07-31"),

  # ---- Claim types treated as "admission" -----------------------------
  # DPC (Diagnosis Procedure Combination) is Japan's bundled inpatient
  # payment scheme, confirmed with study lead to count as admission
  # alongside plain Inpatient claims.
  admission_claim_types = c("Inpatient", "DPC"),
  outpatient_claim_types = c("Outpatient"),

  # ---- Outcome ICD-10 codes (admission claims only) -------------------
  outcome_icd10 = list(
    stroke = c("I61", "I62", "I63"),
    heart_failure = c("I50"),
    chd = c("I20", "I21", "I22", "I23", "I24", "I25")
    # cardiovascular death is not code-based; see 05_outcomes.R, it is
    # Diagnosis.outcome_death_flag == 1 joined to an icd10_level1_code
    # of "I" (any cardiovascular chapter code), per design doc.
  ),

  # Washout excludes anyone with prior history of ANY of the primary
  # outcomes (any claim type, not admission-restricted -- see
  # 02_build_cohort.R). Cardiovascular death is deliberately NOT included
  # here: a CV-death record can't precede a living member's own index date.
  washout_icd10 = c("I61", "I62", "I63", "I50",
                     "I20", "I21", "I22", "I23", "I24", "I25"),

  # ---- Comorbidity ICD-10 codes ----------------------------------------
  # Rule (per design doc): >=2 outpatient claims OR >=1 inpatient/DPC
  # claim with a matching code, ascertained during the lookback period.
  #
  # NOTE: dyslipidemia is E78 (disorders of lipoprotein metabolism and
  # other lipidemias), not D78 -- D78 is an unrelated code and produced a
  # 0-member dyslipidemia count against the real data, which is what
  # caught this.
  #
  # NOTE: "I14" is dropped -- it's not an assigned ICD-10 category (the
  # hypertensive-diseases block is I10-I13 and I15; I14 is an unused gap
  # in the classification).
  comorbidity_icd10 = list(
    hypertension = c("I10", "I11", "I12", "I13", "I15"),
    diabetes = c("E10", "E11", "E12", "E13", "E14"),
    dyslipidemia = c("E78"),
    ckd = c("N18", "N19")
  ),

  # ---- Medication flags (HTN_MED / DM_MED), claims-based ---------------
  # WHO ATC level-2 prefixes matched against Drug_master.who_atc_code
  # during the lookback period. Confirmed with study lead: claims-based,
  # not the checkup self-report fields.
  htn_med_atc_prefix = c("C02", "C03", "C07", "C08", "C09"),
  dm_med_atc_prefix  = c("A10"),

  # ---- CKM stage 2/3 cohort thresholds ----------------------------------
  # Field mapping from design-doc SAS names to actual JMDC schema:
  #   SEX_TYPE    -> gender_of_member ("Male"/"Female")
  #   G1E_WSTC    -> abdominal_circumference
  #   G1E_HDL     -> hdl_cholesterol
  #   G1E_TG      -> triglyceride
  #   G1E_BP_SYS  -> systolic_bp
  #   G1E_BP_DIA  -> diastolic_bp
  #   G1E_FBS     -> fasting_blood_sugar
  #   G1E_CRTN    -> serum_creatinine
  #   G1E_URN_PROT-> uric_protein_qualitative
  #   G1E_GFR     -> derived (see core/R/clinical_calcs.R's egfr())
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
  # uric_protein_qualitative only has codes 1-5 (-, +/-, +, ++, +++) in
  # this JMDC extract; design doc's "grade in (3,4,5,6)" is mapped to the
  # top three existing severity codes (+, ++, +++), confirmed with study
  # lead.
  ckm_urine_protein_grades = c(3, 4, 5),

  # ---- Drug exposure / PDC ----------------------------------------------
  # Exposure is evaluated per WHO ATC code (Drug_master.who_atc_code),
  # confirmed with study lead -- not per specific jmdc_drug_code and not
  # aggregated by ingredient (general_name).
  pdc_threshold = 0.8,
  min_exposed_n = 10,

  # ---- Candidate drug screening rule ------------------------------------
  candidate_hr_min = 1.00,
  candidate_p_max = 0.05,

  # ---- Covariates used in the Cox model ----------------------------------
  covariate_vars = c("age", "sex", "egfr", "systolic_bp",
                      "smoking_habit", "drinking_habit", "exercise_habit",
                      "hdl_cholesterol", "ldl_cholesterol")
)
