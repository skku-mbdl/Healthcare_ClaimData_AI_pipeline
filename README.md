# Healthcare_Data_AI_pipeline

One integrated, automatic pharmacoepidemiology pipeline for Korea (NHIS Sample
Cohort) and Japan (JMDC claims), merging what used to be two independent
projects per country:

- **Stage 1 -- drug screen**: build the CKM (Cardiovascular-Kidney-Metabolic)
  stage 2/3 cohort from raw claims, then run a single-drug-at-a-time Cox
  screen over every prescribed ATC code to flag candidate cardiovascular
  side-effect drugs.
- **Stage 2 -- side-effect model**: take Stage 1's candidate drugs, run
  Boruta multivariate selection, fit cause-specific Cox models per outcome
  (MACE/ASCVD/HF/CV death), validate on a held-out test set, and explain the
  model with SHAP.

Previously, Stage 1 and Stage 2 were separate projects with a **manual** file
hand-off (a candidate-drug CSV copied by hand between them). Here, both
stages run in one command per country, and the hand-off is automatic.

## Running it

```r
Rscript 00_setup.R                              # once per machine
Rscript countries/korea/run_pipeline.R          # Stage 1 -> Stage 2, Korea, real NHIS/CKM_DRUG
Rscript countries/japan/run_pipeline.R          # Stage 1 -> Stage 2, Japan, real JMDC/CKM_DRUG

Rscript countries/korea/run_pipeline.R --demo   # same, but against synthetic data -- no DB needed
Rscript countries/japan/run_pipeline.R --demo   # (see data/metadata.json, data/generate_synthetic_data.R)
```

Each run writes to `countries/<country>/results/`, `logs/`, and `figures/`
(all gitignored). 

## What's NOT included

`CKM_PREVENT/JAPAN/PREVENT` (the AHA PREVENT risk-equation validation
project) is deliberately excluded -- it answers a different research
question (external-validating a published risk score) and is not part of
this drug-side-effect pipeline. 

## Status

Both countries' full pipelines (`--demo` mode) have been **actually run end to
end and finished successfully** -- cohort/outcome extraction through Boruta
selection, Cox model fitting, held-out validation, SHAP explanation, and the
Boruta Excel export -- against synthetic data with genuine embedded signal
(3 of 15 simulated drugs per country carry a real excess hazard). Both
countries' models correctly recovered that signal: Korea confirmed all 3 of
its true "risky" drugs by name across every outcome; Japan confirmed 1 of 3,
also as SHAP's top driver. That pilot run caught and fixed four real bugs
(two introduced by the merge, two pre-existing in the original pipelines but
never previously exercised) -- 

What's still unverified: any actual SQL query against real NHIS/JMDC/
CKM_DRUG, since neither a SQL Server nor real patient data is reachable from
the environment this was built in. Treat the DB-facing SQL in each country's
`stage1_drug_screen`/`stage2_sideeffect_model` adapters as
reviewed-but-unverified until run against a live server -- everything else
has now actually been executed, not just reviewed.
