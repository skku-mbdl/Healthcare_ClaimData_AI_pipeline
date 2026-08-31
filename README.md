# Healthcare_ClaimData_AI_pipeline
Automatic pipeline to develop AI prediction model using medical claim data

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
Rscript 00_setup.R                       # once per machine
Rscript countries/korea/run_pipeline.R   # Stage 1 -> Stage 2, Korea
Rscript countries/japan/run_pipeline.R   # Stage 1 -> Stage 2, Japan
```

Each run writes to `countries/<country>/results/`, `logs/`, and `figures/`
(all gitignored). See [CLAUDE.md](CLAUDE.md) for the full architecture,
methodology notes, and known caveats -- read it before changing anything
under `core/`.

## What's NOT included

`CKM_PREVENT/JAPAN/PREVENT` (the AHA PREVENT risk-equation validation
project) is deliberately excluded -- it answers a different research
question (external-validating a published risk score) and is not part of
this drug-side-effect pipeline. See CLAUDE.md for the one place this
exclusion mattered architecturally (Japan's Stage 2 cohort source).

## Status

Written and reviewed, with the shared statistical/data-transform logic in
`core/R/` unit-tested against synthetic data (see CLAUDE.md's "Testing"
section) -- but **never run against real NHIS/JMDC data or a live CKM_DRUG
database**, since neither is reachable from the environment this was
built in. Treat every script as reviewed-but-unverified until run for real.
