# Reproducibility Verification

**Date:** 2026-08-11
**Method:** Clean-room test — fresh `git clone` into a separate directory
(`ibd_crosscohort_cleanroom`), both environments (R/Bioconductor container,
Python ML env) rebuilt entirely from tracked recipe files, full pipeline
(`snakemake --workflow-profile config/slurm_profile`, no target — the
complete `rule all`) run unattended from raw `curatedMetagenomicData` pull
through combined ML sweeps and LOCO validation.

## Result

- **42/42 rules completed successfully**, unattended, from a fresh clone.
- **R environment**: rebuilt via `envs/r_bioconductor.def` +
  `scripts/r/setup/install_r_packages.R` — all 15 packages installed and
  verified loadable in a single, non-interactive run (previously only ever
  built piecemeal, interactively).
- **Python environment**: rebuilt via `envs/python_ml.yml` (pinned versions)
  — every package version confirmed to match exactly (pandas 2.3.3,
  scikit-learn 1.9.0, xgboost 3.2.0, lightgbm 4.7.0, optuna 4.9.0,
  mlflow 3.15.1, shap 0.51.0).
- **Output comparison**: 6 key result tables (species/pathway/combined ML
  sweep summaries, LOCO summary, statistical comparison, SHAP contribution)
  compared byte-for-byte between the original build and the clean-room
  rebuild — **all 6 identical**.

## Bug found and fixed via this process

`results/maaslin2/` (parent output directory for all 4 MaAsLin2 rules) had
originally been created via a one-time manual `mkdir` during early project
scaffolding, never through a tracked script. Since it's an empty,
gitignored directory, `git clone` had no way to recreate it, and R's
`dir.create()` doesn't create missing parent directories by default —
causing all 4 MaAsLin2 rules to fail on a genuinely fresh clone despite
having worked throughout the entire original build. Fixed by adding
explicit `dir.create(..., recursive = TRUE)` calls to all 4 affected
scripts (see commit history). This is a concrete example of a class of
bug — "works because I set something up once by hand and forgot" — that
only a genuine clean-room test can catch; code review or normal usage of
the original directory would never have surfaced it.

## What "reproducible" means for this project

Given this pipeline runs on HPC infrastructure with a specific SLURM
allocation, full zero-edit portability to any machine isn't a meaningful
standard. The verified claim is narrower and more honest: **given the
documented preconditions (HPC access, your own SLURM account, cloning +
building the two pinned environments, one manual edit to
`config/slurm_profile/config.yaml` for site-specific values), the pipeline
produces identical results, every time.** `config/slurm_profile/config.yaml`
is the one file expected to require manual adjustment per deployment; every
other path in the codebase is relative/config-driven.
