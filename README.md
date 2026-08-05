# IBD Cross-Cohort Microbiome Classification Pipeline

Snakemake-orchestrated, R + Python pipeline classifying IBD vs.
control from gut microbiome species/pathway abundance across
4 curatedMetagenomicData cohorts (NielsenHB_2014 [ESP only],
HMP_2019_ibdmdb, HallAB_2017, IjazUZ_2017).

## Status
Scaffold only — rules being built and verified one at a time.

## Structure
- `config/config.yaml` — all thresholds/params/paths
- `config/slurm_profile/` — cluster execution config (site-specific)
- `rules/` — one .smk file per pipeline stage
- `scripts/r/`, `scripts/python/` — the actual analysis code
- `envs/` — conda env specs + Apptainer container definitions
- `results/` — outputs (mostly gitignored; `results/tables/` tracked)
