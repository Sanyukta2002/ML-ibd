configfile: "config/config.yaml"

include: "rules/species.smk"
include: "rules/pathway.smk"
include: "rules/mlflow_test.smk"
include: "rules/ml.smk"

# ============================================================
# rule all: the single entry point for the ENTIRE pipeline --
# raw curatedMetagenomicData pull through combined/LOCO/SHAP.
# `snakemake --workflow-profile config/slurm_profile` with no
# target runs everything below, unattended, in correct dependency
# order across both R and Python stages.
#
# PRECONDITION (not part of the DAG, documented in README):
# the MLflow tracking server must be running first --
#   bash scripts/start_mlflow_server.sh
# This is a standing background service, not a per-rule dependency,
# same category as the R/Python environments themselves needing to
# be built before any rule can run.
# ============================================================

rule all:
    input:
        # species ML sweep: 4 models x 3 panels
        expand("results/ml_results/species/{model}/{panel}.done",
               model=SPECIES_MODELS, panel=SPECIES_PANELS),
        # pathway ML sweep: 4 models x 3 panels
        expand("results/ml_results/pathway/{model}/{panel}.done",
               model=PATHWAY_MODELS, panel=PATHWAY_PANELS),
        # combined (species+pathway) sweep: 4 models, 1 panel pairing
        expand("results/ml_results/combined/{model}/rfe_full_pool.done",
               model=COMBINED_MODELS),
        # LOCO validation: 3 feature types (species/pathway/combined)
        expand("results/loco/{feature_type}/loco_summary.csv",
               feature_type=LOCO_FEATURE_TYPES),
