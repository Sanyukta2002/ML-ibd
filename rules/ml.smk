# ============================================================
# rules/ml.smk
# Species-only model sweep (Phase 1 of the ML build, per plan):
# 4 models x 6 panels = 24 combos, each its own rule instance/
# SLURM job -- genuinely concurrent, not sequential. Resume is
# automatic: Snakemake only reruns combos missing their .done file.
# Pathway-only and combined sweeps follow once this is validated.
# ============================================================

SPECIES_MODELS = ["random_forest", "xgboost", "lightgbm", "elasticnet"]
SPECIES_PANELS = ["rfe_optimum", "rfe_parsimony", "rfe_full_pool"]

rule species_model_sweep:
    input:
        species_csv = f"{config['paths']['python_export_dir']}/species_full_filtered.csv",
        membership = f"{config['paths']['python_export_dir']}/species_panel_membership.csv",
        config = "config/config.yaml"
    output:
        done = "results/ml_results/species/{model}/{panel}.done",
        summary = "results/ml_results/species/{model}/{panel}_summary.csv"
    log:
        "logs/ml_species_{model}_{panel}.log"
    resources:
        mem_mb = 2000,       # real usage observed: MaxRSS 401MB for RF -- padded, not the earlier
                              # unverified guess of 8000
        cpus_per_task = 4,
        runtime = 120        # real: 70min for RF via Snakemake/SLURM (properly 4-core-matched) --
                              # padded for other models, which may differ; revisit after full batch
    shell:
        "SLURM_CPUS_PER_TASK={resources.cpus_per_task} "
        "mamba run -p /N/project/BacInteraction/schapag_cowrumen/ibd_crosscohort/envs/python_ml "
        "python scripts/python/02_species_model_sweep.py --model {wildcards.model} --panel {wildcards.panel} "
        "> {log} 2>&1"

rule species_model_sweep_all:
    input:
        expand("results/ml_results/species/{model}/{panel}.done",
               model=SPECIES_MODELS, panel=SPECIES_PANELS)
