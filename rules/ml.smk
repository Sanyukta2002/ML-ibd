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
        "PYTHONUNBUFFERED=1 PYTHONWARNINGS=ignore /N/project/BacInteraction/schapag_cowrumen/ibd_crosscohort/envs/python_ml/bin/python "
        "scripts/python/02_species_model_sweep.py --model {wildcards.model} --panel {wildcards.panel} "
        "> {log} 2>&1"

rule species_model_sweep_all:
    input:
        expand("results/ml_results/species/{model}/{panel}.done",
               model=SPECIES_MODELS, panel=SPECIES_PANELS)

# ============================================================
# Pathway-only model sweep (Phase 2). Same structure as species,
# fixes already baked in.
# ============================================================

PATHWAY_MODELS = ["random_forest", "xgboost", "lightgbm", "elasticnet"]
PATHWAY_PANELS = ["rfe_optimum", "rfe_parsimony", "rfe_full_pool"]

rule pathway_model_sweep:
    input:
        pathway_csv = f"{config['paths']['python_export_dir']}/pathway_full_filtered.csv",
        membership = f"{config['paths']['python_export_dir']}/pathway_panel_membership.csv",
        config = "config/config.yaml"
    output:
        done = "results/ml_results/pathway/{model}/{panel}.done",
        summary = "results/ml_results/pathway/{model}/{panel}_summary.csv"
    log:
        "logs/ml_pathway_{model}_{panel}.log"
    resources:
        mem_mb = 2000,
        cpus_per_task = 4,
        runtime = 120
    shell:
        "SLURM_CPUS_PER_TASK={resources.cpus_per_task} "
        "PYTHONUNBUFFERED=1 PYTHONWARNINGS=ignore /N/project/BacInteraction/schapag_cowrumen/ibd_crosscohort/envs/python_ml/bin/python "
        "scripts/python/03_pathway_model_sweep.py --model {wildcards.model} --panel {wildcards.panel} "
        "> {log} 2>&1"

rule pathway_model_sweep_all:
    input:
        expand("results/ml_results/pathway/{model}/{panel}.done",
               model=PATHWAY_MODELS, panel=PATHWAY_PANELS)

# ============================================================
# Combined species+pathway sweep (Phase 3). Fixed at each
# modality's rfe_full_pool panel (each modality's actual best
# performer, per the completed standalone sweeps) -- 37 features
# total. All 4 models, deliberately not assuming ElasticNet wins
# again with a different feature-space structure.
# ============================================================

COMBINED_MODELS = ["random_forest", "xgboost", "lightgbm", "elasticnet"]

rule combined_model_sweep:
    input:
        species_csv = f"{config['paths']['python_export_dir']}/species_full_filtered.csv",
        pathway_csv = f"{config['paths']['python_export_dir']}/pathway_full_filtered.csv",
        species_membership = f"{config['paths']['python_export_dir']}/species_panel_membership.csv",
        pathway_membership = f"{config['paths']['python_export_dir']}/pathway_panel_membership.csv",
        config = "config/config.yaml"
    output:
        done = "results/ml_results/combined/{model}/rfe_full_pool.done",
        summary = "results/ml_results/combined/{model}/rfe_full_pool_summary.csv"
    log:
        "logs/ml_combined_{model}.log"
    resources:
        mem_mb = 2000,
        cpus_per_task = 4,
        runtime = 150   # 37 features vs species' max 19 -- padding above species' RF ceiling (~100min)
    shell:
        "SLURM_CPUS_PER_TASK={resources.cpus_per_task} PYTHONUNBUFFERED=1 PYTHONWARNINGS=ignore "
        "/N/project/BacInteraction/schapag_cowrumen/ibd_crosscohort/envs/python_ml/bin/python "
        "scripts/python/04_combined_model_sweep.py --model {wildcards.model} "
        "> {log} 2>&1"

rule combined_model_sweep_all:
    input:
        expand("results/ml_results/combined/{model}/rfe_full_pool.done", model=COMBINED_MODELS)

# ============================================================
# LOCO (leave-one-cohort-out) validation (Phase 4). Tests the
# 3 winning ElasticNet configs against genuine cross-cohort
# generalization -- directly stress-tests the PERMANOVA finding
# that cohort explains far more variance (R^2~0.13) than IBD
# status (R^2~0.01), which pooled nested-CV can't rule out.
# ============================================================

LOCO_FEATURE_TYPES = ["species", "pathway", "combined"]

rule loco_validation:
    input:
        species_csv = f"{config['paths']['python_export_dir']}/species_full_filtered.csv",
        pathway_csv = f"{config['paths']['python_export_dir']}/pathway_full_filtered.csv",
        config = "config/config.yaml"
    output:
        summary = "results/loco/{feature_type}/loco_summary.csv"
    log:
        "logs/loco_{feature_type}.log"
    resources:
        mem_mb = 2000,
        cpus_per_task = 4,
        runtime = 30
    shell:
        "SLURM_CPUS_PER_TASK={resources.cpus_per_task} PYTHONUNBUFFERED=1 PYTHONWARNINGS=ignore "
        "/N/project/BacInteraction/schapag_cowrumen/ibd_crosscohort/envs/python_ml/bin/python "
        "scripts/python/06_loco_validation.py --feature-type {wildcards.feature_type} "
        "> {log} 2>&1"

rule loco_validation_all:
    input:
        expand("results/loco/{feature_type}/loco_summary.csv", feature_type=LOCO_FEATURE_TYPES)
