# ============================================================
# rules/ml.smk
# Python ML stage: species/pathway/combined nested-CV sweeps +
# LOCO validation. Every rule references the Python environment
# via config['python_env_path'] (relative to repo root) rather
# than a hardcoded absolute path -- portable across clones.
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
        mem_mb = 2000,
        cpus_per_task = 4,
        runtime = 120
    shell:
        "SLURM_CPUS_PER_TASK={resources.cpus_per_task} "
        "PYTHONUNBUFFERED=1 PYTHONWARNINGS=ignore {config[python_env_path]}/bin/python "
        "scripts/python/02_species_model_sweep.py --model {wildcards.model} --panel {wildcards.panel} "
        "> {log} 2>&1"

rule species_model_sweep_all:
    input:
        expand("results/ml_results/species/{model}/{panel}.done",
               model=SPECIES_MODELS, panel=SPECIES_PANELS)

# ============================================================
# Pathway-only model sweep (Phase 2). Same structure as species.
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
        "PYTHONUNBUFFERED=1 PYTHONWARNINGS=ignore {config[python_env_path]}/bin/python "
        "scripts/python/03_pathway_model_sweep.py --model {wildcards.model} --panel {wildcards.panel} "
        "> {log} 2>&1"

rule pathway_model_sweep_all:
    input:
        expand("results/ml_results/pathway/{model}/{panel}.done",
               model=PATHWAY_MODELS, panel=PATHWAY_PANELS)

# ============================================================
# Combined species+pathway sweep (Phase 3). Fixed at each
# modality's rfe_full_pool panel -- 37 features total.
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
        runtime = 150
    shell:
        "SLURM_CPUS_PER_TASK={resources.cpus_per_task} "
        "PYTHONUNBUFFERED=1 PYTHONWARNINGS=ignore {config[python_env_path]}/bin/python "
        "scripts/python/04_combined_model_sweep.py --model {wildcards.model} "
        "> {log} 2>&1"

rule combined_model_sweep_all:
    input:
        expand("results/ml_results/combined/{model}/rfe_full_pool.done", model=COMBINED_MODELS)

# ============================================================
# LOCO (leave-one-cohort-out) validation (Phase 4).
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
        "SLURM_CPUS_PER_TASK={resources.cpus_per_task} "
        "PYTHONUNBUFFERED=1 PYTHONWARNINGS=ignore {config[python_env_path]}/bin/python "
        "scripts/python/06_loco_validation.py --feature-type {wildcards.feature_type} "
        "> {log} 2>&1"

rule loco_validation_all:
    input:
        expand("results/loco/{feature_type}/loco_summary.csv", feature_type=LOCO_FEATURE_TYPES)
