# 04_combined_model_sweep.py
# Combined species+pathway panel: each modality's rfe_full_pool
# (species: 19 features, pathway: 18 features -- each modality's
# best-performing panel, per the completed standalone sweeps),
# CLR'd SEPARATELY (already done per-table), then concatenated --
# NOT re-CLR'd after concatenation (would wrongly treat species and
# pathway as one shared composition; they're different simplexes).
# 37 features total. All 4 models run (not assuming ElasticNet wins
# again -- combining modalities can shift which model performs best).

import argparse
import os
import warnings
from pathlib import Path

import pandas as pd
import yaml
import mlflow

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent / "utils"))
from cv_utils import clr_transform, run_nested_cv_optuna

warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser()
parser.add_argument("--model", required=True)
parser.add_argument("--smoke-test", action="store_true")
args = parser.parse_args()

PANEL_NAME = "rfe_full_pool"  # fixed -- each modality's best-performing panel

REPO_ROOT = Path(__file__).resolve().parents[2]
with open(REPO_ROOT / "config" / "config.yaml") as f:
    cfg = yaml.safe_load(f)

export_dir = REPO_ROOT / "results" / "python_export"
results_dir = REPO_ROOT / "results" / "ml_results" / "combined" / args.model
results_dir.mkdir(parents=True, exist_ok=True)

mlflow.set_tracking_uri(f"http://{cfg['mlflow']['tracking_host']}:{cfg['mlflow']['tracking_port']}")
mlflow.set_experiment(cfg["mlflow"]["experiment_name"])

meta_cols = ["sample_id", "cohort", "group"]

# ---- Species block: CLR on full table, then subset ----
species_df = pd.read_csv(export_dir / "species_full_filtered.csv")
species_feature_cols = [c for c in species_df.columns if c not in meta_cols]
species_clr, species_pc = clr_transform(species_df[species_feature_cols])
print(f"Species CLR applied to full {len(species_feature_cols)}-feature table, pseudocount={species_pc}")

species_membership = pd.read_csv(export_dir / "species_panel_membership.csv")
species_panel_features = species_membership[species_membership["panel_name"] == PANEL_NAME]["feature"].tolist()
species_subset = species_clr[species_panel_features].copy()
species_subset.columns = [f"species__{c}" for c in species_subset.columns]  # prefix to avoid any name collision
species_subset["sample_id"] = species_df["sample_id"].values

# ---- Pathway block: CLR on full table, then subset ----
pathway_df = pd.read_csv(export_dir / "pathway_full_filtered.csv")
pathway_feature_cols = [c for c in pathway_df.columns if c not in meta_cols]
pathway_clr, pathway_pc = clr_transform(pathway_df[pathway_feature_cols])
print(f"Pathway CLR applied to full {len(pathway_feature_cols)}-feature table, pseudocount={pathway_pc}")

pathway_membership = pd.read_csv(export_dir / "pathway_panel_membership.csv")
pathway_panel_features = pathway_membership[pathway_membership["panel_name"] == PANEL_NAME]["feature"].tolist()
pathway_subset = pathway_clr[pathway_panel_features].copy()
pathway_subset.columns = [f"pathway__{c}" for c in pathway_subset.columns]
pathway_subset["sample_id"] = pathway_df["sample_id"].values

# ---- Concatenate on sample_id (species=353 samples, pathway=352 --
# inner join, matching the R-side combined-panel handling from way
# back, which found exactly 1 dropped sample) ----
merged = species_subset.merge(pathway_subset, on="sample_id", how="inner")
dropped = set(species_subset["sample_id"]) - set(merged["sample_id"])
print(f"Combined merge: {merged.shape[0]} samples (dropped: {dropped})")

group_lookup = species_df.set_index("sample_id")["group"]
y = (group_lookup.loc[merged["sample_id"]] == "IBD").astype(int).reset_index(drop=True)

feature_cols_combined = [c for c in merged.columns if c != "sample_id"]
X = merged[feature_cols_combined]
print(f"Combined panel: {X.shape[1]} features ({len(species_panel_features)} species + {len(pathway_panel_features)} pathway), "
      f"{X.shape[0]} samples, {y.sum()} IBD / {(y==0).sum()} control")

run_name = f"combined_{args.model}_{PANEL_NAME}"
existing = mlflow.search_runs(
    experiment_names=[cfg["mlflow"]["experiment_name"]],
    filter_string=f"tags.mlflow.runName = '{run_name}' AND status = 'FINISHED'"
)
if len(existing) > 0 and not args.smoke_test:
    print(f"Run '{run_name}' already completed in MLflow (run_id={existing.iloc[0]['run_id']}) -- skipping.")
    row = existing.iloc[0]
    summary_row = pd.DataFrame([{
        "run_name": run_name, "model": args.model, "panel": PANEL_NAME,
        "n_features": X.shape[1],
        "auroc_mean": row.get("metrics.auroc_mean"), "auroc_std": row.get("metrics.auroc_std"),
        "auprc_mean": row.get("metrics.auprc_mean"),
    }])
    summary_row.to_csv(results_dir / f"{PANEL_NAME}_summary.csv", index=False)
    Path(results_dir / f"{PANEL_NAME}.done").touch()
    sys.stdout.flush()
    os._exit(0)

n_repeats = 1 if args.smoke_test else cfg["nested_cv"]["n_repeats"]
n_splits = 3 if args.smoke_test else cfg["nested_cv"]["n_splits"]
n_trials = 3 if args.smoke_test else cfg["nested_cv"]["n_trials"]

per_repeat, fold_results, best_params = run_nested_cv_optuna(
    X, y, model_name=args.model,
    n_repeats=n_repeats, n_splits=n_splits,
    n_inner_splits=cfg["nested_cv"]["n_inner_splits"], n_trials=n_trials,
    random_state=cfg["nested_cv"]["random_state"],
    run_name=run_name + ("_SMOKETEST" if args.smoke_test else ""),
    extra_params={"feature_type": "combined", "panel": PANEL_NAME, "n_features": X.shape[1],
                  "n_species_features": len(species_panel_features), "n_pathway_features": len(pathway_panel_features),
                  "smoke_test": args.smoke_test},
)

print(f"\n{run_name}: AUROC = {per_repeat['auroc'].mean():.4f} +/- {per_repeat['auroc'].std():.4f}")

summary_row = pd.DataFrame([{
    "run_name": run_name, "model": args.model, "panel": PANEL_NAME, "n_features": X.shape[1],
    "auroc_mean": per_repeat["auroc"].mean(), "auroc_std": per_repeat["auroc"].std(),
    "auprc_mean": per_repeat["auprc"].mean(),
}])
summary_row.to_csv(results_dir / f"{PANEL_NAME}_summary.csv", index=False)

if not args.smoke_test:
    Path(results_dir / f"{PANEL_NAME}.done").touch()
print("Done.")
