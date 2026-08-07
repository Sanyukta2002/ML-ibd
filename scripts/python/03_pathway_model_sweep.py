# 03_pathway_model_sweep.py
# Mirrors 02_species_model_sweep.py exactly, for pathway data.
# Fixes baked in from species' real bugs: nested directory output
# paths (not concatenated filenames), resume-skip writes the
# summary CSV too (not just .done), direct python binary (no
# mamba run lock contention under concurrent SLURM jobs).

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
parser.add_argument("--panel", required=True)
parser.add_argument("--smoke-test", action="store_true")
args = parser.parse_args()

REPO_ROOT = Path(__file__).resolve().parents[2]
with open(REPO_ROOT / "config" / "config.yaml") as f:
    cfg = yaml.safe_load(f)

export_dir = REPO_ROOT / "results" / "python_export"
results_dir = REPO_ROOT / "results" / "ml_results" / "pathway" / args.model
results_dir.mkdir(parents=True, exist_ok=True)

mlflow.set_tracking_uri(f"http://{cfg['mlflow']['tracking_host']}:{cfg['mlflow']['tracking_port']}")
mlflow.set_experiment(cfg["mlflow"]["experiment_name"])

# ---- CLR on the FULL filtered pathway table first, subset to panel after ----
pathway_df = pd.read_csv(export_dir / "pathway_full_filtered.csv")
meta_cols = ["sample_id", "cohort", "group"]
feature_cols = [c for c in pathway_df.columns if c not in meta_cols]

clr_df, pseudocount_used = clr_transform(pathway_df[feature_cols])
print(f"CLR applied to full {len(feature_cols)}-pathway table, pseudocount={pseudocount_used}")

membership = pd.read_csv(export_dir / "pathway_panel_membership.csv")
panel_features = membership[membership["panel_name"] == args.panel]["feature"].tolist()
assert len(panel_features) > 0, f"Panel '{args.panel}' not found in membership file"
assert all(f in clr_df.columns for f in panel_features), "Panel feature missing from CLR table"

X = clr_df[panel_features]
y = (pathway_df["group"] == "IBD").astype(int)
print(f"Panel '{args.panel}': {X.shape[1]} features, {X.shape[0]} samples, {y.sum()} IBD / {(y==0).sum()} control")

run_name = f"pathway_{args.model}_{args.panel}"
existing = mlflow.search_runs(
    experiment_names=[cfg["mlflow"]["experiment_name"]],
    filter_string=f"tags.mlflow.runName = '{run_name}' AND status = 'FINISHED'"
)
if len(existing) > 0 and not args.smoke_test:
    print(f"Run '{run_name}' already completed in MLflow (run_id={existing.iloc[0]['run_id']}) -- skipping.")
    row = existing.iloc[0]
    summary_row = pd.DataFrame([{
        "run_name": run_name, "model": args.model, "panel": args.panel,
        "n_features": len(panel_features),
        "auroc_mean": row.get("metrics.auroc_mean"), "auroc_std": row.get("metrics.auroc_std"),
        "auprc_mean": row.get("metrics.auprc_mean"),
    }])
    summary_row.to_csv(results_dir / f"{args.panel}_summary.csv", index=False)
    Path(results_dir / f"{args.panel}.done").touch()
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
    extra_params={"feature_type": "pathway", "panel": args.panel, "n_features": X.shape[1],
                  "smoke_test": args.smoke_test},
)

print(f"\n{run_name}: AUROC = {per_repeat['auroc'].mean():.4f} +/- {per_repeat['auroc'].std():.4f}")

summary_row = pd.DataFrame([{
    "run_name": run_name, "model": args.model, "panel": args.panel, "n_features": X.shape[1],
    "auroc_mean": per_repeat["auroc"].mean(), "auroc_std": per_repeat["auroc"].std(),
    "auprc_mean": per_repeat["auprc"].mean(),
}])
summary_row.to_csv(results_dir / f"{args.panel}_summary.csv", index=False)

if not args.smoke_test:
    Path(results_dir / f"{args.panel}.done").touch()
print("Done.")
