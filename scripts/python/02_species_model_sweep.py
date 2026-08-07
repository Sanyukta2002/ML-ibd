# 02_species_model_sweep.py
# One model x one panel, per invocation. Called by Snakemake with
# --model and --panel args; each combo is its own rule instance with
# its own output file (resume: Snakemake only reruns missing combos).

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
parser.add_argument("--smoke-test", action="store_true",
                     help="reduced trials/repeats for a quick correctness check, not a real result")
args = parser.parse_args()

REPO_ROOT = Path(__file__).resolve().parents[2]
with open(REPO_ROOT / "config" / "config.yaml") as f:
    cfg = yaml.safe_load(f)

export_dir = REPO_ROOT / "results" / "python_export"
results_dir = REPO_ROOT / "results" / "ml_results" / "species" / args.model
results_dir.mkdir(parents=True, exist_ok=True)

mlflow.set_tracking_uri(f"http://{cfg['mlflow']['tracking_host']}:{cfg['mlflow']['tracking_port']}")
mlflow.set_experiment(cfg["mlflow"]["experiment_name"])

# ---- Load full filtered table + panel membership, CLR on the FULL
# table first, then subset to this panel (sub-compositional coherence) ----
species_df = pd.read_csv(export_dir / "species_full_filtered.csv")
meta_cols = ["sample_id", "cohort", "group"]
feature_cols = [c for c in species_df.columns if c not in meta_cols]

clr_df, pseudocount_used = clr_transform(species_df[feature_cols])
print(f"CLR applied to full {len(feature_cols)}-species table, pseudocount={pseudocount_used}")

membership = pd.read_csv(export_dir / "species_panel_membership.csv")
panel_features = membership[membership["panel_name"] == args.panel]["feature"].tolist()
assert len(panel_features) > 0, f"Panel '{args.panel}' not found in membership file"
assert all(f in clr_df.columns for f in panel_features), "Panel feature missing from CLR table"

X = clr_df[panel_features]
y = (species_df["group"] == "IBD").astype(int)
print(f"Panel '{args.panel}': {X.shape[1]} features, {X.shape[0]} samples, {y.sum()} IBD / {(y==0).sum()} control")

# ---- Resume check: skip if this exact combo already has a completed MLflow run ----
run_name = f"species_{args.model}_{args.panel}"
existing = mlflow.search_runs(
    experiment_names=[cfg["mlflow"]["experiment_name"]],
    filter_string=f"tags.mlflow.runName = '{run_name}' AND status = 'FINISHED'"
)
if len(existing) > 0 and not args.smoke_test:
    print(f"Run '{run_name}' already completed in MLflow (run_id={existing.iloc[0]['run_id']}) -- skipping.")
    row = existing.iloc[0]
    # regenerate the summary CSV too, not just the .done marker -- a real gap
    # found in practice: the resume-skip path used to only touch .done, leaving
    # _summary.csv missing (aggregation later fails loudly on this, by design)
    summary_row = pd.DataFrame([{
        "run_name": run_name, "model": args.model, "panel": args.panel,
        "n_features": len(panel_features),
        "auroc_mean": row.get("metrics.auroc_mean"), "auroc_std": row.get("metrics.auroc_std"),
        "auprc_mean": row.get("metrics.auprc_mean"),
    }])
    summary_row.to_csv(results_dir / f"{args.panel}_summary.csv", index=False)
    Path(results_dir / f"{args.panel}.done").touch()
    sys.stdout.flush()
    os._exit(0)  # bypass Python's normal shutdown/atexit machinery entirely --
                  # sys.exit(0) was observed returning exit code 1 anyway (likely
                  # an atexit hook, possibly from mlflow's HTTP client, raising
                  # during interpreter teardown and overriding the requested code).
                  # Nothing left to clean up at this point, safe to hard-exit.

n_repeats = 1 if args.smoke_test else cfg["nested_cv"]["n_repeats"]
n_splits = 3 if args.smoke_test else cfg["nested_cv"]["n_splits"]
n_trials = 3 if args.smoke_test else cfg["nested_cv"]["n_trials"]

per_repeat, fold_results, best_params = run_nested_cv_optuna(
    X, y, model_name=args.model,
    n_repeats=n_repeats, n_splits=n_splits,
    n_inner_splits=cfg["nested_cv"]["n_inner_splits"], n_trials=n_trials,
    random_state=cfg["nested_cv"]["random_state"],
    run_name=run_name + ("_SMOKETEST" if args.smoke_test else ""),
    extra_params={"feature_type": "species", "panel": args.panel, "n_features": X.shape[1],
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
