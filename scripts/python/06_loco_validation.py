# 06_loco_validation.py
# Leave-one-cohort-out validation on the 3 winning ElasticNet
# configs (species/pathway/combined, each modality's best panel per
# the completed sweeps). Directly tests what our own PERMANOVA
# result warned about: cohort explains far more variance (R^2~0.13)
# than IBD status (R^2~0.01) in this data, so pooled nested-CV can't
# rule out the model leaning on cohort-specific signal. LOCO can:
# train on 3 cohorts, test entirely on the 4th, held-out cohort.

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
import optuna
import mlflow

from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import roc_auc_score, average_precision_score

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent / "utils"))
from cv_utils import clr_transform

warnings.filterwarnings("ignore")
optuna.logging.set_verbosity(optuna.logging.WARNING)

parser = argparse.ArgumentParser()
parser.add_argument("--feature-type", required=True, choices=["species", "pathway", "combined"])
args = parser.parse_args()

REPO_ROOT = Path(__file__).resolve().parents[2]
with open(REPO_ROOT / "config" / "config.yaml") as f:
    cfg = yaml.safe_load(f)

export_dir = REPO_ROOT / "results" / "python_export"
results_dir = REPO_ROOT / "results" / "loco" / args.feature_type
results_dir.mkdir(parents=True, exist_ok=True)

mlflow.set_tracking_uri(f"http://{cfg['mlflow']['tracking_host']}:{cfg['mlflow']['tracking_port']}")
mlflow.set_experiment(cfg["mlflow"]["experiment_name"])

meta_cols = ["sample_id", "cohort", "group"]
PANEL_NAME = "rfe_full_pool"  # winning panel in every modality per the completed sweeps

# ---- Build the feature matrix for the requested modality (mirrors
# 02/03/04_*_model_sweep.py exactly) ----
if args.feature_type in ("species", "combined"):
    species_df = pd.read_csv(export_dir / "species_full_filtered.csv")
    species_feature_cols = [c for c in species_df.columns if c not in meta_cols]
    species_clr, _ = clr_transform(species_df[species_feature_cols])
    species_membership = pd.read_csv(export_dir / "species_panel_membership.csv")
    species_panel_features = species_membership[species_membership["panel_name"] == PANEL_NAME]["feature"].tolist()

if args.feature_type in ("pathway", "combined"):
    pathway_df = pd.read_csv(export_dir / "pathway_full_filtered.csv")
    pathway_feature_cols = [c for c in pathway_df.columns if c not in meta_cols]
    pathway_clr, _ = clr_transform(pathway_df[pathway_feature_cols])
    pathway_membership = pd.read_csv(export_dir / "pathway_panel_membership.csv")
    pathway_panel_features = pathway_membership[pathway_membership["panel_name"] == PANEL_NAME]["feature"].tolist()

if args.feature_type == "species":
    X_full = species_clr[species_panel_features].copy()
    X_full["sample_id"] = species_df["sample_id"].values
    meta_full = species_df[meta_cols]
elif args.feature_type == "pathway":
    X_full = pathway_clr[pathway_panel_features].copy()
    X_full["sample_id"] = pathway_df["sample_id"].values
    meta_full = pathway_df[meta_cols]
else:  # combined
    sp = species_clr[species_panel_features].copy()
    sp.columns = [f"species__{c}" for c in sp.columns]
    sp["sample_id"] = species_df["sample_id"].values
    pw = pathway_clr[pathway_panel_features].copy()
    pw.columns = [f"pathway__{c}" for c in pw.columns]
    pw["sample_id"] = pathway_df["sample_id"].values
    X_full = sp.merge(pw, on="sample_id", how="inner")
    meta_full = species_df[species_df["sample_id"].isin(X_full["sample_id"])][meta_cols]

data = X_full.merge(meta_full, on="sample_id", how="inner")
feature_cols = [c for c in X_full.columns if c != "sample_id"]
cohorts = sorted(data["cohort"].unique())
print(f"Feature type: {args.feature_type} | {len(feature_cols)} features | cohorts: {cohorts}")

loco_results = []

for held_out_cohort in cohorts:
    train_mask = data["cohort"] != held_out_cohort
    test_mask = data["cohort"] == held_out_cohort

    X_train_raw = data.loc[train_mask, feature_cols]
    y_train = (data.loc[train_mask, "group"] == "IBD").astype(int)
    X_test_raw = data.loc[test_mask, feature_cols]
    y_test = (data.loc[test_mask, "group"] == "IBD").astype(int)

    n_train, n_test = len(y_train), len(y_test)
    n_test_ibd, n_test_control = y_test.sum(), (y_test == 0).sum()

    # scaler fit on TRAINING cohorts only -- no leakage from the held-out cohort
    scaler = StandardScaler()
    X_train = pd.DataFrame(scaler.fit_transform(X_train_raw), columns=feature_cols)
    X_test = pd.DataFrame(scaler.transform(X_test_raw), columns=feature_cols)

    def objective(trial):
        C = trial.suggest_float("C", 1e-4, 10.0, log=True)
        l1_ratio = trial.suggest_float("l1_ratio", 0.0, 1.0)
        model = LogisticRegression(C=C, l1_ratio=l1_ratio, penalty="elasticnet",
                                     solver="saga", max_iter=5000, random_state=123)
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=123)
        scores = cross_val_score(model, X_train, y_train, cv=cv, scoring="roc_auc")
        return float(scores.mean())

    study = optuna.create_study(direction="maximize", sampler=optuna.samplers.TPESampler(seed=123))
    study.optimize(objective, n_trials=20, show_progress_bar=False)

    final_model = LogisticRegression(**study.best_params, penalty="elasticnet",
                                       solver="saga", max_iter=5000, random_state=123)
    final_model.fit(X_train, y_train)

    y_proba = final_model.predict_proba(X_test)[:, 1]
    auroc = roc_auc_score(y_test, y_proba) if len(set(y_test)) > 1 else float("nan")
    auprc = average_precision_score(y_test, y_proba) if len(set(y_test)) > 1 else float("nan")

    print(f"Held out {held_out_cohort}: train n={n_train}, test n={n_test} "
          f"({n_test_ibd} IBD / {n_test_control} control) | AUROC={auroc:.4f} AUPRC={auprc:.4f}")

    loco_results.append({
        "feature_type": args.feature_type, "held_out_cohort": held_out_cohort,
        "n_train": n_train, "n_test": n_test, "n_test_ibd": int(n_test_ibd), "n_test_control": int(n_test_control),
        "auroc": auroc, "auprc": auprc, "best_C": study.best_params["C"], "best_l1_ratio": study.best_params["l1_ratio"],
    })

    with mlflow.start_run(run_name=f"loco_{args.feature_type}_holdout_{held_out_cohort}"):
        mlflow.log_param("feature_type", args.feature_type)
        mlflow.log_param("held_out_cohort", held_out_cohort)
        mlflow.log_param("n_train", n_train)
        mlflow.log_param("n_test", n_test)
        mlflow.log_metric("auroc", auroc)
        mlflow.log_metric("auprc", auprc)

loco_df = pd.DataFrame(loco_results)
print(f"\n{args.feature_type} LOCO summary:")
print(loco_df[["held_out_cohort", "n_test", "auroc", "auprc"]].to_string(index=False))
print(f"Mean AUROC across held-out cohorts: {loco_df['auroc'].mean():.4f} (vs. pooled nested-CV for comparison)")

loco_df.to_csv(results_dir / "loco_summary.csv", index=False)
print("Done.")
