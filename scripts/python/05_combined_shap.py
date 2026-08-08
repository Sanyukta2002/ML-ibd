# 05_combined_shap.py
# SHAP analysis on the combined (species+pathway) ElasticNet model --
# the winning model/panel from the ML sweep. Answers the open question
# from the sweep: does pathway carry real importance in the combined
# model, or is it effectively being ignored?
#
# Needs a SINGLE model fit on the full dataset (nested CV never
# produces this -- it only ever fits per-fold models for honest
# held-out evaluation). Standard practice: quick CV-based
# hyperparameter search on all data, then one final fit for
# interpretation, kept clearly separate from the nested-CV
# performance numbers already reported.

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
import optuna
import shap
import mlflow
import matplotlib
matplotlib.use("Agg")  # headless -- no display on the HPC
import matplotlib.pyplot as plt

from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, cross_val_score

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent / "utils"))
from cv_utils import clr_transform

warnings.filterwarnings("ignore")
optuna.logging.set_verbosity(optuna.logging.WARNING)

REPO_ROOT = Path(__file__).resolve().parents[2]
with open(REPO_ROOT / "config" / "config.yaml") as f:
    cfg = yaml.safe_load(f)

export_dir = REPO_ROOT / "results" / "python_export"
figures_dir = REPO_ROOT / "results" / "figures"
figures_dir.mkdir(parents=True, exist_ok=True)
tables_dir = REPO_ROOT / "results" / "tables"

mlflow.set_tracking_uri(f"http://{cfg['mlflow']['tracking_host']}:{cfg['mlflow']['tracking_port']}")
mlflow.set_experiment(cfg["mlflow"]["experiment_name"])
mlflow.start_run(run_name="combined_elasticnet_SHAP_interpretation")

meta_cols = ["sample_id", "cohort", "group"]
PANEL_NAME = "rfe_full_pool"

# ---- Rebuild the exact same combined panel as 04_combined_model_sweep.py ----
species_df = pd.read_csv(export_dir / "species_full_filtered.csv")
species_feature_cols = [c for c in species_df.columns if c not in meta_cols]
species_clr, _ = clr_transform(species_df[species_feature_cols])
species_membership = pd.read_csv(export_dir / "species_panel_membership.csv")
species_panel_features = species_membership[species_membership["panel_name"] == PANEL_NAME]["feature"].tolist()
species_subset = species_clr[species_panel_features].copy()
species_subset.columns = [f"species__{c}" for c in species_subset.columns]
species_subset["sample_id"] = species_df["sample_id"].values

pathway_df = pd.read_csv(export_dir / "pathway_full_filtered.csv")
pathway_feature_cols = [c for c in pathway_df.columns if c not in meta_cols]
pathway_clr, _ = clr_transform(pathway_df[pathway_feature_cols])
pathway_membership = pd.read_csv(export_dir / "pathway_panel_membership.csv")
pathway_panel_features = pathway_membership[pathway_membership["panel_name"] == PANEL_NAME]["feature"].tolist()
pathway_subset = pathway_clr[pathway_panel_features].copy()
pathway_subset.columns = [f"pathway__{c}" for c in pathway_subset.columns]
pathway_subset["sample_id"] = pathway_df["sample_id"].values

merged = species_subset.merge(pathway_subset, on="sample_id", how="inner")
group_lookup = species_df.set_index("sample_id")["group"]
y = (group_lookup.loc[merged["sample_id"]] == "IBD").astype(int).reset_index(drop=True)
feature_cols = [c for c in merged.columns if c != "sample_id"]
X = merged[feature_cols].reset_index(drop=True)

print(f"Combined panel rebuilt: {X.shape[1]} features, {X.shape[0]} samples")
assert X.shape[1] == 37 and X.shape[0] == 352, "Panel shape mismatch vs. the ML sweep -- investigate before trusting SHAP"

# ---- Quick hyperparameter search via CV on the FULL dataset (separate
# from, and simpler than, the nested-CV nested-Optuna sweep -- this is
# purely to pick reasonable params for the single interpretation fit) ----
scaler = StandardScaler()
X_scaled = pd.DataFrame(scaler.fit_transform(X), columns=X.columns, index=X.index)

def objective(trial):
    C = trial.suggest_float("C", 1e-4, 10.0, log=True)
    l1_ratio = trial.suggest_float("l1_ratio", 0.0, 1.0)
    model = LogisticRegression(C=C, l1_ratio=l1_ratio, penalty="elasticnet",
                                 solver="saga", max_iter=5000, random_state=123)
    cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=123)
    scores = cross_val_score(model, X_scaled, y, cv=cv, scoring="roc_auc")
    return float(scores.mean())

study = optuna.create_study(direction="maximize", sampler=optuna.samplers.TPESampler(seed=123))
study.optimize(objective, n_trials=30, show_progress_bar=False)
print(f"Best params for final interpretation fit: {study.best_params} (CV AUROC={study.best_value:.4f})")

# ---- Final fit on ALL data (for interpretation only -- NOT a
# performance estimate, that's what the nested-CV sweep already gave us) ----
final_model = LogisticRegression(**study.best_params, penalty="elasticnet",
                                   solver="saga", max_iter=5000, random_state=123)
final_model.fit(X_scaled, y)

# ---- SHAP (exact, via LinearExplainer -- appropriate since this is a
# genuinely linear model, not an approximation) ----
explainer = shap.LinearExplainer(final_model, X_scaled)
shap_values = explainer.shap_values(X_scaled)

mean_abs_shap = pd.DataFrame({
    "feature": X.columns,
    "mean_abs_shap": np.abs(shap_values).mean(axis=0)
}).sort_values("mean_abs_shap", ascending=False)

mean_abs_shap["feature_type"] = mean_abs_shap["feature"].apply(lambda x: "species" if x.startswith("species__") else "pathway")
mean_abs_shap.to_csv(tables_dir / "combined_shap_importance.csv", index=False)

# Save the raw SHAP values + underlying feature matrix too -- separates
# COMPUTING shap values (needs the model/data/HPC environment) from
# RENDERING figures (purely visual, can happen anywhere later, including
# locally, without rerunning anything). Beeswarm plot layout is a known
# SHAP+matplotlib rendering gotcha not worth chasing further here --
# deferred to the dedicated figures phase per project decision.
np.save(tables_dir / "combined_shap_values_raw.npy", shap_values)
X_scaled.to_csv(tables_dir / "combined_shap_input_features.csv", index=False)

# ---- The actual question: species vs pathway total contribution ----
contribution_by_type = mean_abs_shap.groupby("feature_type")["mean_abs_shap"].sum()
contribution_pct = (contribution_by_type / contribution_by_type.sum() * 100).round(1)
print("\nTotal SHAP contribution by feature type:")
print(contribution_by_type)
print("\nAs percentage:")
print(contribution_pct)

summary_df = pd.DataFrame({
    "feature_type": contribution_by_type.index,
    "total_mean_abs_shap": contribution_by_type.values,
    "pct_of_total": contribution_pct.values
})
summary_df.to_csv(tables_dir / "combined_shap_contribution_by_type.csv", index=False)

# ---- Figures ----
top20 = mean_abs_shap.head(20)
colors = ["#0072B2" if t == "species" else "#D55E00" for t in top20["feature_type"]]
clean_names = top20["feature"].str.replace("species__", "").str.replace("pathway__", "")

fig, ax = plt.subplots(figsize=(8, 7))
ax.barh(clean_names[::-1], top20["mean_abs_shap"][::-1], color=colors[::-1])
ax.set_xlabel("Mean |SHAP value|")
ax.set_title("Top 20 features -- combined species+pathway ElasticNet model")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color="#0072B2", label="Species"), Patch(color="#D55E00", label="Pathway")], loc="lower right")
plt.tight_layout()
plt.savefig(figures_dir / "combined_shap_top20.png", dpi=300)
plt.savefig(figures_dir / "combined_shap_top20.pdf")
plt.close()
X_scaled_display = X_scaled.copy()
X_scaled_display.columns = [c.replace("species__", "").replace("pathway__", "") for c in X_scaled_display.columns]
plt.figure(figsize=(10, 12))
shap.summary_plot(shap_values, X_scaled_display, show=False)
plt.tight_layout()
plt.savefig(figures_dir / "combined_shap_beeswarm.png", dpi=300, bbox_inches="tight")
plt.close()

mlflow.log_params(study.best_params)
mlflow.log_metric("interpretation_fit_cv_auroc", study.best_value)
for _, row in contribution_pct.items():
    pass
mlflow.log_metric("species_pct_shap_contribution", float(contribution_pct.get("species", 0)))
mlflow.log_metric("pathway_pct_shap_contribution", float(contribution_pct.get("pathway", 0)))
mlflow.log_artifact(str(tables_dir / "combined_shap_importance.csv"), artifact_path="shap")
mlflow.log_artifact(str(tables_dir / "combined_shap_contribution_by_type.csv"), artifact_path="shap")
mlflow.log_artifact(str(figures_dir / "combined_shap_top20.png"), artifact_path="shap")
mlflow.log_artifact(str(figures_dir / "combined_shap_beeswarm.png"), artifact_path="shap")
mlflow.end_run()

print("\nFigures saved: combined_shap_top20.{png,pdf}, combined_shap_beeswarm.png")
print("Logged to MLflow: combined_elasticnet_SHAP_interpretation")
print("Done.")
