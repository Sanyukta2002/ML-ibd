# 07_statistical_comparison.py
# Paired statistical comparison of models within each panel/feature-
# type, using the 5 per-repeat AUROC values (genuine paired samples --
# same CV repeats underneath each model). Wilcoxon signed-rank test
# (matches the R-side's own non-parametric convention throughout this
# project) plus paired t-test for reference.
#
# Species' MLflow runs predate the per_repeat_metrics.csv fix -- their
# per-repeat AUROC is reconstructed from fold_predictions.csv (group
# 50 folds into 5 repeats of 10, recompute AUROC per repeat) so all
# three feature types are compared on identical footing.

import ast
import warnings
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
import mlflow
from scipy.stats import wilcoxon, ttest_rel
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")

REPO_ROOT = Path(__file__).resolve().parents[2]
with open(REPO_ROOT / "config" / "config.yaml") as f:
    cfg = yaml.safe_load(f)

mlflow.set_tracking_uri(f"http://{cfg['mlflow']['tracking_host']}:{cfg['mlflow']['tracking_port']}")
client = mlflow.tracking.MlflowClient()
N_SPLITS = cfg["nested_cv"]["n_splits"]


def get_per_repeat_auroc(run_id, run_name):
    """
    Returns a length-5 array of per-repeat AUROC for a given run.
    Tries per_repeat_metrics.csv first (pathway/combined); falls back
    to reconstructing from fold_predictions.csv (species).
    """
    artifacts = client.list_artifacts(run_id)
    artifact_paths = {a.path for a in artifacts}

    if "per_repeat_metrics" in artifact_paths:
        local_path = client.download_artifacts(run_id, "per_repeat_metrics/per_repeat_metrics.csv")
        df = pd.read_csv(local_path)
        return df["auroc"].values

    elif "predictions" in artifact_paths:
        local_path = client.download_artifacts(run_id, "predictions/fold_predictions.csv")
        df = pd.read_csv(local_path)
        df["repeat_id"] = df["fold_num"] // N_SPLITS
        per_repeat = []
        for repeat_id, group in df.groupby("repeat_id"):
            y_true_all, y_proba_all = [], []
            for _, row in group.iterrows():
                y_true_all.extend(ast.literal_eval(row["y_true"]))
                y_proba_all.extend(ast.literal_eval(row["y_proba"]))
            per_repeat.append(roc_auc_score(y_true_all, y_proba_all))
        return np.array(per_repeat)

    else:
        raise ValueError(f"Run {run_name} ({run_id}) has neither per_repeat_metrics nor predictions artifacts")


# ---- Gather all real (non-smoke-test) sweep runs ----
runs = mlflow.search_runs(experiment_names=[cfg["mlflow"]["experiment_name"]])
sweep_runs = runs[
    runs["tags.mlflow.runName"].str.match(r"^(species|pathway|combined)_(random_forest|xgboost|lightgbm|elasticnet)_rfe_", na=False)
    & (runs["status"] == "FINISHED")
]
print(f"Found {len(sweep_runs)} sweep runs to compare\n")

per_repeat_cache = {}
for _, row in sweep_runs.iterrows():
    name = row["tags.mlflow.runName"]
    per_repeat_cache[name] = get_per_repeat_auroc(row["run_id"], name)
    print(f"  Loaded {name}: {per_repeat_cache[name].round(4).tolist()}")

# ---- Parse run names into (feature_type, model, panel) ----
def parse_run_name(name):
    for ft in ["species", "pathway", "combined"]:
        if name.startswith(ft + "_"):
            rest = name[len(ft) + 1:]
            for model in ["random_forest", "xgboost", "lightgbm", "elasticnet"]:
                if rest.startswith(model + "_"):
                    panel = rest[len(model) + 1:]
                    return ft, model, panel
    return None, None, None

parsed = {name: parse_run_name(name) for name in per_repeat_cache}

# ---- Pairwise comparisons: within each (feature_type, panel), compare all model pairs ----
comparison_rows = []
groups = {}
for name, (ft, model, panel) in parsed.items():
    if ft is None:
        continue
    groups.setdefault((ft, panel), {})[model] = name

for (ft, panel), model_runs in groups.items():
    models_here = list(model_runs.keys())
    for model_a, model_b in combinations(models_here, 2):
        vals_a = per_repeat_cache[model_runs[model_a]]
        vals_b = per_repeat_cache[model_runs[model_b]]

        try:
            wstat, wp = wilcoxon(vals_a, vals_b)
        except ValueError:
            wstat, wp = np.nan, np.nan  # identical values -- wilcoxon undefined
        tstat, tp = ttest_rel(vals_a, vals_b)

        comparison_rows.append({
            "feature_type": ft, "panel": panel,
            "model_a": model_a, "model_b": model_b,
            "mean_a": vals_a.mean(), "mean_b": vals_b.mean(), "mean_diff": vals_a.mean() - vals_b.mean(),
            "wilcoxon_stat": wstat, "wilcoxon_p": wp,
            "paired_t_stat": tstat, "paired_t_p": tp,
        })

comparison_df = pd.DataFrame(comparison_rows).sort_values(["feature_type", "panel", "wilcoxon_p"])
# NOTE: Wilcoxon signed-rank with n=5 (our CV repeats) has a hard
# floor of p=0.0625 for its smallest achievable two-sided p-value --
# it can NEVER reach p<0.05 at this sample size, regardless of effect
# strength. Flagged here rather than silently used as the significance
# criterion (an earlier version of this script did exactly that and
# found "0 significant results" -- not a real finding, a test-choice
# error). Paired t-test (continuous) is the primary significance
# criterion at this sample size; Wilcoxon is still reported for
# reference/robustness but not used to gate "significant".
comparison_df["wilcoxon_underpowered_at_n5"] = True  # always true here, documented above
comparison_df["significant_paired_t_p05"] = comparison_df["paired_t_p"] < 0.05

pd.set_option("display.width", 160)
print("\n" + "=" * 100)
print("PAIRWISE MODEL COMPARISONS (within each feature_type x panel)")
print("=" * 100)
print(comparison_df.to_string(index=False))

tables_dir = REPO_ROOT / "results" / "tables"
comparison_df.to_csv(tables_dir / "statistical_comparison_pairwise.csv", index=False)
print(f"\nSaved to {tables_dir / 'statistical_comparison_pairwise.csv'}")

# ---- Specifically: ElasticNet vs. each tree model, per feature_type,
# averaged across panels -- the actual headline claim being tested ----
print("\n" + "=" * 100)
print("ELASTICNET vs. TREE MODELS -- significant differences only (p < 0.05, Wilcoxon)")
print("=" * 100)
en_comparisons = comparison_df[
    (comparison_df["model_a"] == "elasticnet") | (comparison_df["model_b"] == "elasticnet")
]
sig = en_comparisons[en_comparisons["significant_paired_t_p05"]]
print(sig[["feature_type", "panel", "model_a", "model_b", "mean_a", "mean_b", "paired_t_p"]].to_string(index=False))
print(f"\n{len(sig)} of {len(en_comparisons)} ElasticNet-vs-tree comparisons are significant at p<0.05 (paired t-test)")
print("(Wilcoxon signed-rank cannot reach p<0.05 at n=5 CV repeats -- floor is p=0.0625 -- reported for reference only, not used as the significance criterion)")
