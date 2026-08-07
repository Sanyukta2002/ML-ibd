# ============================================================
# cv_utils.py
# Shared CLR transform + nested-CV + Optuna tuning + MLflow logging
# engine. Every model/panel run in this project calls
# run_nested_cv_optuna() rather than reimplementing CV/tuning/
# logging per model.
# ============================================================

import os
import time
import tempfile
import warnings

# Respect whatever SLURM actually allocated for this job's CPUs --
# never request more than granted. Falls back to 4 for interactive/
# non-SLURM use (e.g. manual testing on a login node).
N_JOBS = int(os.environ.get("SLURM_CPUS_PER_TASK", 4))

import numpy as np
import pandas as pd
import optuna
import mlflow
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from sklearn.metrics import (
    roc_auc_score, average_precision_score,
    precision_score, f1_score, accuracy_score, confusion_matrix
)
from sklearn.exceptions import ConvergenceWarning

optuna.logging.set_verbosity(optuna.logging.WARNING)  # suppress per-trial console spam

# ---- Warning suppression (per project decision -- these are known,
# expected warnings from saga solver convergence / lightgbm verbosity,
# not signals of a real problem, and were flooding the console) ----
warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")


def clr_transform(df, pseudocount_method="half_min_nonzero"):
    """
    Sample-wise (row-wise) centered log-ratio transform on a
    features-only DataFrame (samples x features, all columns numeric).

    Convention (confirmed against compositional-data-analysis
    literature and Orchestrating Microbiome Analysis guidance):
      1. pseudocount added first (half the minimum nonzero value
         across the whole table -- matches compute_gfc()'s convention
         on the R side, for consistency across the pipeline)
      2. log taken
      3. row-centered (subtract each row's own mean of logs) --
         this is what makes it "centered" log-ratio

    MUST be applied to the FULL filtered table, BEFORE subsetting to
    any panel -- CLR is not sub-compositionally coherent (the
    geometric mean must reflect the full detected community, not an
    artificially-closed panel subset). Panel subsetting happens
    AFTER calling this function, not before.

    Returns (clr_df, pseudocount_used).
    """
    if pseudocount_method != "half_min_nonzero":
        raise ValueError(f"Unsupported pseudocount_method: {pseudocount_method}")

    values = df.values.astype(float)
    nonzero_vals = values[values > 0]
    if nonzero_vals.size == 0:
        raise ValueError("No nonzero values found in the table -- cannot compute a pseudocount")
    pseudocount = nonzero_vals.min() / 2

    log_vals = np.log(values + pseudocount)
    row_means = log_vals.mean(axis=1, keepdims=True)
    clr_vals = log_vals - row_means

    clr_df = pd.DataFrame(clr_vals, index=df.index, columns=df.columns)

    # sanity check: CLR's defining property is that each row sums to ~0
    row_sums = clr_df.sum(axis=1).abs()
    assert row_sums.max() < 1e-6, f"CLR row sums not ~0 (max abs: {row_sums.max()}) -- something is wrong"

    return clr_df, pseudocount


def compute_fold_metrics(y_true, y_pred, y_proba):
    """Standard metric set for one outer fold."""
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else np.nan
    specificity = tn / (tn + fp) if (tn + fp) > 0 else np.nan
    return {
        "auroc": roc_auc_score(y_true, y_proba),
        "auprc": average_precision_score(y_true, y_proba),
        "sensitivity": sensitivity,
        "specificity": specificity,
        "precision": precision_score(y_true, y_pred, zero_division=0),
        "f1": f1_score(y_true, y_pred, zero_division=0),
        "accuracy": accuracy_score(y_true, y_pred),
    }


def aggregate_per_repeat(fold_results_df, n_splits):
    """
    Collapse fold-level metrics to per-repeat metrics (mean across
    the n_splits folds within each repeat). Matches the R-side
    discipline: never pool raw fold rows for confidence intervals
    (pseudo-replication) -- one row per genuine CV repeat instead.
    """
    df = fold_results_df.copy()
    df["repeat_id"] = df["fold_num"] // n_splits
    metric_cols = ["auroc", "auprc", "sensitivity", "specificity", "precision", "f1", "accuracy"]
    return df.groupby("repeat_id")[metric_cols].mean().reset_index()


def get_model_for_trial(trial_or_params, model_name):
    """
    Build a model instance from either an active Optuna trial
    (during search) or a plain dict of already-chosen params
    (final refit) -- same search-space definition serves both.

    ElasticNet is wrapped in a Pipeline(StandardScaler -> LogisticRegression)
    so that when sklearn's CV machinery calls .fit() on a training fold,
    the scaler is fit ONLY on that fold's training data -- fit fresh
    inside each fold automatically, no possibility of leaking test-fold
    statistics into the scaler. Tree-based models (RF/XGBoost/LightGBM)
    are scale-invariant and stay as bare estimators; wrapping them would
    be harmless but pointless.
    """
    def suggest(name, kind, *args, **kwargs):
        if isinstance(trial_or_params, dict):
            return trial_or_params[name]
        fn = getattr(trial_or_params, f"suggest_{kind}")
        return fn(name, *args, **kwargs)

    if model_name == "random_forest":
        from sklearn.ensemble import RandomForestClassifier
        return RandomForestClassifier(
            n_estimators=suggest("n_estimators", "int", 100, 800),
            max_depth=suggest("max_depth", "int", 2, 20),
            min_samples_split=suggest("min_samples_split", "int", 2, 20),
            min_samples_leaf=suggest("min_samples_leaf", "int", 1, 10),
            max_features=suggest("max_features", "categorical", ["sqrt", "log2", None]),
            random_state=123, n_jobs=N_JOBS,
        )

    elif model_name == "xgboost":
        from xgboost import XGBClassifier
        return XGBClassifier(
            n_estimators=suggest("n_estimators", "int", 100, 800),
            max_depth=suggest("max_depth", "int", 2, 10),
            learning_rate=suggest("learning_rate", "float", 0.005, 0.3, log=True),
            subsample=suggest("subsample", "float", 0.5, 1.0),
            colsample_bytree=suggest("colsample_bytree", "float", 0.5, 1.0),
            reg_alpha=suggest("reg_alpha", "float", 1e-8, 10.0, log=True),
            reg_lambda=suggest("reg_lambda", "float", 1e-8, 10.0, log=True),
            random_state=123, eval_metric="logloss", n_jobs=N_JOBS,
        )

    elif model_name == "lightgbm":
        from lightgbm import LGBMClassifier
        return LGBMClassifier(
            n_estimators=suggest("n_estimators", "int", 100, 800),
            max_depth=suggest("max_depth", "int", 2, 12),
            learning_rate=suggest("learning_rate", "float", 0.005, 0.3, log=True),
            num_leaves=suggest("num_leaves", "int", 7, 127),
            subsample=suggest("subsample", "float", 0.5, 1.0),
            colsample_bytree=suggest("colsample_bytree", "float", 0.5, 1.0),
            reg_alpha=suggest("reg_alpha", "float", 1e-8, 10.0, log=True),
            reg_lambda=suggest("reg_lambda", "float", 1e-8, 10.0, log=True),
            random_state=123, n_jobs=N_JOBS, verbose=-1,
        )

    elif model_name == "elasticnet":
        from sklearn.linear_model import LogisticRegression
        from sklearn.preprocessing import StandardScaler
        from sklearn.pipeline import Pipeline
        return Pipeline([
            ("scaler", StandardScaler()),
            ("model", LogisticRegression(
                C=suggest("C", "float", 1e-4, 10.0, log=True),
                l1_ratio=suggest("l1_ratio", "float", 0.0, 1.0),
                penalty="elasticnet", solver="saga", max_iter=5000, random_state=123,
            )),
        ])

    else:
        raise ValueError(f"Unknown model_name: {model_name}")


def run_nested_cv_optuna(X, y, model_name, n_repeats=5, n_splits=10,
                          n_inner_splits=5, n_trials=30, random_state=123,
                          log_to_mlflow=True, run_name=None, extra_params=None):
    """
    Nested CV with Optuna hyperparameter search.

    Outer loop: RepeatedStratifiedKFold(n_splits, n_repeats) - held-out
    test folds untouched by tuning.
    Inner loop: per outer-training fold, an Optuna study (TPE sampler,
    StratifiedKFold(n_inner_splits)) searches hyperparameters, optimizing
    mean inner-fold ROC AUC. Best config refit on the full outer-training
    fold, evaluated once on the untouched outer-test fold. No leakage --
    for ElasticNet, the scaler inside its Pipeline is also refit fresh
    on each outer-training fold (and each inner-training fold, during
    the search), never touching test-fold data.

    X is expected to already be CLR-transformed (see clr_transform())
    and already subset to the panel being evaluated -- this function
    does not do either of those steps itself.

    Returns:
        per_repeat_metrics : DataFrame, one row per CV repeat
        fold_results_df    : DataFrame, one row per outer fold
                              (includes y_true/y_proba, for ROC/PR curves)
        best_params_per_fold : list of dicts (logged to MLflow as an artifact)
    """
    run_start_time = time.time()

    outer_cv = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=random_state)
    X_arr = X.values if isinstance(X, pd.DataFrame) else np.asarray(X)
    y_arr = y.values if isinstance(y, pd.Series) else np.asarray(y)

    fold_records, best_params_per_fold = [], []

    for fold_num, (train_idx, test_idx) in enumerate(outer_cv.split(X_arr, y_arr)):
        X_train, X_test = X_arr[train_idx], X_arr[test_idx]
        y_train, y_test = y_arr[train_idx], y_arr[test_idx]

        def objective(trial):
            model = get_model_for_trial(trial, model_name)
            inner_cv = StratifiedKFold(n_splits=n_inner_splits, shuffle=True, random_state=random_state)
            scores = []
            for tr_idx, val_idx in inner_cv.split(X_train, y_train):
                model.fit(X_train[tr_idx], y_train[tr_idx])
                proba = model.predict_proba(X_train[val_idx])[:, 1]
                scores.append(roc_auc_score(y_train[val_idx], proba))
            return float(np.mean(scores))

        study = optuna.create_study(direction="maximize", sampler=optuna.samplers.TPESampler(seed=random_state))
        study.optimize(objective, n_trials=n_trials, show_progress_bar=False)

        best_params_per_fold.append({"fold_num": fold_num, **study.best_params})

        final_model = get_model_for_trial(study.best_params, model_name)
        final_model.fit(X_train, y_train)

        y_proba = final_model.predict_proba(X_test)[:, 1]
        y_pred = (y_proba >= 0.5).astype(int)
        metrics = compute_fold_metrics(y_test, y_pred, y_proba)

        fold_records.append({
            "fold_num": fold_num,
            "y_true": y_test.tolist(),
            "y_proba": y_proba.tolist(),
            **metrics,
        })

    elapsed_seconds = time.time() - run_start_time  # real wall-clock time for the WHOLE fold loop,
                                                       # not just the final logging block

    fold_results_df = pd.DataFrame(fold_records)
    per_repeat_metrics = aggregate_per_repeat(fold_results_df, n_splits)

    if log_to_mlflow:
        with mlflow.start_run(run_name=run_name or model_name):
            mlflow.log_param("model_name", model_name)
            mlflow.log_param("n_repeats", n_repeats)
            mlflow.log_param("n_splits", n_splits)
            mlflow.log_param("n_inner_splits", n_inner_splits)
            mlflow.log_param("n_trials", n_trials)
            if extra_params:
                for k, v in extra_params.items():
                    mlflow.log_param(k, v)

            metric_cols = ["auroc", "auprc", "sensitivity", "specificity", "precision", "f1", "accuracy"]
            for col in metric_cols:
                mlflow.log_metric(f"{col}_mean", float(per_repeat_metrics[col].mean()))
                mlflow.log_metric(f"{col}_std", float(per_repeat_metrics[col].std()))

            mlflow.log_metric("elapsed_seconds", elapsed_seconds)  # real timing, fixes original bug

            with tempfile.TemporaryDirectory() as tmp_dir:
                params_path = os.path.join(tmp_dir, "best_params_per_fold.csv")
                pd.DataFrame(best_params_per_fold).to_csv(params_path, index=False)
                mlflow.log_artifact(params_path, artifact_path="hyperparameters")

                preds_path = os.path.join(tmp_dir, "fold_predictions.csv")
                fold_results_df[["fold_num", "y_true", "y_proba"]].to_csv(preds_path, index=False)
                mlflow.log_artifact(preds_path, artifact_path="predictions")

                # Per-repeat metrics (one row per CV repeat, all 7 metrics) --
                # previously computed only in-memory to derive the mean/std,
                # then discarded. Gap found when reviewing species' results:
                # this is what a fold/repeat-spread figure (boxplot-style,
                # matching the reference papers) actually needs -- the
                # aggregate mean/std alone isn't enough to redraw that.
                # Species sweep does NOT have this (ran before the gap was
                # found) -- pathway/combined do, from this point forward.
                per_repeat_path = os.path.join(tmp_dir, "per_repeat_metrics.csv")
                per_repeat_metrics.to_csv(per_repeat_path, index=False)
                mlflow.log_artifact(per_repeat_path, artifact_path="per_repeat_metrics")

    return per_repeat_metrics, fold_results_df, best_params_per_fold
