# ============================================================
# rules/mlflow_test.smk
# One-off validation rule used to confirm the MLflow tracking
# server correctly handles genuinely concurrent SLURM jobs (see
# commit history). Not part of rule all / the main pipeline --
# kept for reference / re-validation if the server setup ever
# changes.
# ============================================================

rule mlflow_concurrency_test:
    output:
        "results/checkpoints/mlflow_test_{job_id}.done"
    log: "logs/mlflow_test_{job_id}.log"
    params:
        mlflow_host = "h2.quartz.uits.iu.edu"
    resources: mem_mb = 2000, runtime = 5
    shell:
        "{config[python_env_path]}/bin/python "
        "scripts/python/test/mlflow_concurrency_test.py {wildcards.job_id} {params.mlflow_host} > {log} 2>&1 && touch {output}"
