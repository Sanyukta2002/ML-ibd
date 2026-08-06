rule mlflow_concurrency_test:
    output:
        "results/checkpoints/mlflow_test_{job_id}.done"
    log: "logs/mlflow_test_{job_id}.log"
    params:
        mlflow_host = "h2.quartz.uits.iu.edu"
    resources: mem_mb = 2000, runtime = 5
    shell:
        "mamba run -p /N/project/BacInteraction/schapag_cowrumen/ibd_crosscohort/envs/python_ml "
        "python scripts/python/test/mlflow_concurrency_test.py {wildcards.job_id} {params.mlflow_host} > {log} 2>&1 && touch {output}"
