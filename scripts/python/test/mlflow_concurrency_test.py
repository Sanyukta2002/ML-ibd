import sys, time, random
import mlflow

job_id = sys.argv[1]
mlflow_host = sys.argv[2]  # login node hostname, passed in

mlflow.set_tracking_uri(f"http://{mlflow_host}:5050")
mlflow.set_experiment("concurrency_test")

with mlflow.start_run(run_name=f"test_job_{job_id}"):
    mlflow.log_param("job_id", job_id)
    time.sleep(random.uniform(5, 15))
    mlflow.log_metric("dummy_metric", random.random())
    print(f"Job {job_id} logged successfully")
