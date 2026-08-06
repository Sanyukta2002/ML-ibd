#!/bin/bash
# ============================================================
# start_mlflow_server.sh
# Starts the MLflow tracking server as a standing background
# process. Must be running before any Python ML rule that logs
# to MLflow executes -- rules connect to it over HTTP, they do
# NOT open the backend database file directly (see commit history
# for why: file-store is deprecated/unsafe, bare SQLite hits
# "database is locked" under concurrent SLURM jobs).
#
# Run this ONCE per session from the login node, inside tmux, so
# it survives you disconnecting. Check status with:
#   ps aux | grep "mlflow server" | grep -v grep
# Stop it with:
#   pkill -f "mlflow server"
# ============================================================

set -e
cd "$(dirname "$0")/.."   # repo root, regardless of where this is called from
REPO_ROOT=$(pwd)

HOST="h2.quartz.uits.iu.edu"   # must match config.yaml's mlflow.tracking_host
PORT=5050

mkdir -p logs mlartifacts

echo "Starting MLflow tracking server..."
echo "  Backend store: sqlite:///${REPO_ROOT}/mlflow.db"
echo "  Artifact root: ${REPO_ROOT}/mlartifacts"
echo "  Listening on: ${HOST}:${PORT}"

nohup mamba run -p "${REPO_ROOT}/envs/python_ml" \
  mlflow server \
  --backend-store-uri "sqlite:///${REPO_ROOT}/mlflow.db" \
  --default-artifact-root "${REPO_ROOT}/mlartifacts" \
  --host 0.0.0.0 --port "${PORT}" \
  --allowed-hosts "${HOST}:${PORT},localhost:${PORT}" \
  > "${REPO_ROOT}/logs/mlflow_server.log" 2>&1 &

sleep 5
if ps aux | grep "mlflow server" | grep -v grep > /dev/null; then
  echo "Server started successfully. PID(s):"
  ps aux | grep "mlflow server" | grep -v grep | awk '{print $2}'
  echo "Verify with: curl -s http://${HOST}:${PORT}/health"
else
  echo "WARNING: server process not found after startup -- check logs/mlflow_server.log"
  exit 1
fi
