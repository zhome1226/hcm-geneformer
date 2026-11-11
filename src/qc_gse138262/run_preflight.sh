#!/usr/bin/env bash
set -euo pipefail
root="$HOME/xiaoe-internship-data"
work="$root/_processing/all_data_cleaning_20261111/gse138262"
cd "$work"
exec 8>preflight.lock
flock -n 8 || exit 75
echo "$$" > preflight_pid.txt
trap 'code=$?; echo "$code" > preflight_exit_code.txt' EXIT
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=2
bash run_download.sh > download_v2.log 2>&1
bash run_inspect.sh > inspect.log 2>&1
"$root/_processing/SCP1303_20261111/venv/bin/python" -u stage_gse138262.py > stage.log 2>&1
