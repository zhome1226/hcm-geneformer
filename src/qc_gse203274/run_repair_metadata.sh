#!/usr/bin/env bash
set -euo pipefail
root="$HOME/xiaoe-internship-data"
work="$root/_processing/all_data_cleaning_20261111/GSE203274"
cd "$work"
exec 9>repair_metadata.lock
flock -n 9 || exit 75
echo "$$" > repair_metadata_pid.txt
trap 'code=$?; echo "$code" > repair_metadata_exit_code.txt' EXIT
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
"$root/_processing/SCP1303_20261111/venv/bin/python" -u repair_stage_metadata.py
