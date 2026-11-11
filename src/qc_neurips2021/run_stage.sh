#!/usr/bin/env bash
set -euo pipefail
root="$HOME/xiaoe-internship-data"
work="$root/_processing/all_data_cleaning_20261111/NeurIPS2021"
cd "$work"
exec 9>stage.lock
flock -n 9 || exit 75
echo "$$" > stage_pid.txt
trap 'code=$?; echo "$code" > stage_exit_code.txt' EXIT
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
"$root/_processing/SCP1303_20261111/venv/bin/python" -u stage_inspect.py
