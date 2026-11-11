#!/usr/bin/env bash
set -euo pipefail
root="$HOME/xiaoe-internship-data"
work="$root/_processing/all_data_cleaning_20261111/GSE249925"
cd "$work"
exec 9>prepare.lock
flock -n 9 || exit 75
echo "$$" > prepare_pid.txt
trap 'code=$?; echo "$code" > prepare_exit_code.txt' EXIT
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
"$root/_processing/SCP1303_20261111/venv/bin/python" -u prepare_gse249925.py
