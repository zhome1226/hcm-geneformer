#!/usr/bin/env bash
set -euo pipefail
root=/home/zhome/xiaoe-internship-data
work="$root/_processing/all_data_cleaning_20261111/liu2023"
trap 'code=$?; printf "%s\n" "$code" > "$work/cleaning_exit_code.txt"' EXIT
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=2
test "$(cat "$work/preflight_exit_code.txt")" = 0
"$root/_processing/all_data_cleaning_20261111/r_env/bin/Rscript" --vanilla "$work/run_liu.R" > "$work/R_cleaning.log" 2>&1
"$root/_processing/SCP1303_20261111/venv/bin/python" -u "$work/review_mixed_liu.py" > "$work/mixed_review.log" 2>&1
"$root/_processing/SCP1303_20261111/venv/bin/python" -u "$work/finalize_liu.py" > "$work/export.log" 2>&1
"$root/_processing/SCP1303_20261111/venv/bin/python" -u "$work/verify_liu.py" > "$work/verify.log" 2>&1
"$root/_processing/SCP1303_20261111/venv/bin/python" -u "$work/summarize_retention_liu.py" > "$work/attrition.log" 2>&1
