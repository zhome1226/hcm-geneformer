#!/usr/bin/env bash
set -euo pipefail
root=/home/zhome/xiaoe-internship-data
work="$root/_processing/all_data_cleaning_20261111/liu2023"
trap 'code=$?; printf "%s\n" "$code" > "$work/final_v2_exit_code.txt"' EXIT
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=2
test -f "$work/R_SUCCESS.txt"
test -f "$work/MIXED_PROGRAM_REVIEW.json"
"$root/_processing/SCP1303_20261111/venv/bin/python" -u "$work/finalize_liu.py" > "$work/export_v2.log" 2>&1
"$root/_processing/SCP1303_20261111/venv/bin/python" -u "$work/verify_liu.py" > "$work/verify_v2.log" 2>&1
"$root/_processing/SCP1303_20261111/venv/bin/python" -u "$work/summarize_retention_liu.py" > "$work/attrition.log" 2>&1
