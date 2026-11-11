#!/usr/bin/env bash
set -euo pipefail
root=/home/zhome/xiaoe-internship-data
work="$root/_processing/all_data_cleaning_20261111/gse181764"
trap 'code=$?; printf "%s\n" "$code" > "$work/R_exit_code.txt"' EXIT
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=2
test "$(cat "$work/stage_exit_code.txt")" = 0
"$root/_processing/SCP1303_20261111/venv/bin/python" -u "$work/calibrate_mixed.py" > "$work/mixed_reference.log" 2>&1
"$root/_processing/all_data_cleaning_20261111/r_env/bin/Rscript" --vanilla "$work/run_gse181764.R"
