#!/usr/bin/env bash
set -euo pipefail
root=/home/zhome/xiaoe-internship-data
work="$root/_processing/all_data_cleaning_20261111/gse181764"
exec 9>"$work/export.lock"
flock -n 9 || { echo 'Another GSE181764 export holds the lock'; exit 73; }
trap 'code=$?; printf "%s\n" "$code" > "$work/resume_exit_code.txt"' EXIT
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=2
test "$(cat "$work/R_exit_code.txt")" = 0
"$root/_processing/SCP1303_20261111/venv/bin/python" -u "$work/finalize_gse181764.py" > "$work/export_resume.log" 2>&1
"$root/_processing/SCP1303_20261111/venv/bin/python" -u "$work/verify_gse181764.py" > "$work/verify_resume.log" 2>&1
