#!/usr/bin/env bash
set -euo pipefail
root=/home/zhome/xiaoe-internship-data
work="$root/_processing/all_data_cleaning_20261111/chin_pair"
trap 'code=$?; printf "%s\n" "$code" > "$work/export_exit_code.txt"' EXIT
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=2
"$root/_processing/SCP1303_20261111/venv/bin/python" -u "$work/finalize_chin.py" > "$work/export.log" 2>&1
"$root/_processing/SCP1303_20261111/venv/bin/python" -u "$work/label_diagnostic_metadata.py" > "$work/diagnostic_label.log" 2>&1
"$root/_processing/SCP1303_20261111/venv/bin/python" -u "$work/verify_chin.py" > "$work/verify.log" 2>&1
