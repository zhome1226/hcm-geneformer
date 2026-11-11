#!/usr/bin/env bash
set -euo pipefail
root=/home/zhome/xiaoe-internship-data
work="$root/_processing/all_data_cleaning_20261111/chin_pair"
trap 'code=$?; printf "%s\n" "$code" > "$work/exit_code.txt"' EXIT
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=2
"$root/_processing/SCP1303_20261111/venv/bin/python" -u "$work/stage_chin.py" > "$work/stage.log" 2>&1
"$root/_processing/all_data_cleaning_20261111/r_env/bin/Rscript" --vanilla "$work/run_chin.R" > "$work/run.log" 2>&1
bash "$work/export_verify_chin.sh"
