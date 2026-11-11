#!/usr/bin/env bash
set -euo pipefail
work=/home/zhome/xiaoe-internship-data/_processing/all_data_cleaning_20261111/liu2023
trap 'code=$?; printf "%s\n" "$code" > "$work/preflight_exit_code.txt"' EXIT
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=2
/home/zhome/xiaoe-internship-data/_processing/SCP1303_20261111/venv/bin/python -u "$work/preflight_liu.py"
