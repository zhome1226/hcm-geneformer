#!/usr/bin/env bash
set -euo pipefail
root=/home/zhome/xiaoe-internship-data
work="$root/_processing/all_data_cleaning_20261111/gse181764"
trap 'code=$?; printf "%s\n" "$code" > "$work/stage_exit_code.txt"' EXIT
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=2
"$root/_processing/SCP1303_20261111/venv/bin/python" -u "$work/stage_gse181764.py"
