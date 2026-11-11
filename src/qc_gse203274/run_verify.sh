#!/usr/bin/env bash
set -euo pipefail
root="$HOME/xiaoe-internship-data"
work="$root/_processing/all_data_cleaning_20261111/GSE203274"
cd "$work"
exec 9>verify.lock
flock -n 9 || exit 75
echo "$$" > verify_pid.txt
trap 'code=$?; echo "$code" > verify_exit_code.txt' EXIT
test "$(cat finalize_exit_code.txt)" = 0
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 NUMBA_NUM_THREADS=4
"$root/_processing/SCP1303_20261111/venv/bin/python" -u verify_gse203274.py
