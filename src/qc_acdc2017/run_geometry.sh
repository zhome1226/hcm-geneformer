#!/usr/bin/env bash
set -euo pipefail
root="$HOME/xiaoe-internship-data"
work="$root/_processing/all_data_cleaning_20261111/ACDC2017"
cd "$work"
exec 9>geometry.lock
flock -n 9 || exit 75
echo "$$" > geometry_pid.txt
trap 'code=$?; echo "$code" > geometry_exit_code.txt' EXIT
test "$(cat verify_exit_code.txt)" = 0
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
"$root/_processing/SCP1303_20261111/venv/bin/python" -u geometry_views.py
