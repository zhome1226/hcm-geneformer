#!/usr/bin/env bash
set -euo pipefail
root="$HOME/xiaoe-internship-data"
work="$root/_processing/all_data_cleaning_20261111/gse138262"
cd "$work"
exec 9>finalize.lock
flock -n 9 || exit 75
echo "$$" > finalize_pid.txt
trap 'code=$?; echo "$code" > finalize_exit_code.txt' EXIT
test "$(cat stage_exit_code.txt)" = 0
test "$(cat identity_exit_code.txt)" = 0
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=2
py="$root/_processing/SCP1303_20261111/venv/bin/python"
"$py" -u calibrate_mixed.py > mixed_reference.log 2>&1
"$py" -u finalize_gse138262.py
