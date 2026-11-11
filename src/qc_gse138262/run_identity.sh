#!/usr/bin/env bash
set -euo pipefail
root="$HOME/xiaoe-internship-data"
work="$root/_processing/all_data_cleaning_20261111/gse138262"
cd "$work"
exec 9>identity.lock
flock -n 9 || exit 75
echo "$$" > identity_pid.txt
trap 'code=$?; echo "$code" > identity_exit_code.txt' EXIT
test "$(cat stage_exit_code.txt)" = 0
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=2
"$root/_processing/all_data_cleaning_20261111/r_env/bin/Rscript" --vanilla run_identity.R
