#!/usr/bin/env bash
set -euo pipefail
root="$HOME/xiaoe-internship-data"
work="$root/_processing/all_data_cleaning_20261111/GSE203274"
cd "$work"
exec 9>R.lock
flock -n 9 || exit 75
echo "$$" > R_pid.txt
trap 'code=$?; echo "$code" > R_exit_code.txt' EXIT
test "$(cat stage_exit_code.txt)" = 0
test "$(cat repair_metadata_exit_code.txt)" = 0
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 NUMBA_NUM_THREADS=4
py="$root/_processing/SCP1303_20261111/venv/bin/python"
rscript="$root/_processing/all_data_cleaning_20261111/r_env/bin/Rscript"
"$py" -u calibrate_mixed.py > mixed_reference.log 2>&1
"$py" -u prepare_residual_review.py > residual_preparation.log 2>&1
"$rscript" --vanilla run_identity.R > identity.log 2>&1
"$rscript" --vanilla run_residual.R > residual_R.log 2>&1
