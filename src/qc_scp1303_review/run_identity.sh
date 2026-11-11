#!/usr/bin/env bash
set -euo pipefail
root=/home/zhome/xiaoe-internship-data
work="$root/_processing/SCP1303_review_20261111"
out="$root/processed/SCP1303_review_20261111/identity"
mkdir -p "$out"
trap 'code=$?; printf "%s\n" "$code" > "$out/exit_code.txt"' EXIT
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=4
printf 'preparing_reference\n' > "$out/stage.txt"
"$root/_processing/SCP1303_20261111/venv/bin/python" -u "$work/prepare_identity_reference.py" > "$out/prepare.log" 2>&1
printf 'running_singler\n' > "$out/stage.txt"
"$work/r_env/bin/Rscript" --vanilla "$work/run_identity_singler.R" "$out" > "$out/singler.log" 2>&1
printf 'annotation_complete\n' > "$out/stage.txt"
