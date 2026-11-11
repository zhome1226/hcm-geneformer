#!/usr/bin/env bash
set -euo pipefail
root=/home/zhome/xiaoe-internship-data
work="$root/_processing/SCP1303_review_20261111"
out="$root/processed/SCP1303_qc_reviewed_final_20261111"
mkdir -p "$out"
trap 'code=$?; printf "%s\n" "$code" > "$out/exit_code.txt"' EXIT
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=2
"$root/_processing/SCP1303_20261111/venv/bin/python" -u "$work/export_reviewed_scp1303.py" --output-name SCP1303_qc_reviewed_final_20261111 --use-group-evidence > "$out/export.log" 2>&1
"$root/_processing/SCP1303_20261111/venv/bin/python" -u "$work/verify_reviewed_scp1303.py" --output-name SCP1303_qc_reviewed_final_20261111 > "$out/verification.log" 2>&1
