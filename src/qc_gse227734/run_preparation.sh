#!/usr/bin/env bash
set -euo pipefail
work=/home/zhome/xiaoe-internship-data/_processing/all_data_cleaning_20261111/GSE227734
py=/home/zhome/xiaoe-internship-data/_processing/SCP1303_20261111/venv/bin/python
exec 9>"$work/preparation.lock"
flock -n 9 || exit 73
trap 'code=$?; printf "%s\n" "$code" > "$work/preparation_exit_code.txt"' EXIT
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
cd /home/zhome/xiaoe-internship-data/validation/cox7b_GSE227734
md5sum -c MD5SUMS > "$work/original_md5_recheck.log"
"$py" -u "$work/prepare_gse227734.py" > "$work/preparation.log" 2>&1
"$py" -u "$work/verify_gse227734.py" > "$work/preparation_verification.log" 2>&1
