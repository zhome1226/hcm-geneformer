#!/usr/bin/env bash
set -euo pipefail
cd /home/zhome/hcm-geneformer-analysis/hcm_expansion_20261114/nuclear_qc
trap 'rc=$?; printf "%s\n" "$rc" > qc_exit_code.txt' EXIT
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMBA_NUM_THREADS=2
export MPLBACKEND=Agg
PY=/home/zhome/xiaoe-internship-data/_processing/SCP1303_20261111/venv/bin/python
RS=/home/zhome/xiaoe-internship-data/_processing/all_data_cleaning_20261111/r_env/bin/Rscript
"$PY" -u stage.py
"$PY" -u calibrate_mixed.py
"$RS" identity.R
"$PY" -u finalize.py
"$PY" -u verify.py
