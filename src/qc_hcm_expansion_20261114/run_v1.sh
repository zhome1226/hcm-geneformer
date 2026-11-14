#!/usr/bin/env bash
set -euo pipefail
cd /home/zhome/hcm-geneformer-analysis/hcm_expansion_20261114/nuclear_qc
trap 'rc=$?; printf "%s\n" "$rc" > v1_exit_code.txt' EXIT
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMBA_NUM_THREADS=2
export HF_HUB_OFFLINE=1 HF_DATASETS_OFFLINE=1 TOKENIZERS_PARALLELISM=false
PY=/home/zhome/xiaoe-internship-data/_processing/geneformer_v1_tokenization_20261112/venv/bin/python
"$PY" -u prepare_v1.py
"$PY" -u tokenize_v1.py --mode smoke
"$PY" -u verify_v1.py --mode smoke
"$PY" -u tokenize_v1.py --mode full
"$PY" -u verify_v1.py --mode full
