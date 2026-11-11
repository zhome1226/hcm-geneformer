#!/usr/bin/env bash
set -euo pipefail
root=/home/zhome/xiaoe-internship-data
work="$root/_processing/SCP1303_review_20261111"
out="$root/processed/SCP1303_review_20261111/identity/doublets"
mkdir -p "$out"
trap 'code=$?; printf "%s\n" "$code" > "$out/exit_code.txt"' EXIT
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=2
/home/zhome/miniconda3/bin/conda install -y -p "$work/r_env" --solver libmamba \
 --override-channels -c conda-forge -c bioconda --strict-channel-priority \
 bioconductor-scdblfinder > "$out/install.log" 2>&1
/home/zhome/miniconda3/bin/conda list -p "$work/r_env" --explicit > "$work/r_environment.explicit.txt"
"$root/_processing/SCP1303_20261111/venv/bin/python" -u "$work/prepare_secondary_doublets.py" > "$out/prepare.log" 2>&1
"$work/r_env/bin/Rscript" --vanilla "$work/run_secondary_doublets.R" "$out" > "$out/run.log" 2>&1
