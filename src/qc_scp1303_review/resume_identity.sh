#!/usr/bin/env bash
set -euo pipefail
root=/home/zhome/xiaoe-internship-data
work="$root/_processing/SCP1303_review_20261111"
out="$root/processed/SCP1303_review_20261111/identity"
trap 'code=$?; printf "%s\n" "$code" > "$out/exit_code.txt"' EXIT
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=4
printf 'installing_wilcox_dependency\n' > "$out/stage.txt"
/home/zhome/miniconda3/bin/conda install -y -p "$work/r_env" --solver libmamba \
  --override-channels -c conda-forge -c bioconda --strict-channel-priority \
  bioconductor-scrapper > "$out/install_scrapper.log" 2>&1
/home/zhome/miniconda3/bin/conda list -p "$work/r_env" --explicit > "$work/r_environment.explicit.txt"
printf 'running_singler\n' > "$out/stage.txt"
"$work/r_env/bin/Rscript" --vanilla "$work/run_identity_singler.R" "$out" > "$out/singler.log" 2>&1
printf 'annotation_complete\n' > "$out/stage.txt"
