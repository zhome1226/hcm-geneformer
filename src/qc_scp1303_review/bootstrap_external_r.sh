#!/usr/bin/env bash
set -euo pipefail
root=/home/zhome/xiaoe-internship-data
work="$root/_processing/all_data_cleaning_20261111"
source_env="$root/_processing/SCP1303_review_20261111/r_env"
mkdir -p "$work"
trap 'code=$?; printf "%s\n" "$code" > "$work/r_bootstrap_exit_code.txt"' EXIT
if [[ ! -x "$work/r_env/bin/Rscript" ]]; then
 /home/zhome/miniconda3/bin/conda create -y -p "$work/r_env" --clone "$source_env"
fi
/home/zhome/miniconda3/bin/conda install -y -p "$work/r_env" --solver libmamba \
 --override-channels -c conda-forge -c bioconda --strict-channel-priority \
 'r-base=4.5.*' bioconductor-dropletutils r-soupx r-seuratobject
/home/zhome/miniconda3/bin/conda list -p "$work/r_env" --explicit > "$work/r_environment.explicit.txt"
"$work/r_env/bin/Rscript" --vanilla -e 'stopifnot(requireNamespace("DropletUtils"), requireNamespace("SoupX"), requireNamespace("SeuratObject"), requireNamespace("SingleR"), requireNamespace("scDblFinder")); sessionInfo()' > "$work/R_environment_smoke.txt" 2>&1
