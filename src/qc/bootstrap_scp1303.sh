#!/usr/bin/env bash
set -euo pipefail
ROOT="$HOME/xiaoe-internship-data"
WORK="$ROOT/_processing/SCP1303_20261111"
mkdir -p "$WORK" "$WORK/skill"
export MPLBACKEND=Agg
export OPENBLAS_NUM_THREADS=4
export OMP_NUM_THREADS=4
export NUMBA_NUM_THREADS=4
if [ ! -x "$WORK/venv/bin/python" ]; then
  "$HOME/miniconda3/bin/python" -m venv "$WORK/venv"
fi
"$WORK/venv/bin/python" -m pip install --disable-pip-version-check --index-url https://pypi.tuna.tsinghua.edu.cn/simple 'scanpy>=1.11,<1.13' h5py openpyxl scikit-image > "$WORK/install.log" 2>&1
"$WORK/venv/bin/python" -m pip freeze > "$WORK/environment.freeze.txt"
"$WORK/venv/bin/python" "$WORK/prepare_input.py"
"$WORK/venv/bin/python" "$WORK/skill/inspect_data.py" "$WORK/author_input.h5ad" > "$WORK/skill_inspection.txt" 2>&1
"$WORK/venv/bin/python" "$WORK/inspect_metadata.py" > "$WORK/metadata_inspection.log" 2>&1
printf 'Bootstrap and inspection complete\n'
