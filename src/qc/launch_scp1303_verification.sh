#!/usr/bin/env bash
set -euo pipefail
work=/home/zhome/xiaoe-internship-data/_processing/SCP1303_20261111
out=/home/zhome/xiaoe-internship-data/processed/SCP1303_qc_20261111
test -f "$out/SUCCESS.json"
if [[ -f "$out/verification.pid" ]] && kill -0 "$(cat "$out/verification.pid")" 2>/dev/null; then
  echo 'Verification is already running'
  exit 0
fi
cd "$work"
nohup env MPLBACKEND=Agg OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 NUMBA_NUM_THREADS=4 \
  "$work/venv/bin/python" -u "$work/verify_scp1303_qc.py" \
  --input "$work/author_input.h5ad" --output "$out" \
  > "$out/verification.log" 2>&1 < /dev/null &
echo "$!" > "$out/verification.pid"
echo "Started independent verification: PID $(cat "$out/verification.pid")"
