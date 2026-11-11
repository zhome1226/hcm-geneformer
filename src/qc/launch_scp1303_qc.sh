#!/usr/bin/env bash
set -euo pipefail
work=/home/zhome/xiaoe-internship-data/_processing/SCP1303_20261111
out=/home/zhome/xiaoe-internship-data/processed/SCP1303_qc_20261111
mkdir -p "$out"
if [[ -f "$out/job.pid" ]] && kill -0 "$(cat "$out/job.pid")" 2>/dev/null; then
  echo "QC job is already running: $(cat "$out/job.pid")"
  exit 0
fi
if [[ -f "$out/SUCCESS.json" ]]; then
  echo 'QC export already complete; refusing to restart'
  exit 0
fi
cd "$work"
nohup env MPLBACKEND=Agg OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 NUMBA_NUM_THREADS=4 \
  "$work/venv/bin/python" -u "$work/run_scp1303_qc.py" \
  --input "$work/author_input.h5ad" --output "$out" \
  --seed 20261111 --light-size 20000 > "$out/run.log" 2>&1 < /dev/null &
echo "$!" > "$out/job.pid"
echo "Started QC: PID $(cat "$out/job.pid"), log $out/run.log"
