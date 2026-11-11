#!/usr/bin/env bash
set -euo pipefail
work="$HOME/xiaoe-internship-data/_processing/all_data_cleaning_20261111/gse138262"
cd "$work"
exec 9>download.lock
flock -n 9 || exit 75
echo "$$" > download_pid.txt
trap 'code=$?; echo "$code" > download_exit_code.txt' EXIT
"$HOME/xiaoe-internship-data/_processing/SCP1303_20261111/venv/bin/python" -u download_corrected.py
