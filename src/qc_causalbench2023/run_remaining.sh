#!/usr/bin/env bash
# Continue the already-running preparation without launching it a second time.
set -euo pipefail
root="$HOME/xiaoe-internship-data"
work="$root/_processing/all_data_cleaning_20261111/CausalBench2023"
cd "$work"
exec 9>remaining.lock
flock -n 9 || exit 75
echo "$$" > remaining_pid.txt
trap 'code=$?; echo "$code" > remaining_exit_code.txt' EXIT
prepare_pid="$(cat prepare_pid.txt)"
while test -d "/proc/$prepare_pid"; do sleep 10; done
test "$(cat prepare_exit_code.txt)" = 0
test ! -e verify_pid.txt
"$root/_processing/SCP1303_20261111/venv/bin/python" -m py_compile verify.py collect_completion.py
bash run_verify.sh > verify.log 2>&1
"$root/_processing/SCP1303_20261111/venv/bin/python" collect_completion.py > collect_completion.log 2>&1
