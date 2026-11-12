#!/usr/bin/env bash
set -u
cd /home/zhome/xiaoe-internship-data/_processing/final_cleaning_audit_20261112 || exit 2
/home/zhome/xiaoe-internship-data/_processing/SCP1303_20261111/venv/bin/python -u audit_server.py > audit.log 2>&1
result=$?
printf '%s\n' "$result" > audit_exit_code.txt
exit "$result"
