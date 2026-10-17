#!/usr/bin/env bash

set -e # fail on errors

folder1=$1
folder2=$2

quantum_files="contr_descr.chk eigen_descr0_1.chk eigen_descr0_2.chk eigen_descr0_3.chk eigen_descr0_4.chk"

if command -v pipenv >/dev/null 2>&1 && pipenv --venv >/dev/null 2>&1; then
  PYTHON_CMD="pipenv run python"
else
  PYTHON_CMD="python3"
fi

$PYTHON_CMD compare_results.py --kind quantum --folder1 "$folder1" --folder2 "$folder2" $quantum_files

$PYTHON_CMD compare_results.py --kind column --column 3 --precision 1e-2 --folder1 "$folder1" --folder2 "$folder2" external.chk

$PYTHON_CMD compare_results.py --kind column --column 2 --precision 1e-8 --folder1 "$folder1" --folder2 "$folder2" potential.chk
