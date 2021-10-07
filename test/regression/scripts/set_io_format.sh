#!/usr/bin/env bash

files=$(ls *.inp)
option=$1

if [[ $option == "enable" ]]; then
  echo "enable"
  for f in $files; do
    if ! grep -qi MPIIO $f; then
      sed -i '/eigenfunc/ a format MPIIO' $f
    fi
  done
elif [[ $option == "disable" ]]; then
  echo "disable"
  for f in $files; do
    if grep -qi MPIIO $f; then
      sed -i '/MPIIO/d' $f
    fi
  done
else
  echo "options are enable or disable"
fi
