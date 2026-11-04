#!/usr/bin/env bash
# Runs the P=2 test of this example (step 1: J=0 energies, basis set and dipole matrix elements;
# step 2: J=0 vibrational transition moments and ExoMol .states/.trans files) in the folder of this
# script (copy the folder to a working directory first).
# TROVE executable: $TROVE_EXE if set, otherwise j-trove.x of the TROVE tree this folder is in,
# otherwise j-trove.x from the PATH.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
EXE=${TROVE_EXE:-}
if [ -z "$EXE" ]; then
  for exe in "$HERE/../../../j-trove.x" "$HERE/../../../../j-trove.x"; do
    if [ -x "$exe" ]; then
      EXE="$exe"
      break
    fi
  done
  EXE=${EXE:-j-trove.x}
fi
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
for inp in C2H6_CBS-TQ_6D_P2_step1_J0 C2H6_CBS-TQ_6D_P2_step2_TM; do
  echo "running $inp"
  "$EXE" < $inp.inp > $inp.out
done
