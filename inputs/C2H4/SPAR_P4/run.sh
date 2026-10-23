#!/usr/bin/env bash
# Runs this example in the folder of this script (copy the folder to a working directory first).
# TROVE executable: $TROVE_EXE if set, otherwise j-trove.x of the TROVE tree this folder is in,
# otherwise j-trove.x from the PATH.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
EXE=${TROVE_EXE:-}
if [ -z "$EXE" ]; then
  if [ -x "$HERE/../../../j-trove.x" ]; then EXE="$HERE/../../../j-trove.x"; else EXE=j-trove.x; fi
fi
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
for inp in C2H4_SPAR_P4_J0 C2H4_SPAR_P4_J0_convert C2H4_SPAR_P4_J0_j0 C2H4_SPAR_P4_J1_j0; do
  echo "running $inp"
  "$EXE" < $inp.inp > $inp.out
done
