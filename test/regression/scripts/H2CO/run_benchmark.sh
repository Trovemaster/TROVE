#!/usr/bin/env bash

set -e

## ARGUMENTS

nproc=$1 # Number of cores
exe=$2

## SYSTEM OPTIONS

export OMP_NUM_THREADS=$nproc

# Ensure stacksize unlimited (for fortran)
ulimit -d unlimited

if [ -n "${USE_MPI}" ]; then
  echo "MPI enabled"
  LAUNCH="mpirun -np 4 --mca opal_warn_on_missing_libcuda 0"
else
  echo "MPI disabled"
  LAUNCH="time"
fi

echo "Time: `date`"
echo "Current directory: `pwd`"

for name in file{1..12} file_intensity; do
  $LAUNCH ./$exe $name.inp > $name.out
done

echo "DONE"
