#!/usr/bin/env bash

set -e

## ARGUMENTS

nproc=$1 # Number of cores
exe=$2

## SYSTEM OPTIONS

export OMP_NUM_THREADS=$nproc

# Ensure stacksize unlimited (for fortran)
ulimit -d unlimited

if [ ${USE_MPI} = 1 ]; then
  LAUNCH="time mpirun -ppn 1 -np $nproc"
  echo "Will run with MPI"
else
  LAUNCH="time"
  echo "Will run without MPI"
fi

echo "Time: `date`"
echo "Current directory: `pwd`"
echo "Using ${nproc} process(es)"

for name in file{1..12} file_intensity; do
  $LAUNCH ./$exe < $name.inp > $name.out
done

echo "DONE"
