#!/usr/bin/env bash

set -e

## ARGUMENTS

nproc=$1 # Number of cores
exe=$2

## SYSTEM OPTIONS

export OMP_NUM_THREADS=$nproc

# Ensure stacksize unlimited (for fortran)
ulimit -d unlimited

if [[ ${USE_MPI} -eq 1 ]]; then
  LAUNCH="time mpirun -ppn 1 -np $nproc"
  ./set_io_format.sh enable
  echo "Will run with MPI"
else
  LAUNCH="time"
  echo "Will run without MPI"
fi

echo "Time: `date`"
echo "Current directory: `pwd`"
echo "Using ${nproc} process(es)"

files_to_check=(file{1..12})
if [[ ${USE_MPI} -ne 1 ]]; then
  # The intensity file does not work with MPI at the moment
  files_to_check+=(file_intensity)
fi
for name in ${files_to_check[@]}; do
  $LAUNCH ./$exe $name.inp -o $name > $name.out
done

echo "DONE"
