#!/usr/bin/env bash

set -e

## ARGUMENTS

nproc=$1 # Number of cores
exe=$2

## SYSTEM OPTIONS

# Ensure stacksize unlimited (for fortran)
ulimit -d unlimited

if [[ ${USE_MPI} == 1 ]]; then
  echo "MPI enabled"
  LAUNCH="time mpirun -np $nproc"
  ./set_io_format.sh enable
  export OMP_NUM_THREADS=1
else
  echo "MPI disabled"
  LAUNCH="time"
  ./set_io_format.sh disable
  export OMP_NUM_THREADS=$nproc
fi

echo "Time: `date`"
echo "Current directory: `pwd`"
echo "Using ${nproc} process(es)"

files_to_check=(file{1..12})
if [[ ${USE_MPI} == 0 ]]; then
  # The intensity file does not work with MPI at the moment
  files_to_check+=(file_intensity)
fi
for name in ${files_to_check[@]}; do
  $LAUNCH ./$exe $name.inp -o $name > $name.out
done

echo "DONE"
