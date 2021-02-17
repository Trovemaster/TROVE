#!/usr/bin/env bash

set -e

## ARGUMENTS

nproc=$1 # Number of cores

## DEFAULT OPTIONS

exe=j-trove.x

## SYSTEM OPTIONS

export OMP_NUM_THREADS=$nproc

# Ensure stacksize unlimited (for fortran)
ulimit -d unlimited

LAUNCH="time"

echo "Time: `date`"
echo "Current directory: `pwd`"

for name in file{1..12} file_intensity; do
  $LAUNCH ./$exe < $name.inp > $name.out
done

echo "DONE"
