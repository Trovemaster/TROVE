#!/usr/bin/env bash

set -e

## ARGUMENTS

nproc=$1 # Number of cores
exe=$2

## SYSTEM OPTIONS

export OMP_NUM_THREADS=$nproc

# Ensure stacksize unlimited (for fortran)
ulimit -d unlimited

LAUNCH="time"

echo "Time: `date`"
echo "Current directory: `pwd`"

$LAUNCH ./$exe < H2CS_J0.inp > H2CS_J0.out

echo "DONE"
