#!/bin/bash
# HCN/HNC (Jacobi KEO + RELAXED two-well basis + SPAR PES + Jacobi-frame DMS) TROVE example.
#
# Usage: ./run.sh [path/to/j-trove.x]
#
# 1) J=0 energies; saves basis, contraction, matrix elements, the DMS field (external.chk)
#    and the J=0 eigenvectors.
# 2) J=0 vibrational transition moments (INTENSITY tm) from the step-1 checkpoints.
# 3) J=0 (kmax=2, MATELEM SAVE), J=1 and J=2 rovibrational energies
#    (B0, D0, nu2 l-doubling of HCN and HNC). NB: the J=0 kmax=2 step writes a
#    ~29 GB contr_matelem.chk and J=2 needs ~27 GB of RAM.
#
# All runs need potential.chk (SPAR form of the PES, POTENTIAL READ WITH_MODES)
# in the working directory. external.chk is provided for reference; set
# "EXTERNAL READ" in the step-1 input to reuse it instead of recomputing it.

TROVE=${1:-../../j-trove.x}

ulimit -s unlimited
export OMP_STACKSIZE=1G
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-12}

set -e

$TROVE < HCN_HNC_J0_step1.inp    > HCN_HNC_J0_step1.out
$TROVE < HCN_HNC_J0_step2_tm.inp > HCN_HNC_J0_step2_tm.out

mkdir -p rot && cp potential.chk HCN_HNC_J0.inp HCN_HNC_J1.inp HCN_HNC_J2.inp rot/
cd rot
for j in 0 1 2; do
  $TROVE < HCN_HNC_J$j.inp > HCN_HNC_J$j.out
done
