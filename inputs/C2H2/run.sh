#!/bin/bash
set -euo pipefail

TROVE_BIN="${TROVE_BIN:-../../j-trove.x}"

echo "=== 1. Running C2H2 J=0, 1, 2 (K <= 2, Pmax=12) on test PES ==="
cp -f potential.chk potential_test_backup.chk
"$TROVE_BIN" < C2H2_sing3D_D2h_J0_kmax2.inp > C2H2_sing3D_D2h_J0_kmax2_run.out
"$TROVE_BIN" < C2H2_sing3D_D2h_J1.inp > C2H2_sing3D_D2h_J1_run.out
"$TROVE_BIN" < C2H2_sing3D_D2h_J2.inp > C2H2_sing3D_D2h_J2_run.out

echo "=== 2. Running C2H2 J=0, 1, 2 (K <= 2, Pmax=12) on 4,331-pt CCSD(T)/CBS PES ==="
cp -f potential_CBS.chk potential.chk
"$TROVE_BIN" < C2H2_sing3D_D2h_CBS_PES_J0.inp > C2H2_sing3D_D2h_CBS_PES_J0_run.out
"$TROVE_BIN" < C2H2_sing3D_D2h_CBS_PES_J1.inp > C2H2_sing3D_D2h_CBS_PES_J1_run.out
"$TROVE_BIN" < C2H2_sing3D_D2h_CBS_PES_J2.inp > C2H2_sing3D_D2h_CBS_PES_J2_run.out

echo "=== 3. Running 6D DMS contraction and J=0..2 absorption intensities ==="
"$TROVE_BIN" < C2H2_sing3D_D2h_CBS_DMS_J0.inp > C2H2_sing3D_D2h_CBS_DMS_J0_run.out
"$TROVE_BIN" < C2H2_sing3D_D2h_CBS_intens.inp > C2H2_sing3D_D2h_CBS_intens_run.out

mv -f potential_test_backup.chk potential.chk
