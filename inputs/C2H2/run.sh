#!/bin/bash
set -euo pipefail

TROVE_BIN="${TROVE_BIN:-../../j-trove.x}"

echo "=== 1. Running C2H2 J=0 (K <= 2) on test PES ==="
cp -f potential.chk potential_active.chk
ln -sf potential_active.chk potential.chk 2>/dev/null || true
"$TROVE_BIN" < C2H2_sing3D_D2h_J0_kmax2.inp > C2H2_sing3D_D2h_J0_kmax2_run.out

echo "=== 2. Running C2H2 J=1 (K <= 1) rovibrational step ==="
"$TROVE_BIN" < C2H2_sing3D_D2h_J1.inp > C2H2_sing3D_D2h_J1_run.out

echo "=== 3. Running C2H2 J=0 (K <= 2) on new CCSD(T)/CBS PES ==="
cp -f potential_CBS.chk potential.chk
"$TROVE_BIN" < C2H2_sing3D_D2h_CBS_PES_J0.inp > C2H2_sing3D_D2h_CBS_PES_J0_run.out
