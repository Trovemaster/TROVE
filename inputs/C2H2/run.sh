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

# Optional (RUN_P12=1 ./run.sh): converged Pmax=12 basis with K <= 5 on the CCSD(T)/CBS PES.
# The J=0 step writes ~4 GB of matelem*.chk (matelem save split); use OpenMP and an optimised BLAS/LAPACK.
if [ "${RUN_P12:-0}" = "1" ]; then
  echo "=== 4. Running converged Pmax=12 (K <= 5) CCSD(T)/CBS PES examples: J=0, 1, 5 ==="
  cp -f potential_CBS.chk potential.chk
  "$TROVE_BIN" < C2H2_sing3D_D2h_CBS_P12_J0_kmax5.inp > C2H2_sing3D_D2h_CBS_P12_J0_kmax5_run.out
  "$TROVE_BIN" < C2H2_sing3D_D2h_CBS_P12_J1.inp > C2H2_sing3D_D2h_CBS_P12_J1_run.out
  "$TROVE_BIN" < C2H2_sing3D_D2h_CBS_P12_J5.inp > C2H2_sing3D_D2h_CBS_P12_J5_run.out

  echo "=== 5. Running converged Pmax=12 (K <= 5) Empirically Refined PES examples: J=0, 1, 2, 5 ==="
  cp -f potential_Refined_P12.chk potential.chk
  "$TROVE_BIN" < C2H2_sing3D_D2h_Refined_P12_J0_kmax5.inp > C2H2_sing3D_D2h_Refined_P12_J0_kmax5_run.out
  "$TROVE_BIN" < C2H2_sing3D_D2h_Refined_P12_J1.inp > C2H2_sing3D_D2h_Refined_P12_J1_run.out
  "$TROVE_BIN" < C2H2_sing3D_D2h_Refined_P12_J2.inp > C2H2_sing3D_D2h_Refined_P12_J2_run.out
  "$TROVE_BIN" < C2H2_sing3D_D2h_Refined_P12_J5.inp > C2H2_sing3D_D2h_Refined_P12_J5_run.out
fi

# Optional (RUN_P12_DMS=1 ./run.sh): 6th-order CCSD(T)/CBS DMS (external_CBS_6th.chk) intensities
# on the refined Pmax=12 (K <= 5) basis for J = 0..5. All energy runs (EIGENFUNC save) come first;
# the DMS J=0 step writes its primitive basis to prim_bset_dms.chk (chk_file), so prim_bset.chk is kept.
if [ "${RUN_P12_DMS:-0}" = "1" ]; then
  echo "=== 6. Refined Pmax=12 (K <= 5) J=0..5 + 6th-order CBS DMS contraction and J=0..5 intensities ==="
  cp -f potential_Refined_P12.chk potential.chk
  cp -f external.chk external_benchmark_backup.chk
  cp -f external_CBS_6th.chk external.chk
  "$TROVE_BIN" < C2H2_sing3D_D2h_Refined_P12_J0_kmax5.inp > C2H2_sing3D_D2h_Refined_P12_J0_kmax5_run.out
  for J in 1 2 3 4 5; do
    "$TROVE_BIN" < C2H2_sing3D_D2h_Refined_P12_J${J}.inp > C2H2_sing3D_D2h_Refined_P12_J${J}_run.out
  done
  "$TROVE_BIN" < C2H2_sing3D_D2h_Refined_P12_DMS_J0_kmax5.inp > C2H2_sing3D_D2h_Refined_P12_DMS_J0_kmax5_run.out
  "$TROVE_BIN" < C2H2_sing3D_D2h_Refined_P12_intens.inp > C2H2_sing3D_D2h_Refined_P12_intens_run.out
  mv -f external_benchmark_backup.chk external.chk
fi

# Optional (RUN_P14=1 ./run.sh): converged Pmax=14 (lmax=7, kmax=7) Empirically Refined PES examples.
if [ "${RUN_P14:-0}" = "1" ]; then
  echo "=== 7. Running converged Pmax=14 (lmax=7, kmax=7) Empirically Refined PES examples: J=0, 1, 2, 5, 10, 15, 20 ==="
  cp -f potential_Refined_P14.chk potential.chk
  "$TROVE_BIN" < C2H2_sing3D_D2h_Refined_P14_J0_kmax7.inp > C2H2_sing3D_D2h_Refined_P14_J0_kmax7_run.out
  for J in 1 2 5 10 15 20; do
    "$TROVE_BIN" < C2H2_sing3D_D2h_Refined_P14_J${J}.inp > C2H2_sing3D_D2h_Refined_P14_J${J}_run.out
  done
fi

# Optional (RUN_P14_DMS=1 ./run.sh): 6th-order CCSD(T)/CBS DMS (external_CBS_6th.chk) intensities
# on the refined Pmax=14 (lmax=7, kmax=7) basis for J = 0..5.
if [ "${RUN_P14_DMS:-0}" = "1" ]; then
  echo "=== 8. Refined Pmax=14 (lmax=7, kmax=7) J=0..5 + 6th-order CBS DMS contraction and J=0..5 intensities ==="
  cp -f potential_Refined_P14.chk potential.chk
  cp -f external.chk external_benchmark_backup.chk
  cp -f external_CBS_6th.chk external.chk
  "$TROVE_BIN" < C2H2_sing3D_D2h_Refined_P14_J0_kmax7.inp > C2H2_sing3D_D2h_Refined_P14_J0_kmax7_run.out
  for J in 1 2 3 4 5; do
    "$TROVE_BIN" < C2H2_sing3D_D2h_Refined_P14_J${J}.inp > C2H2_sing3D_D2h_Refined_P14_J${J}_run.out
  done
  "$TROVE_BIN" < C2H2_sing3D_D2h_Refined_P14_DMS_J0_kmax7.inp > C2H2_sing3D_D2h_Refined_P14_DMS_J0_kmax7_run.out
  "$TROVE_BIN" < C2H2_sing3D_D2h_Refined_P14_intens.inp > C2H2_sing3D_D2h_Refined_P14_intens_run.out
  mv -f external_benchmark_backup.chk external.chk
fi

mv -f potential_test_backup.chk potential.chk
