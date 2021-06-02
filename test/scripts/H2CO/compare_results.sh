#!/usr/bin/env bash

set -e # fail on errors

folder1=$1
folder2=$2

# These are descriptions of quantum states; the energies must be compared
# Quantum states can be equivalent without being equal so cannot be directly compared
quantum_files="eigen_descr0_1.chk eigen_descr0_2.chk eigen_descr0_3.chk eigen_descr0_4.chk j0contr_descr.chk j0eigen_descr0_1.chk j0eigen_descr0_2.chk j0eigen_descr0_3.chk j0eigen_descr0_4.chk j0eigen_descr1_1.chk j0eigen_descr1_2.chk j0eigen_descr1_3.chk j0eigen_descr1_4.chk j0eigen_descr2_1.chk j0eigen_descr2_2.chk j0eigen_descr2_3.chk j0eigen_descr2_4.chk j0eigen_descr3_1.chk j0eigen_descr3_2.chk j0eigen_descr3_3.chk j0eigen_descr3_4.chk j0eigen_descr4_1.chk j0eigen_descr4_2.chk j0eigen_descr4_3.chk j0eigen_descr4_4.chk j0eigen_descr10_1.chk j0eigen_descr10_2.chk j0eigen_descr10_3.chk j0eigen_descr10_4.chk j0eigen_descr11_1.chk j0eigen_descr11_2.chk j0eigen_descr11_3.chk j0eigen_descr11_4.chk j0eigen_descr12_1.chk j0eigen_descr12_2.chk j0eigen_descr12_3.chk j0eigen_descr12_4.chk" 

python compare_results.py --kind quantum --folder1 "$folder1" --folder2 "$folder2" $quantum_files

python compare_results.py --kind intensity --precision 1e-6 --folder1 "$folder1" --folder2 "$folder2" file_intensity.out

python compare_results.py --kind column --column 3 --precision 5e-3 --folder1 "$folder1" --folder2 "$folder2" external.chk

python compare_results.py --kind column --column 2 --precision 1e-8 --folder1 "$folder1" --folder2 "$folder2" potential.chk

python compare_results.py --kind column --column 4 --precision 1e-8 --folder1 "$folder1" --folder2 "$folder2" kinetic.chk
