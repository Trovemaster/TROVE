#!/usr/bin/env bash

folder1=$1
folder2=$2

process_quantum_energies () {
  filename=$1

  fingerprints=$(awk '/^Start Fingerprints/,/^End Fingerprints/ {print $0}' $filename)
  printf "$fingerprints"
  echo -e "\n"

  quantum_numbers=$(awk '/^Start Quantum/,/^End Quantum/ {print $0}' $filename)
  printf "$quantum_numbers" | awk '!/number|Class/ {printf("%0.6e\n", $5);}'
}

process_simple_files () {
  filename=$1
  cat $filename | awk '{printf("%0.10e\n", $2);}'
}

simple_files="external.chk kinetic.chk potential.chk"

complex_files="contr_descr.chk eigen_descr0_1.chk eigen_descr0_2.chk eigen_descr0_3.chk eigen_descr0_4.chk j0contr_descr.chk j0eigen_descr0_1.chk j0eigen_descr0_2.chk j0eigen_descr0_3.chk j0eigen_descr0_4.chk j0eigen_descr1_1.chk j0eigen_descr1_2.chk j0eigen_descr1_3.chk j0eigen_descr1_4.chk j0eigen_descr2_1.chk j0eigen_descr2_2.chk j0eigen_descr2_3.chk j0eigen_descr2_4.chk j0eigen_descr3_1.chk j0eigen_descr3_2.chk j0eigen_descr3_3.chk j0eigen_descr3_4.chk j0eigen_descr4_1.chk j0eigen_descr4_2.chk j0eigen_descr4_3.chk j0eigen_descr4_4.chk j0eigen_descr10_1.chk j0eigen_descr10_2.chk j0eigen_descr10_3.chk j0eigen_descr10_4.chk j0eigen_descr11_1.chk j0eigen_descr11_2.chk j0eigen_descr11_3.chk j0eigen_descr11_4.chk j0eigen_descr12_1.chk j0eigen_descr12_2.chk j0eigen_descr12_3.chk j0eigen_descr12_4.chk" 

let return_sum=0

for f in $simple_files $complex_files; do
  if grep -q 'Start Quantum' $folder1/$f; then
    diff -u <(process_quantum_energies $folder1/$f) <(process_quantum_energies $folder2/$f)
  else
    diff -u <(process_simple_files $folder1/$f) <(process_simple_files $folder2/$f)
  fi
  ret=$?
  if [[ $ret -ne 0 ]]; then
    echo "Differences in $f"
  fi
  let return_sum+=$ret
done

if [ $return_sum -eq 0 ]; then
  exit 0
else
  exit 1
fi

