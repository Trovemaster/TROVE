#!/usr/bin/env bash

folder1=$1
folder2=$2

precision=5 # round to this many significant figures

extract_column () {
  # Extract column for comparison
  filename=$1
  column_no=$2

  # print only single column to certain precision and round 
  cat $filename | awk '{
  if(sqrt(($'$column_no')^2) > 1e-10)
    if($'$column_no' > 10)
      printf("%0.4e\n", $'$column_no');
    else
      printf("%0.2f\n", $'$column_no');
  else
    print 0
  }' | awk '{
  if(sqrt(($0)^2) > 1e-10)
      print $0
  else
    print 0
  }'

}

compare_columns () {
  col1=$1
  col2=$2

  merged=$(paste <(echo $1) <(echo $2))
  printf "$merged" | awk '{printf("%s\t%5.1f\n", $0, sqrt(($2-$1)^2))}'
}

# These are files where a single column should be compared
column_files="external.chk kinetic.chk potential.chk"

# These are descriptions of quantum states; the energies must be compared
# Quantum states can be equivalent without being equal so cannot be directly compared
quantum_files="eigen_descr0_1.chk eigen_descr0_2.chk eigen_descr0_3.chk eigen_descr0_4.chk j0contr_descr.chk j0eigen_descr0_1.chk j0eigen_descr0_2.chk j0eigen_descr0_3.chk j0eigen_descr0_4.chk j0eigen_descr1_1.chk j0eigen_descr1_2.chk j0eigen_descr1_3.chk j0eigen_descr1_4.chk j0eigen_descr2_1.chk j0eigen_descr2_2.chk j0eigen_descr2_3.chk j0eigen_descr2_4.chk j0eigen_descr3_1.chk j0eigen_descr3_2.chk j0eigen_descr3_3.chk j0eigen_descr3_4.chk j0eigen_descr4_1.chk j0eigen_descr4_2.chk j0eigen_descr4_3.chk j0eigen_descr4_4.chk j0eigen_descr10_1.chk j0eigen_descr10_2.chk j0eigen_descr10_3.chk j0eigen_descr10_4.chk j0eigen_descr11_1.chk j0eigen_descr11_2.chk j0eigen_descr11_3.chk j0eigen_descr11_4.chk j0eigen_descr12_1.chk j0eigen_descr12_2.chk j0eigen_descr12_3.chk j0eigen_descr12_4.chk" 

let return_sum=0

pipenv run python compare_results.py --kind quantum --folder1 "$folder1" --folder2 "$folder2" $quantum_files

for f in $column_files; do
  if [[ "$f" == "external.chk" ]]; then
    pipenv run python compare_results.py --kind column --column 3 --precision 5e-3 --folder1 "$folder1" --folder2 "$folder2" $f
  elif [[ "$f" == "potential.chk" ]]; then
    diff -u <(extract_column $folder1/$f 3) <(extract_column $folder2/$f 3)
  elif [[ "$f" == "kinetic.chk" ]]; then
    diff -u <(extract_column $folder1/$f 5) <(extract_column $folder2/$f 5)
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
