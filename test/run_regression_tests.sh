#!/usr/bin/env bash

set -e

exe=../j-trove.x
nproc=1

# TODO check exe is present

# TODO quiet input or put it somewhere

# TODO upload to zenodo
#./download_benchmarks.sh
# unzip benchmarks
for benchmark in H2CO; do
  wd=outputs/$benchmark
  mkdir -p $wd
  cp $exe $wd
  cp benchmarks/$benchmark/input/*.inp $wd
  cp scripts/$benchmark/run_benchmark.sh $wd
  pushd $wd
  echo "Running $benchmark"
  ./run_benchmark.sh $nproc
  popd
  bash scripts/$benchmark/compare_results.sh $wd benchmarks/$benchmark/outputs
  # TODO actually process the result ( and handle set -e properly in case it fails)
done
