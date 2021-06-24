#!/usr/bin/env bash

set -e

exe_name=j-trove.x
exe=../$exe_name
# Use 1 process unless we have specified differently (e.g. in CI)
nproc=${nproc:-1}

# Check exe is present
if [ ! -f $exe ]; then
  echo "ERROR: $exe not found."
  exit 2
fi

# If benchmarks haven't been downloaded, download them
if [ ! -d benchmarks ]; then
  echo "Downloading benchmarks"
  ./download_benchmarks.sh
  tar xzf benchmarks.tar.gz
fi

# Run series of benchmarks
for benchmark in H2CO; do
  # Create output folder & move necessary files
  wd=outputs/$benchmark
  mkdir -p $wd
  cp $exe $wd
  cp benchmarks/$benchmark/input/*.inp $wd
  cp scripts/$benchmark/run_benchmark.sh $wd

  # Run benchmark
  pushd $wd
  echo "Running $benchmark"
  ./run_benchmark.sh $nproc $exe_name
  popd

  # Compare results with "truth"
  bash scripts/$benchmark/compare_results.sh $wd benchmarks/$benchmark/outputs
  success=$?
  if [ $success -eq 0 ]; then
    echo "PASS: $benchmark"
  else
    echo "FAIL: $benchmark failed with exit code $success"
    exit $success
  fi
done
