## TROVE H2CO Benchmark

This is a set of sample inputs and outputs for the software [TROVE](https://github.com/Trovemaster/TROVE) and is intended for use as a regression test.

In order to run this manually:
1. `cp input/*inp <output_folder>`
2. `cp run_benchmark.sh <output_folder>`
3. Copy the trove exe to `<output_folder>`
4. `cd <output_folder>; ./run_benchmark.sh`
5. `cd ..`
6. `./compare_results.sh <output_folder> gfortran`
7. `./compare_results.sh <output_folder> intel`
