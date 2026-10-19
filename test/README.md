# TROVE H2CO Regression Benchmark

This directory contains the regression test suite for [TROVE](https://github.com/Trovemaster/TROVE), centred on the **13-step $\text{H}_2\text{CO}$ (formaldehyde) benchmark**.

---

## 1. Does it require multiple calls?

**Yes.** The $\text{H}_2\text{CO}$ benchmark consists of **13 sequential calls** to the TROVE executable (`j-trove.x`) executed in the **same working directory**. Later steps read `.chk` checkpoint files produced by earlier steps:

| Step | Input File | Purpose | Key Checkpoints Read / Written | Runtime (1 core) | Peak Memory |
| :---: | :--- | :--- | :--- | :---: | :---: |
| **1** | `file1.inp` | Builds 1D Numerov primitive basis, kinetic (`KinOrder 2`), potential (`PotOrder 4`), and dipole expansions, contracted basis (`Npolyads 4`), and $J=0$ vibrational states & matrix elements. | **Writes:** `kinetic.chk`, `potential.chk`, `external.chk`, `prim_bset.chk`, `numerov_bset.chk`, `contr_*.chk`, `contr_matelem.chk`, `contr_extfield.chk`, `eigen_*0*.chk` | ~8.3 s | ~21 MB |
| **2** | `file2.inp` | Converts the contracted Hamiltonian and dipole matrix elements into the **$J=0$ representation** (`model j=0`, `matelem convert`, `extmatelem convert`). | **Reads:** Step 1 `.chk` files<br>**Writes:** `j0contr_*.chk`, `j0_matelem.chk`, `j0_extfield.chk` | ~0.5 s | ~15 MB |
| **3–7** | `file3.inp` – `file7.inp` | Diagonalizes the ro-vibrational Hamiltonian in the $J=0$ basis for **$J = 0, 1, 2, 3, 4$** respectively. | **Reads:** `j0_matelem.chk`, `j0contr_*.chk`<br>**Writes:** `j0eigen_quanta{0..4}.chk`, `j0eigen_descr{0..4}_*.chk`, `j0eigen_vectors{0..4}_*.chk` | ~0.5–0.7 s each | ~15–17 MB |
| **8–9** | `file8.inp`, `file9.inp` | Diagonalizes the ro-vibrational Hamiltonian for **$J = 10$**. | **Writes:** `j0eigen_*10*.chk` | ~1.7 s each | ~17 MB |
| **10** | `file10.inp` | Builds and saves the Hamiltonian matrix for **$J = 20$** ($\Gamma = 1$, `DIAGONALIZER save`). | **Writes:** `matrix20_1.chk`, `j0eigen_*20*.chk` | ~2.4 s | ~25 MB |
| **11–12** | `file11.inp`, `file12.inp` | Diagonalizes the ro-vibrational Hamiltonian for **$J = 11$** and **$J = 12$**. | **Writes:** `j0eigen_*11*.chk`, `j0eigen_*12*.chk` | ~2.0–2.3 s each | ~20–21 MB |
| **13** | `file_intensity.inp` | Computes ExoMol absorption line intensities for **$J = 0 \dots 4$** ($T = 1000\text{ K}$, $0\text{--}2000\text{ cm}^{-1}$). | **Reads:** `j0_extfield.chk` + `j0eigen_*` ($J=0..4$) | ~0.6 s | ~1.2 MB |

---

## 2. Runtime, Memory, and Disk Requirements

- **Total Wall-Clock Time:** **~24 seconds** on 1 CPU core (`OMP_NUM_THREADS=1`).
- **Peak RAM (Resident Set Size):** **~57 MB RSS** (TROVE internal array allocation peaks at **~25 MB** during `file10.inp`).
  - *Note:* The first line of `file2.inp`–`file_intensity.inp` specifies `mem 64 gb` (`mem 2 gb` in `file1.inp`). This is only an upper-bound limit for TROVE's internal memory counter—TROVE does **not** allocate 64 GB of RAM.
- **Disk Space:** **~71 MB** in the run directory (`outputs/H2CO/`) for all generated `.chk` and `.out` files.

---

## 3. How to Run

### Option A: Automated (Recommended)

1. Compile `j-trove.x` in the top-level `TROVE/` directory:
   ```bash
   make
   ```
2. Ensure Python 3 + `pytest`, `curl`, and `jq` are installed, then run from `test/`:
   ```bash
   cd test
   ./run_regression_tests.sh
   ```
   This script automatically:
   - Downloads and extracts `benchmarks.tar.gz` from Zenodo (`10.5281/zenodo.4546263`) if `benchmarks/H2CO` is not yet present.
   - Copies `j-trove.x`, `benchmarks/H2CO/input/*.inp`, and `scripts/H2CO/run_benchmark.sh` into `outputs/H2CO/`.
   - Executes all 13 steps in sequence via `run_benchmark.sh`.
   - Compares the resulting energies and expansion coefficients against `benchmarks/H2CO/outputs/` using `compare_results.py` and prints `PASS: H2CO`.

### Option B: Manual Step-by-Step Execution

To run the 13 calls manually in a scratch directory (e.g. `run_H2CO/`):

```bash
mkdir -p run_H2CO
cp benchmarks/H2CO/input/*.inp run_H2CO/
cp ../j-trove.x run_H2CO/
cd run_H2CO

export OMP_NUM_THREADS=1
ulimit -d unlimited

# Run the 13 sequential steps in order:
for name in file{1..12} file_intensity; do
  ./j-trove.x < ${name}.inp > ${name}.out
done
```

To verify the outputs against the reference benchmark data afterward:
```bash
cd ..
bash scripts/H2CO/compare_results.sh run_H2CO benchmarks/H2CO/outputs
```
