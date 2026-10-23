# C2H4, SPAR kinetic energy operator and potential

Full-dimensional (12D) C2H4 in the `C2H4_SPAR` coordinates, D2h(M), polyad cut-off
4 for the primitive and contracted basis sets. Both the kinetic energy operator
and the potential energy function are read in SPAR (sum-of-products) form:

- `kinetic.chk`: KEO (`BASIC-FUNCTION mass`, `kinetic read with_modes`);
- `potential.chk`: potential (`BASIC-FUNCTION POTEN`, `potential read with_modes`),
  a SPAR representation of `poten_c2h4_886666`; the `POTEN` parameters in the
  inputs are not used for the matrix elements.

## Files

| File | Purpose |
| :--- | :--- |
| `C2H4_SPAR_P4_J0.inp` | Step 1, J = 0: primitive and contracted basis sets, vibrational matrix elements (`matelem save split`) and J = 0 energies. |
| `C2H4_SPAR_P4_J0_convert.inp` | Step 2: converts the matrix elements to the J = 0 contracted representation (`model j=0`, `matelem convert split`). |
| `C2H4_SPAR_P4_J0_j0.inp` | Step 3, J = 0 in the J = 0 representation; it provides the zero-point energy to which the J > 0 energies are referred. |
| `C2H4_SPAR_P4_J1_j0.inp` | Step 4, J = 1 in the J = 0 representation. For other J change `Jrot`. |
| `kinetic.chk`, `potential.chk` | SPAR KEO and potential. |
| `run.sh` | Runs the four steps. |

The polyad cut-off is set by `Npolyads` in `PRIMITIVES` and `CONTRACTION`; the same
inputs with `Npolyads 6` give a larger calculation (about 7 min for step 1 on 24
cores and several GB of checkpoint files).

## Running

Copy the folder to a working directory and run the steps in the order listed
above (or `./run.sh`; set `TROVE_EXE` to the TROVE executable, otherwise `j-trove.x`
is taken from the PATH):

```bash
cp -r inputs/C2H4/SPAR_P4 /path/to/work && cd /path/to/work/SPAR_P4
./run.sh
```

Step 1 takes about 53 s on 24 cores; the other steps take seconds.

## Reference energies

Zero-point energy 11011.285848 cm-1. Lowest levels in cm-1, relative to the
zero-point energy:

| J | Symmetry | # | Energy |
| :---: | :---: | ---: | ---: |
| 0 | Ag | 1 | 0.0000 |
| 0 | B2u | 1 | 838.9611 |
| 0 | B2g | 1 | 944.2412 |
| 0 | B3u | 1 | 954.7722 |
| 0 | Au | 1 | 1029.8399 |
| 0 | B3g | 1 | 1248.4054 |
| 0 | Ag | 2 | 1361.1891 |
| 0 | B1u | 1 | 1453.8578 |
| 0 | Ag | 3 | 1640.2279 |
| 0 | Ag | 4 | 1742.9784 |
| 0 | Au | 2 | 1839.7559 |
| 0 | B1g | 1 | 1853.0669 |

| J | Symmetry | # | Energy |
| :---: | :---: | ---: | ---: |
| 1 | B1g | 1 | 1.8344 |
| 1 | B2g | 1 | 5.6903 |
| 1 | B3g | 1 | 5.8577 |
| 1 | B3u | 1 | 840.8391 |
| 1 | Au | 1 | 844.5206 |
| 1 | B1u | 1 | 844.6732 |
| 1 | B3g | 2 | 946.1830 |
| 1 | Ag | 1 | 949.8620 |
