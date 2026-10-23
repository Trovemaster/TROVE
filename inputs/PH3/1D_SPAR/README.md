# PH3, 1D SPAR kinetic energy operator (Taylor series in s4, s5)

Full-dimensional (6D) PH3 in the SPAR coordinates `R-SYMPHI-BETA` (r1, r2, r3, s4,
s5, beta), C3v(M), polyad cut-off 6. The kinetic energy operator is read from
`kinetic.chk`: exact in r1, r2, r3 and beta, and a Taylor series in the bending
pair (s4, s5) to total degree 8. Every term is a product of 1D functions, so all
matrix elements are products of 1D integrals. The potential (MORBID,
`poten_xy3_morbid_10`) is expanded by TROVE to `PotOrder 8` in step 1.

The same calculation with the exact KEO is in `../2D_exact`; with N = 8 the two
agree to about 0.002 cm-1 below 3500 cm-1 and 0.1 cm-1 near 7000 cm-1.

## Files

| File | Purpose |
| :--- | :--- |
| `PH3_P6_J0.inp` | Step 1, J = 0: potential expansion, primitive and contracted basis sets, vibrational matrix elements (`matelem save split`). |
| `PH3_P6_J1.inp` | Step 2, J = 1: reads the checkpoints of step 1 (`matelem read split`). For other J change `Jrot`. |
| `kinetic.chk` | 1D SPAR KEO, Taylor degree N = 8 in (s4, s5) (`BASIC-FUNCTION mass` format). |
| `generate_xy3_spar_keo.py`, `make_spar_kinetic.py` | Regenerate `kinetic.chk`: `python3 make_spar_kinetic.py 8 kinetic.chk` (other N for convergence tests). |
| `run.sh` | Runs the two steps. |

Do not put a `keo_2d_s4s5.chk` into this folder: TROVE would then use 2D integrals.

## Running

Copy the folder to a working directory and run the steps in the order listed
above (or `./run.sh`; set `TROVE_EXE` to the TROVE executable, otherwise `j-trove.x`
is taken from the PATH):

```bash
cp -r inputs/PH3/1D_SPAR /path/to/work && cd /path/to/work/1D_SPAR
./run.sh
```

## Reference energies

In cm-1, relative to the zero-point energy of 5214.173918 cm-1:

| J | Symmetry | State | Energy |
| :---: | :---: | :--- | ---: |
| 0 | A1 | nu2 | 992.4200 |
| 0 | E | nu4 | 1118.4932 |
| 0 | A1 | 2nu2 | 1973.8687 |
| 0 | E | nu2 + nu4 | 2109.0115 |
| 0 | A1 | 2nu4 | 2227.2677 |
| 0 | E | 2nu4 | 2235.4942 |
| 0 | A1 | nu1 | 2321.6121 |
| 0 | E | nu3 | 2327.2029 |
| 1 | E | ground state, K = 1 | 8.3712 |
| 1 | A2 | ground state, K = 0 | 8.9045 |
| 1 | E | nu2, K = 1 | 1000.6400 |
| 1 | A2 | nu2, K = 0 | 1000.9661 |
| 1 | E | nu4 | 1123.3477 |
| 1 | E | nu4 | 1127.6002 |
| 1 | A1 | nu4 | 1130.3633 |
| 1 | A2 | nu4 | 1130.7114 |
