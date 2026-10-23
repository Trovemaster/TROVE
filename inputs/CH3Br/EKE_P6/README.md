# CH3Br, exact kinetic energy operator

Full-dimensional (9D) CH3Br in the `X-R-BETA-SYM` coordinates, C3v(M), polyad
cut-off 6. The exact kinetic energy operator is read from `kinetic.chk`
(`BASIC-FUNCTION mass`, mass-factorised SPAR form). The potential
(`poten_zxy3_sym`, parameters in the inputs) is expanded by TROVE to
`PotOrder 4` in step 1 and saved in `potential.chk`, which the later steps read.

## Files

| File | Purpose |
| :--- | :--- |
| `CH3Br_EKE_P6_J0.inp` | Step 1, J = 0: potential expansion, primitive and contracted basis sets, vibrational matrix elements (`matelem save split`) and J = 0 energies. |
| `CH3Br_EKE_P6_J0_convert.inp` | Step 2: converts the matrix elements to the J = 0 contracted representation (`model j=0`, `matelem convert split`). |
| `CH3Br_EKE_P6_J0_j0.inp` | Step 3, J = 0 in the J = 0 representation; it provides the zero-point energy to which the J > 0 energies are referred. |
| `CH3Br_EKE_P6_J1_j0.inp` | Step 4, J = 1 in the J = 0 representation. For other J change `Jrot`. |
| `kinetic.chk` | Exact KEO. |
| `run.sh` | Runs the four steps. |

## Running

Copy the folder to a working directory and run the steps in the order listed
above (or `./run.sh`; set `TROVE_EXE` to the TROVE executable, otherwise `j-trove.x`
is taken from the PATH):

```bash
cp -r inputs/CH3Br/EKE_P6 /path/to/work && cd /path/to/work/EKE_P6
./run.sh
```

Step 1 takes about 2 min and step 2 about 2 min on 24 cores; the other steps take seconds.

## Reference energies

Zero-point energy 8063.099944 cm-1. Lowest levels in cm-1, relative to the
zero-point energy:

| J | Symmetry | # | Energy |
| :---: | :---: | ---: | ---: |
| 0 | A1 | 1 | 0.0000 |
| 0 | A1 | 2 | 614.5575 |
| 0 | E | 1 | 957.0555 |
| 0 | A1 | 3 | 1229.0720 |
| 0 | A1 | 4 | 1311.1521 |
| 0 | E | 2 | 1445.7313 |
| 0 | E | 3 | 1575.1246 |
| 0 | A1 | 5 | 1874.1061 |
| 0 | A1 | 6 | 1918.6975 |
| 0 | E | 4 | 1925.6918 |
| 0 | A1 | 7 | 1932.2845 |
| 0 | E | 5 | 2065.5022 |

| J | Symmetry | # | Energy |
| :---: | :---: | ---: | ---: |
| 1 | A2 | 1 | 0.6324 |
| 1 | E | 1 | 5.4833 |
| 1 | A2 | 2 | 615.1865 |
| 1 | E | 2 | 620.0334 |
| 1 | E | 3 | 957.6873 |
| 1 | A2 | 3 | 960.2412 |
| 1 | A1 | 1 | 960.2412 |
| 1 | E | 4 | 964.8982 |
