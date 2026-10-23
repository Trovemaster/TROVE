# CH3OH, exact kinetic energy operator and SPAR potential (reduced 3D model)

CH3OH in the `R-ALPHA132-THETA-S-TAU` coordinates (frame
`R-ALPHA-THETA-TAU-OH-CLOCKWISE`), symmetry group `C3v(M)-2`, non-rigid reference
configuration along the torsion. Only modes 10 and 11 (harmonic basis, 0-4) and
the torsion, mode 12 (Fourier basis, 0-2), are active; the other modes are kept
in their ground states (`range 0, 0`). Both operators are read in SPAR form:

- `kinetic.chk`: exact KEO (`BASIC-FUNCTION mass`);
- `potential.chk`: potential (`potential read with_modes`) for `poten_ch3oh_sym_IV`.

## Files

| File | Purpose |
| :--- | :--- |
| `CH3OH_EKE_J0.inp` | Step 1, J = 0: primitive and contracted basis sets, vibrational matrix elements (`matelem save split`) and J = 0 energies. |
| `CH3OH_EKE_J1.inp` | Step 2, J = 1: reads the checkpoints of step 1 (`matelem read split`). For other J change `Jrot`. |
| `kinetic.chk`, `potential.chk` | Exact KEO and SPAR potential. |
| `run.sh` | Runs the two steps. |

## Running

Copy the folder to a working directory and run the steps in the order listed
above (or `./run.sh`; set `TROVE_EXE` to the TROVE executable, otherwise `j-trove.x`
is taken from the PATH):

```bash
cp -r inputs/CH3OH/EKE /path/to/work && cd /path/to/work/EKE
./run.sh
```

## Reference energies

Zero-point energy 11270.712541 cm-1. Lowest levels in cm-1, relative to the
zero-point energy:

| J | Symmetry | # | Energy |
| :---: | :---: | ---: | ---: |
| 0 | A1 | 1 | 0.0000 |
| 0 | E | 1 | 10.4978 |
| 0 | E | 2 | 1467.1315 |
| 0 | A1 | 2 | 1478.8038 |
| 0 | E | 3 | 1481.0782 |
| 0 | A2 | 1 | 1481.8182 |
| 0 | A1 | 3 | 2929.5762 |
| 0 | E | 4 | 2938.1577 |
| 0 | E | 5 | 2951.5329 |
| 0 | A2 | 2 | 2956.6163 |
| 0 | E | 6 | 2958.5293 |
| 0 | A1 | 4 | 2959.2970 |

| J | Symmetry | # | Energy |
| :---: | :---: | ---: | ---: |
| 1 | A2 | 1 | 2.1007 |
| 1 | E | 1 | 12.6350 |
| 1 | E | 2 | 16.5270 |
| 1 | A1 | 1 | 23.3407 |
| 1 | A2 | 2 | 23.4012 |
| 1 | E | 3 | 51.2281 |
| 1 | E | 4 | 1469.2377 |
| 1 | A2 | 3 | 1480.9772 |
