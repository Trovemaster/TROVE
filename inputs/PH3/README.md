# PH3 with the exact (s4, s5) kinetic energy operator

This example computes ro-vibrational energies of PH3 with TROVE in the SPAR
coordinates `R-SYMPHI-BETA` (r1, r2, r3, s4, s5, beta) using an exact kinetic
energy operator (KEO). In these coordinates every term of G_vib, G_rot, G_cor
and of the pseudo-potential U is a product of 1D functions of r1, r2, r3 and
beta and of a 2D function of the bending pair (s4, s5). The 2D functions are not
expanded in powers of s4 and s5: their matrix elements in the harmonic-oscillator
basis of modes 4 and 5 are computed by 2D quadrature and read by TROVE from
`keo_2d_s4s5.chk`.

## Files

| File | Purpose |
| :--- | :--- |
| `exact_2d_keo.py` | Writes `kinetic.chk` (exact KEO in the `BASIC-FUNCTION mass` format) and `keo_2d_s4s5.chk` (2D integrals). Needs `numpy` and `scipy`. |
| `PH3_exact2D_P6_J0.inp` | Step 1, J = 0: computes the potential expansion, the primitive and contracted basis sets and saves the vibrational matrix elements (`matelem save split`). |
| `PH3_exact2D_P6_J1.inp` | Step 2, J = 1: reads the checkpoints of step 1 (`matelem read split`). For other J, change `Jrot`. |

## Running

Run the three commands in the same directory:

```bash
python3 exact_2d_keo.py --mode exact --vmax 6 --out-dir .
j-trove.x < PH3_exact2D_P6_J0.inp > PH3_exact2D_P6_J0.out
j-trove.x < PH3_exact2D_P6_J1.inp > PH3_exact2D_P6_J1.out
```

TROVE uses the 2D integrals whenever `keo_2d_s4s5.chk` is present in the working
directory and then prints

```
 Loaded exact 2D (s4,s5) KEO integrals from keo_2d_s4s5.chk: kmax=  8, vmax=  6
```

## Keeping the script and the TROVE input consistent

- `--vmax` must be at least the basis-set range of modes 4 and 5 (`range 0, 6`).
- `--rmax` is the radius of the 2D quadrature; use the `borders` of modes 4 and 5
  (2.45).
- `--ft` is the length scale of the harmonic-oscillator functions of modes 4 and
  5. It must be the first number of the `params` line that TROVE prints for the
  HARMONIC basis of modes 4 and 5:

  ```
   params  =         0.28658334     1145.16544179        0.00000000
  ```

  The value follows from the potential and the KEO at equilibrium (`AUTO-HARMON`).
  The default, 0.2865833376039929, belongs to this input. With the printed
  8-digit value the energies change by less than 1e-4 cm-1.
- Keep the `BASIC-FUNCTION` section, `KinOrder` and the classes of the `BASIS`
  section as in the example. Modes 4 and 5 must form a class of their own.
- In the exact `kinetic.chk` the indices (i4, i5) of the 2D terms label 2D basis
  functions, not powers of s4 and s5, so this `kinetic.chk` is meaningful only
  together with its `keo_2d_s4s5.chk`.

TROVE stops with an error message if modes 4 and 5 are not a separate class, or
if `keo_2d_s4s5.chk` does not cover the basis set of modes 4 and 5 or the (i4, i5)
indices of the G-matrix terms.

## Format of keo_2d_s4s5.chk

The file holds the integrals

    me2d(iop4, iop5, i4, i5, v4, v4', v5, v5') =
        Int ds4 ds5 D(iop4)[v4,v4'](s4) D(iop5)[v5,v5'](s5) phi(i4,i5)(s4,s5)

where phi(i4,i5) are the 2D functions of the KEO terms (phi(0,0) = 1), and
D(0) = psi_v psi_v', D(1) = psi_v psi_v'', D(2) = -psi_v' psi_v', D(3) = -psi_v' psi_v''
(psi_v are the harmonic-oscillator functions of modes 4 and 5; a prime on psi denotes its
derivative, so iop = 1, 2, 3 places the momentum on the right, left or both functions).

By default the file is ASCII: a keyword line starting with `keo_2d_s4s5`, a line with
kmax, vmax and the number of entries N, a line describing how the table was made, a
column-header line, N lines `iop4 iop5 i4 i5 v4 v4' v5 v5' value` (the non-zero
elements, 17 significant digits) and the line `End of keo_2d_s4s5`. The ASCII file
for this example is about 83 MB. `--binary` writes a smaller binary stream file
instead (4-byte integers kmax, vmax followed by the array in Fortran order, about
25 MB). TROVE recognises either format.

`--mode taylor --kmax N` writes the 2D integrals of the monomials s4^i4 s5^i5
(i4 + i5 <= N) instead. With a Taylor-expanded (SPAR) `kinetic.chk` of order N
this reproduces the standard 1D x 1D integration of TROVE exactly, which is a
useful check of the 2D integrals.

## Reference energies

With the inputs as given (polyad cut-off 6), in cm-1, relative to the zero-point
energy of 5214.173894 cm-1:

| J | Symmetry | State | Energy |
| :---: | :---: | :--- | ---: |
| 0 | A1 | nu2 | 992.4200 |
| 0 | E | nu4 | 1118.4931 |
| 0 | A1 | 2nu2 | 1973.8687 |
| 0 | E | nu2 + nu4 | 2109.0113 |
| 0 | A1 | 2nu4 | 2227.2674 |
| 0 | E | 2nu4 | 2235.4937 |
| 0 | A1 | nu1 | 2321.6121 |
| 0 | E | nu3 | 2327.2029 |
| 1 | E | ground state, K = 1 | 8.3712 |
| 1 | A2 | ground state, K = 0 | 8.9045 |
| 1 | E | nu2, K = 1 | 1000.6400 |
| 1 | A2 | nu2, K = 0 | 1000.9660 |
| 1 | E | nu4 | 1123.3475 |
| 1 | E | nu4 | 1127.6001 |
| 1 | A1 | nu4 | 1130.3632 |
| 1 | A2 | nu4 | 1130.7113 |
