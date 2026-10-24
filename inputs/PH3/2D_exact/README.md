# PH3 with the exact 2D (s4, s5) SPAR kinetic energy operator

This example computes ro-vibrational energies of PH3 with TROVE in the SPAR
coordinates `R-SYMPHI-BETA` (r1, r2, r3, s4, s5, beta) using an **exact** kinetic
energy operator (KEO) defined entirely through `BASIC-FUNCTION mass` and
`kinetic.chk`.

## How the 2D `BASIC-FUNCTION` (`DELTA_XY3`) works

In the `R-SYMPHI-BETA` coordinates, every term of $G_{\text{vib}}$, $G_{\text{rot}}$,
$G_{\text{cor}}$, and the pseudopotential $U$ is an exact finite linear combination
of products of:
- **Modes 1–3 ($r_1, r_2, r_3$)**: 1D inverse-bond-length functions $r_i^{-1}, r_i^{-2}$,
- **Mode 6 ($\beta$)**: 1D trigonometric functions $\csc^2\beta, \cot\beta, \sin\beta\cos\beta, \sin^2\beta$,
- **Mode 4 ($s_4$)**: 1D trigonometric harmonics $\cos(k_4 a_4 s_4)$ and $\sin(k_4 a_4 s_4)$ ($a_4 = 1/\sqrt{6}$, $k_4 = 1\dots 10$),
- **Mode 5 ($s_5, s_4$)**: 2D functions $\cos(k_5 a_5 s_5)\,\Delta(s_4, s_5)^{-d}$ and $\sin(k_5 a_5 s_5)\,\Delta(s_4, s_5)^{-d}$ ($a_5 = 1/\sqrt{2}$, $d = 0\dots 4$), using the 2D elementary constructor `DELTA_XY3`:
  $$\Delta(s_4, s_5) = \cos\left(\frac{s_4}{\sqrt{6}} + \frac{s_5}{\sqrt{2}}\right) + \cos\left(\frac{s_4}{\sqrt{6}} - \frac{s_5}{\sqrt{2}}\right) + \cos\left(\frac{2 s_4}{\sqrt{6}}\right)$$

Because Mode 5's `BASIC-FUNCTION` entries contain `DELTA_XY3`, TROVE automatically
detects the 2D $(s_4, s_5)$ coupling (`molec%basic_2d_s4s5 = .true.`) and computes
all required 2D harmonic-oscillator matrix elements (`me2d_s4s5`) on the fly via
OpenMP-parallelized polar Gauss–Legendre $\times$ Fourier quadrature in `PTcompute_me2d_s4s5`
(~3.8 s wall time for $v_{\max} = 6$). No external `keo_2d_s4s5.chk` file is required.

The same calculation with the 1D Taylor-expanded SPAR KEO is in `../1D_SPAR`.

## Files

| File | Purpose |
| :--- | :--- |
| `PH3_P6_J0.inp` | Step 1, $J = 0$: potential expansion, primitive and contracted basis sets, vibrational and rotational reduced matrix elements (`matelem save split`). |
| `PH3_P6_J1.inp` | Step 2, $J = 1$: reads the checkpoints of step 1 (`matelem read split`). For other $J$, change `Jrot`. |
| `kinetic.chk` | Exact algebraic 6D KEO (`10,293` records, 700 KB, `BASIC-FUNCTION mass` format). |
| `exact_2d_keo.py` | Regenerates `kinetic.chk` and the `BASIC-FUNCTION mass` block algebraically in ~1.3 s (pure Python `math`, no external dependencies). |
| `run.sh` | Runs the $J = 0$ and $J = 1$ steps. |

## Running

Copy the folder to a working directory and run `./run.sh` (set `TROVE_EXE` to the
TROVE executable, otherwise `j-trove.x` is taken from the repository root or `PATH`):

```bash
cp -r inputs/PH3/2D_exact /path/to/work && cd /path/to/work/2D_exact
./run.sh
```

During the primitive matrix-element step, TROVE prints:

```
 Computed exact 2D (s4,s5) KEO integrals from BASIC-FUNCTION: kmax=( 20, 54), active pairs= 265, vmax=  6
```

## Reference energies

With the inputs as given (polyad cut-off $P = 6$), in $\text{cm}^{-1}$, relative to the zero-point
energy of `5214.173894 cm-1`:

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
