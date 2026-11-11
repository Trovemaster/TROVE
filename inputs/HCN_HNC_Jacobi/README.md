# HCN/HNC isomerising system: exact Jacobi KEO, RELAXED two-well basis, SPAR PES and Jacobi-frame DMS

A single TROVE model covering both HCN and HNC, including isomerisation, built on a new
ab initio PES and DMS: frozen-core RHF-CCSD(T)/aug-cc-pVQZ, about 9000 geometries in 1D/2D/3D
grids along the isomerisation path, plus large-distortion "guard" points; the DMS uses finite
fields (h = 0.002 a.u.).

## Model

* **Coordinates:** Jacobi r = r(CN), R = |H − CM(CN)| and γ (0 = HCN, 180° = HNC), with
  `TRANSFORM/FRAME R1-Z-R2-ALPHA` and `REFER-CONF NON-RIGID`. The frame has z along N→C and
  H in the xz plane on the +x side.
* **KEO:** `KINETIC_XYZ_EKE_JACOBI_SINRHO`, the exact Jacobi KEO with the γ = 0 and γ = π
  singularities removed analytically (sin ρ formalism).
* **PES:** `POTEN_XYZ_JACOBI_MORSE_COS_MEP`, a Morse × cos(mγ) expansion referenced to the
  minimum-energy path (MEP):
  V = Σ F(p1,p2,m) u1^p1 u2^p2 cos(mγ), with u_k = 1 − exp(−a_k (r_k − r_k^MEP(γ))),
  a = (1.3, 1.0) Å⁻¹, p1 + p2 ≤ 8 and m ≤ 30. That is 747 parameters, plus cos-series for
  r^MEP(γ) and R^MEP(γ) (rows −1/−2).
  * Fit RMS: 0.16 / 0.13 / 0.21 / 5.9 cm⁻¹ below 5000 / 10000 / 20000 / 30000 cm⁻¹ (8823 points).
  * HNC lies 5171.86 cm⁻¹ above HCN; the isomerisation barrier is 16 660.4 cm⁻¹ at γ = 76°.
  * There are no holes in the region covered by the TROVE stretching basis.
  * `potential.chk` is the exact SPAR conversion (sum of 1D products in
    y_k = 1 − exp(−a_k (r_k − r_k^ref)) times functions of γ). It is read with
    `POTENTIAL READ WITH_MODES`, and the `EQUILIBRIUM` values in the inputs are its reference
    geometry.
* **DMS:** `DIPOLE_XYZ_JACOBI_MEP`, built in the same u variables:
  μ_z = Σ Z u1^i u2^j cos(kγ) and μ_x = Σ X u1^i u2^j sin(kγ), with i + j ≤ 6
  (400 + 372 parameters). The RMS is about 2×10⁻⁴ D below 20 000 cm⁻¹. Rows 91, 92 and 93
  carry r^MEP, R^MEP and the Morse exponents. The block needs `COORDS linear linear linear`.
* **Basis (new):**
  * `numerov` for r (class 1).
  * `numerov` for R and `sinrho-2xlaguerre-k` for γ, contracted together as class 2.
  * The new `BASIS` keyword `RELAXED` on R and γ:
    * on a stretch, the Numerov basis is generated on the envelope min over ρ of V(x; ρ);
    * on the non-rigid mode, the 1D bending problem is solved along the minimum-energy path:
      at each ρ the stretches are moved to the minimum of V, and V, g_vib, g_rot(z,z) and the
      pseudo-potential are evaluated there.
  * Use ≥ 3000 γ grid points; with 1000 points, bending levels are about 0.05 cm⁻¹ high per quantum.

Without `RELAXED` and the merged (R, γ) class, the HNC levels are badly unconverged, by 35 to
several hundred cm⁻¹, because the R and γ functions are then built around the HCN geometry.
With them, TROVE reproduces an independent exact J = 0 DVR (`dvr_j0_check.py`, sinc × sinc ×
Legendre DVR) to better than 0.1 cm⁻¹ for all 61 HCN and HNC states below 10 400 cm⁻¹. It also
reproduces the J = 0 dipole matrix elements to 2×10⁻⁴ D.

## Files

| File | Content |
|---|---|
| `HCN_HNC_J0_step1.inp/.out` | J = 0 energies. Saves basis, contraction, matrix elements, the DMS field (`external.chk`) and eigenvectors |
| `HCN_HNC_J0_step2_tm.inp/.out` | J = 0 vibrational transition moments (`INTENSITY tm`), read from the step-1 checkpoints |
| `HCN_HNC_J0.inp`, `HCN_HNC_J1.inp`, `HCN_HNC_J2.inp` (+ `.out`) | J = 0 (kmax = 2, `MATELEM SAVE`), 1 and 2: rotational and ℓ-doubling constants of both isomers |
| `potential.chk` | SPAR form of the PES (`POTENTIAL READ WITH_MODES`) |
| `external.chk` | DMS expansion written by step 1 (`EXTERNAL SAVE`); use `EXTERNAL READ` to reuse it |
| `HCN_CCSDT_AVQZ_poten.inp_block`, `HCN_CCSDT_AVQZ_dipole.inp_block` | Ab initio PES and DMS parameter blocks (already included in the inputs) |
| `HCN_CCSDT_AVQZ_TROVE_refined_k8_staged_poten.inp_block`, `potential_refined_k8_staged.chk` | Staged $J \le 20$ ($K \le 8$) empirically refined PES block and SPAR `potential.chk` (2,653 MARVEL/Mellau HCN + HNC levels, $\text{RMS} = 0.193\text{ cm}^{-1}$) |
| `HCN_HNC_k8_J0_20_stage4_147.inp/.out` | $J = 0..20$ ($K \le 8$) empirical PES refinement input and output |
| `HCN_HNC_J0_k4_ref_step1.inp`, `HCN_HNC_J0_k4_ref_step2_tm.inp/.out`, `trove_j0_ref_G10.xsec` | $J = 0$ ($K \le 4$) parallel ($\Delta\ell = 0$) and perpendicular ($\Delta\ell = \pm 1$) transition moment inputs/outputs and $T = 300\text{--}2000\text{ K}$ cross sections (`G10`) using the refined PES + ab initio DMS |
| `dvr_j0_check.py`, `dvr_J0_reference.out` | Independent J = 0 DVR check (NumPy/SciPy) and its converged output |
| `run.sh` | Runs the whole sequence |

Resources:
- Step 1 takes about 2.5 min on 12 cores.
- The J = 0 (kmax = 2) run takes about 9 min and writes a 29 GB `contr_matelem.chk`.
- J = 1 takes about 2 min; J = 2 takes about 5.5 min and needs about 27 GB of RAM.

## Reference results for this PES

J = 0 term values (cm⁻¹; HNC states relative to the HNC ground state):

| State | TROVE | Experiment |
|---|---|---|
| ZPE | 3466.33 | |
| HCN 2ν₂⁰ | 1407.84 | 1411.41 |
| HCN ν(C≡N) | 2091.31 | 2096.85 |
| HCN 4ν₂⁰ | 2797.96 | 2802.96 |
| HCN ν(C–H) | 3305.13 | 3311.48 |
| HNC ground state (above HCN ground) | 5063.89 | |
| HNC 2ν₂⁰ | 921.87 | |
| HNC 4ν₂⁰ | 1859.02 | |
| HNC ν(N≡C) | 2015.09 | 2023.86 |
| HNC ν(N–H) | 3650.13 | 3652.66 |

J = 0 dipole matrix elements (D):

| Quantity | Value |
|---|---|
| HCN ⟨0\|μ_z\|0⟩ | 2.9783 (exp. 2.985) |
| HCN ν(C≡N) | 0.00094 |
| HCN ν(C–H) | 0.0834 |
| HNC ⟨0\|μ_z\|0⟩ | −2.9807 |
| HNC ν(N≡C) | 0.1096 |
| HNC ν(N–H) | 0.1497 |

Rotational and ℓ-doubling constants from `HCN_HNC_J{0,1,2}.out` (cm⁻¹):

| Quantity | HCN (TROVE) | HCN (exp.) | HNC (TROVE) | HNC (exp.) |
|---|---|---|---|---|
| B₀ | 1.47004 | 1.47822 | 1.50333 | 1.51210 |
| D₀ | 2.88×10⁻⁶ | 2.91×10⁻⁶ | 3.25×10⁻⁶ | 3.3×10⁻⁶ |
| G(ν₂¹) | 710.95 | 711.98 | 461.26 | 462.72 |
| B(ν₂) − B₀ | +0.00354 | +0.00355 | +0.00502 | |
| q(ν₂) | 0.007420 | 0.007488 | 0.010740 | |

Frozen-core aug-cc-pVQZ gives stretching fundamentals 2.5–9 cm⁻¹ low and B₀ about 0.55% low,
because r(CN) is 0.0035 Å too long. The empirically refined PES (`potential_refined_k8_staged.chk`,
fitted to 2,653 MARVEL/Mellau levels up to $J = 20, K = 8$) reduces the overall ro-vibrational RMS
to $0.193\text{ cm}^{-1}$ ($0.196\text{ cm}^{-1}$ for 1,543 HCN levels and $0.189\text{ cm}^{-1}$
for 1,110 HNC levels) and reproduces ExoMol `1H-12C-14N__Harris` integrated opacities across
$T = 300\text{--}1500\text{ K}$ to within $0.1\text{--}7.6\%$ while correcting the $+492.8\text{ cm}^{-1}$
HNC isomerisation energy offset ($5212.22\text{ cm}^{-1}$ vs. $5705.00\text{ cm}^{-1}$ in `Harris`).
