# C2H2 6D `singular_3D` Rovibrational & Transition-Intensity Examples in TROVE

This directory contains input files, pre-computed checkpoint files, and reference output files for 6D rovibrational energy and transition-intensity calculations of acetylene ($\text{C}_2\text{H}_2$) in `TROVE` using:
- **Exact polynomial $t$-$s$ kinetic energy operator** (`TROVE_EXT_T_S`, 4th-order expansion) with the singular $1/\sin^2\rho_i$ ($i=1,2$) terms marked by trailing underscores (`_`) in `kinetic.chk`.
- **Uncontracted 3D coupled singular bend–bend–torsion basis (`singular_3D`)** on $(s_4=\rho_1, s_5=\rho_2, s_6=\tau)$ constructed from associated Legendre basis functions (`sinrho-legendre-k-sing` for modes 4 and 5) and unoptimised half-integer Fourier torsional functions (`fourier-unoptimised` for mode 6 on $\tau \in [0, 4\pi]$, `borders 0., 720.0 deg`), enforcing the exact angular-momentum constraint $\{l_1, l_2\} = \{|k_\tau - K/2|,\, k_\tau + K/2\}$ at every basis function and matrix element.
- **Extended $\mathcal{D}_{2h}(\text{M})$ molecular symmetry** on $\tau \in [0, 4\pi]$ (`SYMGROUP D2h(M)`), where physical single-valued rovibrational wavefunctions ($E'\Psi = +\Psi$) span the four irreps `gamma 1 2 3 4` ($A_g, A_u, B_{1g}, B_{1u}$), while the double-valued irreps ($B_{2g}, B_{2u}, B_{3g}, B_{3u}$) are used only at $J=0$ to generate the $K=\text{odd}$ contracted vibrational blocks ($B_{3g}$ for $\Pi_g$, $B_{2u}$ for $\Pi_u$).
- **$K$-parametrised $J=0$ contracted basis (`kmax 2`)** that builds and stores all $K \le K_{\max}$ contracted 3D and 6D basis blocks during the $J=0$ step (`contr_*.chk`, `matelem*.chk`) for subsequent $J>0$ Coriolis-coupled rovibrational diagonalizations (`J=1, 2`) and dipole transition-moment contractions (`extmatelem*.chk`).
- **Cylindrically regular 6D Potential Energy Function (`POTEN_C2H2_MORSE_COSALPHA_COSNTAU-SIN-SIN`)**:
  $$V(\mathbf{q}) = \sum_{ijklmn} c_{ijklmn} f_1^i f_2^j f_3^k f_4^l f_5^m \cos(n\tau) \times (\sin\rho_1 \sin\rho_2)^n,$$
  where $f_1 = 1 - e^{-a_1(R_{\text{CC}}-R_{\text{CC}}^{\text{eq}})}$, $f_{2,3} = 1 - e^{-a_2(r_{1,2}-r_{\text{CH}}^{\text{eq}})}$, and $f_{4,5} = 1 - \cos\rho_{1,2}$ ($n = 0, \dots, 4$). The exact factor $(\sin\rho_1\sin\rho_2)^n$ guarantees cylindrical smoothness and resolution of torsional degeneracies whenever either $\text{C}-\text{C}-\text{H}$ angle approaches linearity ($\rho_1 \to 0$ or $\rho_2 \to 0$).
- **6D Bisector-Frame Dipole Moment Surface (`external.chk`) & `singular_3D` `is_extf` Intensities**:
  Supports $\Delta K = 0$ parallel ($\mu_z$, `imu = 3`, $B_{1u}$) and $\Delta K = \pm 1$ perpendicular ($\mu_x, \mu_y$, `imu = 1, 2`, $B_{2u} / B_{3u}$ with $\cos(\tau/2) / \sin(\tau/2)$) dipole matrix elements in `singular_3D`, reproducing exact Hönl–London $\times g_{\text{ns}}$ line-strength ratios (`1 : 3` para/ortho nuclear-spin weights `GNS 1.0 1.0 3.0 3.0 0.0 0.0 0.0 0.0`).

## Files

- **Test Morse-$\cos\alpha$ PES ($P_{\max}=12$, $K \le 2$)**:
  - `C2H2_sing3D_D2h_J0_kmax2.inp` / `.out`: $J=0$ ($K \le 2$, $P_{\max}=12$) calculation on `potential.chk` ($\text{ZPE} = 5765.537473\text{ cm}^{-1}$).
  - `C2H2_sing3D_D2h_J1.inp` / `.out`: $J=1$ ($K \le 1$, `gamma 1 2 3 4`) rovibrational calculation ($2B_0 = 2.346311\text{ cm}^{-1}$, $\Delta_{fe}(\nu_4, J=1) = +0.010538\text{ cm}^{-1}$, $\Delta_{fe}(\nu_5, J=1) = +0.009184\text{ cm}^{-1}$).
  - `C2H2_sing3D_D2h_J2.inp` / `.out`: $J=2$ ($K \le 2$, `gamma 1 2 3 4`) rovibrational calculation ($6B_0 - 24D_0 = 7.038840\text{ cm}^{-1}$, $\Delta_{fe}(J=2) / \Delta_{fe}(J=1) = 3.000$, $\Delta_{g/u}$ $K=2$ levels).
- **4,331-Point `CCSD(T)/CBS(aug-cc-pVTZ, aug-cc-pVQZ)` Ab Initio PES ($P_{\max}=12$, $K \le 2$)**:
  - `C2H2_CBS_PES_poten.inp_block` & `potential_CBS.chk`: 257-parameter (`NPARAM 264`, 466 expanded terms) 6D `POTEN_C2H2_MORSE_COSALPHA_COSNTAU-SIN-SIN` fit to 4,331 `CCSD(T)/CBS` geometries ($\text{RMSD} = 0.611\text{ cm}^{-1}$ for $E \le 2,000\text{ cm}^{-1}$ and $0.943\text{ cm}^{-1}$ for $E \le 15,000\text{ cm}^{-1}$).
  - `C2H2_sing3D_D2h_CBS_PES_J0.inp` / `.out`: $J=0$ ($K \le 2$, $P_{\max}=12$) calculation ($\text{ZPE} = 5759.948073\text{ cm}^{-1}$).
  - `C2H2_sing3D_D2h_CBS_PES_J1.inp` / `.out`: $J=1$ ($K \le 1$, `gamma 1 2 3 4`) rovibrational calculation ($2B_0 = 2.346535\text{ cm}^{-1}$).
  - `C2H2_sing3D_D2h_CBS_PES_J2.inp` / `.out`: $J=2$ ($K \le 2$, `gamma 1 2 3 4`) rovibrational calculation ($6B_0 - 24D_0 = 7.039563\text{ cm}^{-1}$).
- **6D Benchmark DMS & Transition Intensities ($P_{\max}=12$, $J = 0, 1, 2$)**:
  - `external.chk`: 6D $\mathcal{D}_{\infty h}$-invariant analytical benchmark bisector-frame DMS ($\mu_x, \mu_y, \mu_z$).
  - `C2H2_sing3D_D2h_CBS_DMS_J0.inp` / `.out`: Computes contracted dipole matrix elements `extmatelem{1,2,3}.chk` in `singular_3D`.
  - `C2H2_sing3D_D2h_CBS_intens.inp` / `.out`: Computes 296 K absorption intensities (`J 0 2`) across parallel ($\Sigma_u^+ \leftarrow \Sigma_g^+$), perpendicular ($\Pi_u \leftarrow \Sigma_g^+$), and hot ($\Pi_g \to \Sigma_u^\pm, \Delta_u$) bands.
- `kinetic.chk`: Pre-computed 4th-order `TROVE_EXT_T_S` kinetic energy expansion for $\text{C}_2\text{H}_2$.
- `run.sh`: Helper script to run the full $J=0, 1, 2$ energy, DMS contraction, and intensity pipeline.
