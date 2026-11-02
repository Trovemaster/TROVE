# C2H2 6D `singular_3D` Rovibrational Examples in TROVE

This directory contains input files, pre-computed checkpoint files, and reference output files for 6D rovibrational calculations of acetylene ($\text{C}_2\text{H}_2$) in `TROVE` using:
- **Exact polynomial $t$-$s$ kinetic energy operator** (`TROVE_EXT_T_S`, 4th-order expansion) with the singular $1/\sin^2\rho_i$ ($i=1,2$) terms marked by trailing underscores (`_`) in `kinetic.chk`.
- **Uncontracted 3D coupled singular bend–bend–torsion basis (`singular_3D`)** on $(s_4=\rho_1, s_5=\rho_2, s_6=\tau)$ constructed from associated Legendre basis functions (`sinrho-legendre-k-sing` for modes 4 and 5) and unoptimised half-integer Fourier torsional functions (`fourie-unoptimised` for mode 6 on $\tau \in [0, 4\pi]$), enforcing the exact angular-momentum constraint $K = |k_1 + s\cdot k_2|$ ($k_{\text{tor}} = K/2$) at every basis function and matrix element.
- **Extended $\mathcal{D}_{2h}(\text{M})$ molecular symmetry** on $\tau \in [0, 4\pi]$ (`SYMGROUP D2H(M)`), which separates the physical integer-$K$ and half-integer-$K_{\text{tor}}$ representations into ungerade/gerade blocks ($A_g, A_u, B_{1g}, B_{1u}$ for even $K=0,2,\dots$ and $B_{3g}, B_{2u}$ for odd $K=1,3,\dots$; $B_{2g}$ and $B_{3u}$ are unphysical for integer $K$).
- **$K$-parametrised $J=0$ contracted basis (`KMAX`)** that builds and stores all $K \le K_{\max}$ contracted 3D and 6D basis blocks during the $J=0$ step (`bset_contr.chk`, `matelem*.chk`) for subsequent $J>0$ Coriolis-coupled rovibrational diagonalizations.
- **Linear-geometry-resolved 6D Potential Energy Function (`POTEN_C2H2_MORSE_COSALPHA_COSNTAU-SIN-SIN`)**:
  $$V(\mathbf{q}) = \sum_{ijklmn} c_{ijklmn} f_1^i f_2^j f_3^k f_4^l f_5^m \cos(n\tau) \times \begin{cases} \sin\rho_1 \sin\rho_2 & n > 0, \\ 1 & n = 0, \end{cases}$$
  where $f_1 = 1 - e^{-a_1(R_{\text{CC}}-R_{\text{CC}}^{\text{eq}})}$, $f_{2,3} = 1 - e^{-a_2(r_{1,2}-r_{\text{CH}}^{\text{eq}})}$, and $f_{4,5} = 1 - \cos\rho_{1,2}$. The factor $\sin\rho_1\sin\rho_2$ for $n>0$ guarantees single-valuedness and exact resolution of torsional degeneracies at the linear geometry ($\rho_1=0$ or $\rho_2=0$).

## Files

- `C2H2_sing3D_D2h_J0_kmax2.inp` / `.out`: $J=0$ ($K \le 2$, $P_{\max}=6$) calculation on the test Morse-$\cos\alpha$ PES (`potential.chk`), saving `bset_contr.chk` and `matelem*.chk` for $J>0$.
- `C2H2_sing3D_D2h_J1.inp` / `.out`: $J=1$ ($K \le 1$) rovibrational calculation reading the $K$-parametrised contracted matrix elements and reproducing $E(J=1, K=0) = 2B_0 = 2.3463\text{ cm}^{-1}$ and $\Pi_g(\nu_4) / \Pi_u(\nu_5)$ $\ell$-type doubling.
- `C2H2_sing3D_D2h_CBS_PES_J0.inp` / `.out`: $J=0$ ($K \le 2$, $P_{\max}=6$) calculation on the new 6D `CCSD(T)/CBS(aug-cc-pVTZ, aug-cc-pVQZ)` ab initio PES (`potential_CBS.chk` and `C2H2_CBS_PES_poten.inp_block`) fitted to 4,046 Borg DFA geometries with `POTEN_C2H2_MORSE_COSALPHA_COSNTAU-SIN-SIN`.
- `kinetic.chk`: Pre-computed 4th-order `TROVE_EXT_T_S` kinetic energy expansion for $\text{C}_2\text{H}_2$.
- `potential.chk`: Pre-computed potential expansion (`with_modes`) for the test PES.
- `potential_CBS.chk`: Pre-computed potential expansion (`with_modes`) for the new 6D `CCSD(T)/CBS` PES.
- `run.sh`: Helper script to run the $J=0$ and $J=1$ examples.
