#!/usr/bin/env python3
"""Independent J=0 DVR check of HCN/HNC vibrational levels on a fitted PES.

(Companion of the TROVE HCN/HNC Jacobi example; reads the POTEN_XYZ_JACOBI_MORSE_COS_MEP
parameter block and, optionally, the DIPOLE_XYZ_JACOBI_MEP block used in the TROVE inputs.)

Jacobi coordinates r = r(CN), R = |H - CM(CN)|, gamma (0 = HCN, 180 = HNC),
exact J=0 KEO with radial functions phi = r R Psi:
  H = -h2/(2 mu_r) d2/dr2 - h2/(2 mu_R) d2/dR2
      + h2 (1/(2 mu_r r^2) + 1/(2 mu_R R^2)) j^2 + V(r, R, gamma)
Sinc DVR (Colbert-Miller) in r and R, Gauss-Legendre DVR in cos(gamma).
Lowest eigenvalues by Lanczos (ARPACK) on a matrix-free operator.

Usage:
  dvr_j0_check.py <poten_block.txt> <a1> <a2> [--nr 20 --nR 48 --ng 90
      --rmin 0.98 --rmax 1.45 --Rmin 1.05 --Rmax 2.70 --nev 160 --vcap 6e4]
"""
import argparse
import json
import time

import numpy as np
from numpy.polynomial import legendre as npleg
from scipy.sparse.linalg import LinearOperator, eigsh

H2_2 = 16.857629206  # hbar^2/(2 * 1 amu * 1 A^2) in cm^-1
M_C, M_N, M_H = 12.0, 14.0030740048, 1.00782503207
MU_R = M_C * M_N / (M_C + M_N)
MU_RR = M_H * (M_C + M_N) / (M_H + M_C + M_N)


def load_block(path):
  c1, c2, terms = {}, {}, []
  for line in open(path):
    p = line.split()
    if len(p) < 5:
      continue
    try:
      i, j, m, v = int(p[1]), int(p[2]), int(p[3]), float(p[4])
    except ValueError:
      continue
    if i == -1:
      c1[m] = v
    elif i == -2:
      c2[m] = v
    else:
      terms.append((i, j, m, v))
  return (np.array([c1[k] for k in sorted(c1)]),
          np.array([c2[k] for k in sorted(c2)]), terms)


def pes(blk, a1, a2, r, R, g):
  c1, c2, terms = blk
  rm = np.cos(np.multiply.outer(g, np.arange(len(c1)))) @ c1
  Rm = np.cos(np.multiply.outer(g, np.arange(len(c2)))) @ c2
  u1 = 1.0 - np.exp(-a1 * (r - rm))
  u2 = 1.0 - np.exp(-a2 * (R - Rm))
  pmax1 = max(t[0] for t in terms)
  pmax2 = max(t[1] for t in terms)
  mmax = max(t[2] for t in terms)
  p1 = [np.ones_like(u1)]
  for _ in range(pmax1):
    p1.append(p1[-1] * u1)
  p2 = [np.ones_like(u2)]
  for _ in range(pmax2):
    p2.append(p2[-1] * u2)
  cm = [np.cos(m * g) for m in range(mmax + 1)]
  v = np.zeros_like(r)
  for i, j, m, c in terms:
    v += c * p1[i] * p2[j] * cm[m]
  return v


def load_dms(path):
  """Reads a DIPOLE_XYZ_JACOBI_MEP block: returns (bm1, bm2, c1, c2, X, Z)."""
  xs, zs, c1, c2, bm = [], [], {}, {}, {}
  section = None
  nx = nz = None
  rows = []
  for line in open(path):
    p = line.split()
    if not p:
      continue
    if p[0].upper() == "NPARAM":
      nx, nz = int(p[1]), int(p[3])
      continue
    if len(p) == 5:
      try:
        rows.append((p[0], int(p[1]), int(p[2]), int(p[3]), float(p[4])))
      except ValueError:
        pass
  for idx, (name, i, j, k, v) in enumerate(rows):
    section = "x" if idx < nx else "z"
    if i == 93:
      bm[j] = v
    elif i == 91:
      c1[k] = v
    elif i == 92:
      c2[k] = v
    elif section == "x":
      xs.append((i, j, k, v))
    else:
      zs.append((i, j, k, v))
  return (bm[1], bm[2], np.array([c1[m] for m in sorted(c1)]),
          np.array([c2[m] for m in sorted(c2)]), xs, zs)


def dms(d, r, R, g):
  bm1, bm2, c1, c2, xs, zs = d
  rm = np.cos(np.multiply.outer(g, np.arange(len(c1)))) @ c1
  Rm = np.cos(np.multiply.outer(g, np.arange(len(c2)))) @ c2
  u1 = 1.0 - np.exp(-bm1 * (r - rm))
  u2 = 1.0 - np.exp(-bm2 * (R - Rm))
  mx = np.zeros_like(r)
  mz = np.zeros_like(r)
  for i, j, k, c in xs:
    mx += c * u1**i * u2**j * np.sin(k * g)
  for i, j, k, c in zs:
    mz += c * u1**i * u2**j * np.cos(k * g)
  return mx, mz


def sinc_t(n, dx, mu):
  i = np.arange(n)
  d = i[:, None] - i[None, :]
  with np.errstate(divide="ignore"):
    t = 2.0 * (-1.0) ** d / (dx**2 * d.astype(float) ** 2)
  t[np.diag_indices(n)] = np.pi**2 / (3.0 * dx**2)
  return H2_2 / mu * t


def main():
  ap = argparse.ArgumentParser()
  ap.add_argument("block")
  ap.add_argument("a1", type=float)
  ap.add_argument("a2", type=float)
  ap.add_argument("--nr", type=int, default=20)
  ap.add_argument("--nR", type=int, default=48)
  ap.add_argument("--ng", type=int, default=90)
  ap.add_argument("--rmin", type=float, default=0.98)
  ap.add_argument("--rmax", type=float, default=1.45)
  ap.add_argument("--Rmin", type=float, default=1.05)
  ap.add_argument("--Rmax", type=float, default=2.70)
  ap.add_argument("--nev", type=int, default=160)
  ap.add_argument("--vcap", type=float, default=6.0e4)
  ap.add_argument("--json")
  ap.add_argument("--dms", help="DIPOLE_XYZ_JACOBI_MEP block for J=0 TMs")
  a = ap.parse_args()

  t0 = time.time()
  blk = load_block(a.block)
  rg = np.linspace(a.rmin, a.rmax, a.nr)
  Rg = np.linspace(a.Rmin, a.Rmax, a.nR)
  x, w = npleg.leggauss(a.ng)          # x = cos(gamma), ascending
  x, w = x[::-1], w[::-1]               # gamma ascending (0 -> pi)
  gg = np.arccos(x)
  # Normalised Legendre transformation T[k, n] = sqrt(w_k) Ptilde_n(x_k).
  tr = np.zeros((a.ng, a.ng))
  for n in range(a.ng):
    c = np.zeros(n + 1)
    c[n] = 1.0
    tr[:, n] = np.sqrt(w) * npleg.legval(x, c) * np.sqrt((2 * n + 1) / 2.0)
  jj = np.arange(a.ng) * (np.arange(a.ng) + 1.0)
  lmat = tr @ np.diag(jj) @ tr.T
  tr_r = sinc_t(a.nr, rg[1] - rg[0], MU_R)
  tr_R = sinc_t(a.nR, Rg[1] - Rg[0], MU_RR)
  r3, R3, g3 = np.meshgrid(rg, Rg, gg, indexing="ij")
  v = pes(blk, a.a1, a.a2, r3.ravel(), R3.ravel(), g3.ravel()).reshape(
      r3.shape)
  vmin = v.min()
  v = np.minimum(v, a.vcap)
  b = H2_2 * (1.0 / (MU_R * rg[:, None] ** 2) + 1.0 / (MU_RR * Rg[None, :] ** 2))
  shape = (a.nr, a.nR, a.ng)
  n = a.nr * a.nR * a.ng
  print(f"grid {shape} = {n} points; Vmin on grid {vmin:.2f}; "
        f"PES eval {time.time() - t0:.1f}s", flush=True)

  def matvec(psi):
    p = psi.reshape(shape)
    out = np.tensordot(tr_r, p, axes=(1, 0))
    out += np.einsum("ij,ajk->aik", tr_R, p, optimize=True)
    out += b[:, :, None] * (p @ lmat.T)
    out += v * p
    return out.ravel()

  op = LinearOperator((n, n), matvec=matvec, dtype=float)
  t1 = time.time()
  ev, vec = eigsh(op, k=a.nev, which="SA", tol=1e-9, ncv=min(n, 3 * a.nev))
  order = np.argsort(ev)
  ev, vec = ev[order], vec[:, order]
  print(f"eigsh {time.time() - t1:.1f}s", flush=True)
  cosg = np.cos(gg)
  zpe = ev[0]
  rows = []
  for i, e in enumerate(ev):
    p2 = (vec[:, i].reshape(shape)) ** 2
    cg = float(np.sum(p2 * cosg[None, None, :]))
    rr = float(np.sum(p2 * rg[:, None, None]))
    RR = float(np.sum(p2 * Rg[None, :, None]))
    iso = "HCN" if cg > 0 else "HNC"
    rows.append({"i": i, "E": float(e - zpe), "iso": iso, "cosg": cg,
                 "r": rr, "R": RR})
  e_hnc = min(rw["E"] for rw in rows if rw["iso"] == "HNC") if any(
      rw["iso"] == "HNC" for rw in rows) else None
  print(f"ZPE = {zpe - 0.0:.4f} cm-1 (V_min on grid {vmin:.3f})")
  for rw in rows:
    rel = (f"  HNC+{rw['E'] - e_hnc:9.3f}" if rw["iso"] == "HNC" and
           e_hnc is not None else "")
    print(f"{rw['i']:4d} {rw['E']:11.4f} {rw['iso']} <cos g>={rw['cosg']:+.3f}"
          f" <r>={rw['r']:.4f} <R>={rw['R']:.4f}{rel}")
  if a.dms:
    mx, mz = dms(load_dms(a.dms), r3.ravel(), R3.ravel(), g3.ravel())
    i_hnc = min((rw["i"] for rw in rows if rw["iso"] == "HNC"), default=None)
    print("\nJ=0 dipole matrix elements (Debye), |<f|mu|i>| for i = HCN ground"
          " and HNC ground; diagonal <mu_z> classifies the isomer")
    for rw in rows:
      f = vec[:, rw["i"]]
      rw["mu_z_diag"] = float(f @ (mz * f))
      rw["mu_x_diag"] = float(f @ (mx * f))
      rw["tm_z_0"] = float(abs(f @ (mz * vec[:, 0])))
      rw["tm_x_0"] = float(abs(f @ (mx * vec[:, 0])))
      if i_hnc is not None:
        rw["tm_z_hnc"] = float(abs(f @ (mz * vec[:, i_hnc])))
        rw["tm_x_hnc"] = float(abs(f @ (mx * vec[:, i_hnc])))
      print(f"{rw['i']:4d} {rw['E']:11.4f} {rw['iso']} <mu_z>={rw['mu_z_diag']:+.5f}"
            f" <mu_x>={rw['mu_x_diag']:+.5f} | TM(HCN0) z={rw['tm_z_0']:.6f}"
            f" x={rw['tm_x_0']:.6f}" +
            (f" | TM(HNC0) z={rw['tm_z_hnc']:.6f} x={rw['tm_x_hnc']:.6f}"
             if i_hnc is not None else ""))
  if a.json:
    json.dump({"zpe": float(zpe), "grid": shape, "levels": rows,
               "args": vars(a)}, open(a.json, "w"), indent=1)


if __name__ == "__main__":
  main()
