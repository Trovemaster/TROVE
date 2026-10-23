#!/usr/bin/env python3
"""Vectorized exact 1D(r1)*1D(r2)*1D(r3)*2D(s4,s5)*1D(beta) KEO and 2D IHO quadrature."""

import argparse
from collections import defaultdict
import math
import os
import time
import numpy as np
from scipy.special import roots_legendre


BETA_MODE_MAP = {
    (0, 0): 0,   # 1
    (-2, 0): 1,  # csc^2(beta)
    (-1, 1): 2,  # cot(beta)
    (1, 1): 3,   # sin(beta)*cos(beta)
    (2, 0): 4,   # sin^2(beta)
}


class VecJet2:
  """Vectorized 2nd-order 2D Taylor jet in (s4, s5)."""

  __slots__ = ("v", "d4", "d5", "d44", "d45", "d55")

  def __init__(self, v=0.0, d4=0.0, d5=0.0, d44=0.0, d45=0.0, d55=0.0):
    self.v = v
    self.d4 = d4
    self.d5 = d5
    self.d44 = d44
    self.d45 = d45
    self.d55 = d55

  def __neg__(self):
    return VecJet2(-self.v, -self.d4, -self.d5, -self.d44, -self.d45, -self.d55)

  def __add__(self, o):
    if isinstance(o, (int, float, np.ndarray)):
      return VecJet2(
          self.v + o, self.d4, self.d5, self.d44, self.d45, self.d55
      )
    return VecJet2(
        self.v + o.v,
        self.d4 + o.d4,
        self.d5 + o.d5,
        self.d44 + o.d44,
        self.d45 + o.d45,
        self.d55 + o.d55,
    )

  def __radd__(self, o):
    return self.__add__(o)

  def __sub__(self, o):
    if isinstance(o, (int, float, np.ndarray)):
      return VecJet2(
          self.v - o, self.d4, self.d5, self.d44, self.d45, self.d55
      )
    return VecJet2(
        self.v - o.v,
        self.d4 - o.d4,
        self.d5 - o.d5,
        self.d44 - o.d44,
        self.d45 - o.d45,
        self.d55 - o.d55,
    )

  def __rsub__(self, o):
    return (-self).__add__(o)

  def __mul__(self, o):
    if isinstance(o, (int, float, np.ndarray)):
      return VecJet2(
          self.v * o,
          self.d4 * o,
          self.d5 * o,
          self.d44 * o,
          self.d45 * o,
          self.d55 * o,
      )
    return VecJet2(
        self.v * o.v,
        self.d4 * o.v + self.v * o.d4,
        self.d5 * o.v + self.v * o.d5,
        self.d44 * o.v + 2.0 * self.d4 * o.d4 + self.v * o.d44,
        self.d45 * o.v + self.d4 * o.d5 + self.d5 * o.d4 + self.v * o.d45,
        self.d55 * o.v + 2.0 * self.d5 * o.d5 + self.v * o.d55,
    )

  def __rmul__(self, o):
    return self.__mul__(o)

  def inv(self):
    iv = 1.0 / self.v
    iv2 = iv * iv
    iv3 = iv2 * iv
    return VecJet2(
        iv,
        -self.d4 * iv2,
        -self.d5 * iv2,
        2.0 * self.d4 * self.d4 * iv3 - self.d44 * iv2,
        2.0 * self.d4 * self.d5 * iv3 - self.d45 * iv2,
        2.0 * self.d5 * self.d5 * iv3 - self.d55 * iv2,
    )


def vjet_cos_sin(phi_j):
  c = np.cos(phi_j.v)
  s = np.sin(phi_j.v)
  d4, d5 = phi_j.d4, phi_j.d5
  return (
      VecJet2(c, -s * d4, -s * d5, -c * d4 * d4, -c * d4 * d5, -c * d5 * d5),
      VecJet2(s, c * d4, c * d5, -s * d4 * d4, -s * d4 * d5, -s * d5 * d5),
  )


ZERO_VJ = VecJet2(0.0)
ONE_VJ = VecJet2(1.0)


def vjpoly_add(pA, pB, scaleB=1.0):
  res = dict(pA)
  for k, v in pB.items():
    nv = res.get(k, ZERO_VJ) + (v * scaleB)
    if (
        np.max(np.abs(nv.v)) < 1e-16
        and np.max(np.abs(nv.d4)) < 1e-16
        and np.max(np.abs(nv.d5)) < 1e-16
        and np.max(np.abs(nv.d44)) < 1e-16
        and np.max(np.abs(nv.d45)) < 1e-16
        and np.max(np.abs(nv.d55)) < 1e-16
    ):
      res.pop(k, None)
    else:
      res[k] = nv
  return res


def vjpoly_mul(pA, pB):
  if not pA or not pB:
    return {}
  res = defaultdict(lambda: ZERO_VJ)
  for (p1a, p2a, p3a, ma, ca), va in pA.items():
    for (p1b, p2b, p3b, mb, cb), vb in pB.items():
      prod = va * vb
      p1, p2, p3, m = p1a + p1b, p2a + p2b, p3a + p3b, ma + mb
      if ca + cb < 2:
        res[(p1, p2, p3, m, ca + cb)] = res[(p1, p2, p3, m, ca + cb)] + prod
      else:
        res[(p1, p2, p3, m, 0)] = res[(p1, p2, p3, m, 0)] + prod
        res[(p1, p2, p3, m + 2, 0)] = res[(p1, p2, p3, m + 2, 0)] - prod
  return vjpoly_add({}, res)


def vjpoly_diff_q(pA, q_idx):
  res = defaultdict(lambda: ZERO_VJ)
  for (p1, p2, p3, m, c), v in pA.items():
    if q_idx == 0 and p1 > 0:
      res[(p1 + 1, p2, p3, m, c)] = res[(p1 + 1, p2, p3, m, c)] - v * p1
    elif q_idx == 1 and p2 > 0:
      res[(p1, p2 + 1, p3, m, c)] = res[(p1, p2 + 1, p3, m, c)] - v * p2
    elif q_idx == 2 and p3 > 0:
      res[(p1, p2, p3 + 1, m, c)] = res[(p1, p2, p3 + 1, m, c)] - v * p3
    elif q_idx == 3:
      res[(p1, p2, p3, m, c)] = res[(p1, p2, p3, m, c)] + VecJet2(
          v.d4, v.d44, v.d45, 0.0, 0.0, 0.0
      )
    elif q_idx == 4:
      res[(p1, p2, p3, m, c)] = res[(p1, p2, p3, m, c)] + VecJet2(
          v.d5, v.d45, v.d55, 0.0, 0.0, 0.0
      )
    elif q_idx == 5:
      if m != 0:
        if c == 0:
          res[(p1, p2, p3, m - 1, 1)] = res[(p1, p2, p3, m - 1, 1)] + v * m
        else:
          res[(p1, p2, p3, m - 1, 0)] = res[(p1, p2, p3, m - 1, 0)] + v * m
          res[(p1, p2, p3, m + 1, 0)] = res[(p1, p2, p3, m + 1, 0)] - v * m
      if c == 1:
        res[(p1, p2, p3, m + 1, 0)] = res[(p1, p2, p3, m + 1, 0)] - v
  return vjpoly_add({}, res)


def eval_exact_4d_keo_vec(s4_arr, s5_arr):
  """Evaluates exact 4D-decomposed G_vib, G_rot, G_cor, U on 1D arrays s4_arr, s5_arr."""
  z = np.zeros_like(s4_arr)
  o = np.ones_like(s4_arr)
  s4_j = VecJet2(s4_arr, o, z, z, z, z)
  s5_j = VecJet2(s5_arr, z, o, z, z, z)
  phi = [
      s5_j * (math.sqrt(2.0) / 3.0),
      2.0 * math.pi / 3.0
      - s4_j * (1.0 / math.sqrt(6.0))
      - s5_j * (1.0 / (3.0 * math.sqrt(2.0))),
      4.0 * math.pi / 3.0
      + s4_j * (1.0 / math.sqrt(6.0))
      - s5_j * (1.0 / (3.0 * math.sqrt(2.0))),
  ]
  cs = [vjet_cos_sin(phi[i]) for i in range(3)]
  c = [cs[i][0] for i in range(3)]
  s = [cs[i][1] for i in range(3)]
  Delta = (
      (s[2] * c[1] - c[2] * s[1])
      + (s[0] * c[2] - c[0] * s[2])
      + (s[1] * c[0] - c[1] * s[0])
  )
  inv_Delta = Delta.inv()
  M_inv = [[ZERO_VJ for _ in range(3)] for _ in range(3)]
  for i in range(3):
    j, k = (i + 1) % 3, (i + 2) % 3
    M_inv[0][i] = (c[j] - c[k]) * inv_Delta
    M_inv[1][i] = (s[j] - s[k]) * inv_Delta
    M_inv[2][i] = (s[k] * c[j] - c[k] * s[j]) * inv_Delta
  w_rot_z = [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]
  w_s4 = [0.0, -math.sqrt(6.0) / 2.0, math.sqrt(6.0) / 2.0]
  w_s5 = [math.sqrt(2.0), -math.sqrt(2.0) / 2.0, -math.sqrt(2.0) / 2.0]
  sb_p = {(0, 0, 0, 1, 0): VecJet2(o)}
  cb_p = {(0, 0, 0, 0, 1): VecJet2(o)}
  cscb_p = {(0, 0, 0, -1, 0): VecJet2(o)}
  cotb_p = {(0, 0, 0, -1, 1): VecJet2(o)}
  r_hat, b_hat, p_hat = [], [], []
  for i in range(3):
    r_hat.append(
        [vjpoly_add({}, sb_p, c[i]), vjpoly_add({}, sb_p, s[i]), dict(cb_p)]
    )
    b_hat.append([
        vjpoly_add({}, cb_p, c[i]),
        vjpoly_add({}, cb_p, s[i]),
        vjpoly_add({}, sb_p, -1.0),
    ])
    p_hat.append([
        vjpoly_add({}, {(0, 0, 0, 0, 0): VecJet2(o)}, -s[i]),
        vjpoly_add({}, {(0, 0, 0, 0, 0): VecJet2(o)}, c[i]),
        {},
    ])
  S_hat = [[[{} for _ in range(3)] for _ in range(3)] for _ in range(9)]
  for i in range(3):
    for x in range(3):
      S_hat[i][i][x] = r_hat[i][x]
      S_hat[6][i][x] = vjpoly_add({}, b_hat[i][x], M_inv[0][i])
      S_hat[7][i][x] = vjpoly_add({}, b_hat[i][x], M_inv[1][i])
      S_hat[5][i][x] = vjpoly_add({}, b_hat[i][x], M_inv[2][i])
  for l_row, w_vec in ((3, w_s4), (4, w_s5), (8, w_rot_z)):
    wc = c[0] * w_vec[0] + c[1] * w_vec[1] + c[2] * w_vec[2]
    ws = s[0] * w_vec[0] + s[1] * w_vec[1] + s[2] * w_vec[2]
    for i in range(3):
      for x in range(3):
        term1 = vjpoly_mul(vjpoly_add({}, cscb_p, w_vec[i]), p_hat[i][x])
        comb_rot = vjpoly_add(
            vjpoly_add({}, S_hat[6][i][x], wc), S_hat[7][i][x], ws
        )
        S_hat[l_row][i][x] = vjpoly_add(term1, vjpoly_mul(cotb_p, comb_rot))
  ir_polys = [
      {(1, 0, 0, 0, 0): VecJet2(o)},
      {(0, 1, 0, 0, 0): VecJet2(o)},
      {(0, 0, 1, 0, 0): VecJet2(o)},
  ]
  st = [[[{} for _ in range(3)] for _ in range(4)] for _ in range(9)]
  for l in range(9):
    for i in range(3):
      for x in range(3):
        st[l][i + 1][x] = (
            vjpoly_mul(S_hat[l][i][x], ir_polys[i])
            if l >= 3
            else S_hat[l][i][x]
        )
    for x in range(3):
      st[l][0][x] = vjpoly_add(
          {}, vjpoly_add(vjpoly_add(st[l][1][x], st[l][2][x]), st[l][3][x]), -1.0
      )
  g_vib = [[[{} for _ in range(2)] for _ in range(6)] for _ in range(6)]
  g_rot = [[[{} for _ in range(2)] for _ in range(3)] for _ in range(3)]
  g_cor = [[[{} for _ in range(2)] for _ in range(3)] for _ in range(6)]
  for l1 in range(9):
    for l2 in range(l1, 9):
      sum_X, sum_Y = {}, {}
      for x in range(3):
        sum_X = vjpoly_add(sum_X, vjpoly_mul(st[l1][0][x], st[l2][0][x]))
      for n in (1, 2, 3):
        for x in range(3):
          sum_Y = vjpoly_add(sum_Y, vjpoly_mul(st[l1][n][x], st[l2][n][x]))
      for idx, val in enumerate((sum_X, sum_Y)):
        val_flt = {k: v.v for k, v in val.items() if np.max(np.abs(v.v)) >= 1e-16}
        if l1 < 6 and l2 < 6:
          g_vib[l1][l2][idx] = val_flt
          g_vib[l2][l1][idx] = val_flt
        elif l1 >= 6 and l2 >= 6:
          g_rot[l1 - 6][l2 - 6][idx] = val_flt
          g_rot[l2 - 6][l1 - 6][idx] = val_flt
        elif l1 < 6 and l2 >= 6:
          g_cor[l1][l2 - 6][idx] = val_flt

  def epsil(a, b, c_idx):
    if (a, b, c_idx) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
      return 1
    if (a, b, c_idx) in ((2, 1, 0), (0, 2, 1), (1, 0, 2)):
      return -1
    return 0

  pseudo = [{}, {}]
  for n in range(4):
    idx = 0 if n == 0 else 1
    u1_n, u2_n, u3_n, u4_n = {}, {}, {}, {}
    for x0 in range(3):
      curl_x0 = {}
      for q1 in range(3):
        for x1 in range(3):
          e = epsil(x0, q1, x1)
          if e != 0:
            curl_x0 = vjpoly_add(curl_x0, st[6 + q1][n][x1], float(e))
            for q2 in range(6):
              u2_n = vjpoly_add(
                  u2_n,
                  vjpoly_mul(
                      vjpoly_diff_q(st[6 + q1][n][x0], q2), st[q2][n][x1]
                  ),
                  0.25 * e,
              )
      u1_n = vjpoly_add(u1_n, vjpoly_mul(curl_x0, curl_x0), 0.125)
      div_x0 = {}
      for q2 in range(6):
        ds_q2 = vjpoly_diff_q(st[q2][n][x0], q2)
        div_x0 = vjpoly_add(div_x0, ds_q2)
        for q1 in range(6):
          u3_n = vjpoly_add(
              u3_n,
              vjpoly_mul(st[q1][n][x0], vjpoly_diff_q(ds_q2, q1)),
              -0.25,
          )
      u4_n = vjpoly_add(u4_n, vjpoly_mul(div_x0, div_x0), -0.125)
    tot_u = vjpoly_add(vjpoly_add(u1_n, u2_n), vjpoly_add(u3_n, u4_n))
    for k, v in tot_u.items():
      if np.max(np.abs(v.v)) >= 1e-16:
        pseudo[idx][k] = pseudo[idx].get(k, 0.0) + v.v
  return g_vib, g_rot, g_cor, pseudo


def build_polar_grid(R_max=2.45, n_rho=80, n_theta=96):
  """Constructs C3v(M)-invariant 2D polar Gauss-Legendre x Fourier quadrature grid."""
  xi_leg, w_leg = roots_legendre(n_rho)
  rho_pts = 0.5 * R_max * (xi_leg + 1.0)
  w_rho = 0.5 * R_max * w_leg * rho_pts
  theta_pts = np.linspace(0.0, 2.0 * math.pi, n_theta, endpoint=False)
  w_theta = (2.0 * math.pi / n_theta) * np.ones(n_theta)
  rho_2d, theta_2d = np.meshgrid(rho_pts, theta_pts, indexing="ij")
  s4_grid = (rho_2d * np.cos(theta_2d)).ravel()
  s5_grid = (rho_2d * np.sin(theta_2d)).ravel()
  w_grid = (w_rho[:, None] * w_theta[None, :]).ravel()
  rho_grid = rho_2d.ravel()
  return s4_grid, s5_grid, w_grid, rho_grid


def eval_1d_ho_psi_dpsi(s_arr, vmax, ft):
  """Evaluates 1D harmonic oscillator psi_v(s) and dpsi_v(s)/ds matching me_bnd.f90."""
  n_pts = len(s_arr)
  psi = np.zeros((vmax + 1, n_pts), dtype=np.float64)
  dpsi = np.zeros((vmax + 1, n_pts), dtype=np.float64)
  xi = s_arr / ft
  pi4 = 0.75112554446494248285870300477622 / math.sqrt(ft)
  p1 = np.full(n_pts, pi4, dtype=np.float64)
  p2 = np.zeros(n_pts, dtype=np.float64)
  psi[0] = p1
  dpsi[0] = 0.0
  for n in range(1, vmax + 1):
    v = float(n)
    vm1 = float(n - 1)
    p3 = p2
    p2 = p1
    p1 = xi * math.sqrt(2.0 / v) * p2 - math.sqrt(vm1 / v) * p3
    psi[n] = p1
    dpsi[n] = p2 * math.sqrt(2.0 * v) / ft
  exp_fac = np.exp(-0.5 * xi * xi)
  psi = exp_fac[None, :] * psi
  dpsi = exp_fac[None, :] * dpsi - (s_arr[None, :] / (ft * ft)) * psi
  return psi, dpsi


def build_1d_op_densities(s_arr, vmax, ft):
  """Builds D[iop, vi, vj, p] for iop in 0..3 matching ME_Harmonic_Numeric."""
  psi, dpsi = eval_1d_ho_psi_dpsi(s_arr, vmax, ft)
  n_pts = len(s_arr)
  nv = vmax + 1
  D = np.zeros((4, nv, nv, n_pts), dtype=np.float64)
  # iop=0: psi_vi * psi_vj (matelements(-1, ...))
  D[0] = psi[:, None, :] * psi[None, :, :]
  # iop=1: psi_vi * dpsi_vj (matelements(1, lambda, vl, vr))
  D[1] = psi[:, None, :] * dpsi[None, :, :]
  # iop=2: -dpsi_vi * psi_vj (-matelements(1, lambda, vr, vl))
  D[2] = -dpsi[:, None, :] * psi[None, :, :]
  # iop=3: -dpsi_vi * dpsi_vj (matelements(2, ...))
  D[3] = -dpsi[:, None, :] * dpsi[None, :, :]
  return D


def compute_me2d_s4s5(phi_dict, kmax, vmax, ft, s4_grid, s5_grid, w_grid):
  """Computes me2d_s4s5(0:3, 0:3, 0:kmax, 0:kmax, 0:vmax, 0:vmax, 0:vmax, 0:vmax)."""
  nv = vmax + 1
  nk = kmax + 1
  n_pts = len(w_grid)
  D4 = build_1d_op_densities(s4_grid, vmax, ft).reshape(4, nv * nv, n_pts)
  D5 = build_1d_op_densities(s5_grid, vmax, ft).reshape(4, nv * nv, n_pts)

  # Collect active (i4, i5) pairs and stack their weighted grid values
  pairs = sorted(phi_dict.keys())
  phi_w_mat = np.array([phi_dict[pair] * w_grid for pair in pairs], dtype=np.float64)

  me2d = np.zeros((4, 4, nk, nk, nv, nv, nv, nv), dtype=np.float64)
  active_ops = [
      (0, 0),
      (1, 0),
      (2, 0),
      (0, 1),
      (0, 2),
      (3, 0),
      (0, 3),
      (2, 1),
      (1, 2),
  ]
  for iop4, iop5 in active_ops:
    # Contract over quadrature points p:
    # For each pair m, integral[m, a4, a5] = sum_p (phi_w_mat[m, p] * D4[iop4, a4, p]) * D5[iop5, a5, p]
    d4_op = D4[iop4]  # (nv*nv, n_pts)
    d5_op_T = D5[iop5].T  # (n_pts, nv*nv)
    for m_idx, (i4, i5) in enumerate(pairs):
      pw = phi_w_mat[m_idx]  # (n_pts,)
      block = (d4_op * pw[None, :]) @ d5_op_T  # (nv*nv, nv*nv)
      me2d[iop4, iop5, i4, i5] = block.reshape(nv, nv, nv, nv)
  return me2d


def save_keo_2d_s4s5_chk(filepath, me2d, kmax, vmax, binary=False, info=""):
  """Writes keo_2d_s4s5.chk with the table me2d(iop4, iop5, i4, i5, v4, v4', v5, v5').

  ASCII (default): the keyword line 'keo_2d_s4s5 ...', a line with kmax, vmax and the number of
  entries N, a line describing how the table was made, a column-header line, N lines
  'iop4 iop5 i4 i5 v4 v4' v5 v5' value' for the non-zero elements (17 significant digits, so
  the binary values are reproduced exactly) and the line 'End of keo_2d_s4s5'.
  Binary (binary=True): Fortran stream file with int32 kmax, vmax followed by the float64 array in
  Fortran order.
  """
  if binary:
    header = np.array([kmax, vmax], dtype=np.int32)
    data_f = np.asfortranarray(me2d, dtype=np.float64)
    with open(filepath, "wb") as f:
      f.write(header.tobytes())
      f.write(data_f.tobytes(order="F"))
    return
  idx = np.nonzero(me2d)
  vals = me2d[idx]
  table = np.column_stack([np.asarray(i, dtype=np.float64) for i in idx] + [vals])
  with open(filepath, "w", encoding="utf-8") as f:
    f.write("keo_2d_s4s5   2D (s4,s5) integrals of the KEO terms of modes 4 and 5\n")
    f.write(f"{kmax:6d}{vmax:6d}{len(vals):12d}   <- kmax, vmax, number of entries\n")
    f.write(f" {info}\n")
    f.write(" iop4 iop5  i4  i5  v4 v4'  v5 v5'                   value\n")
    np.savetxt(f, table, fmt="%5d%5d%4d%4d%4d%4d%4d%4d %24.16e")
    f.write("End of keo_2d_s4s5\n")


def get_i4_i5_pairs(rank=63, kmax=8):
  """Returns `rank` non-zero (i4, i5) pairs in {0..kmax}x{0..kmax} \\ {(0,0)}, ordered by total degree."""
  pairs = []
  for deg in range(1, 2 * kmax + 1):
    for i4 in range(min(deg, kmax), max(0, deg - kmax) - 1, -1):
      i5 = deg - i4
      pairs.append((i4, i5))
      if len(pairs) == rank:
        return pairs
  raise ValueError(f"Cannot fit rank={rank} in kmax={kmax}")


def generate_exact_2d_keo(out_dir, vmax=6, kmax=8, rank=63, ft=0.2865833376039929, R_max=2.45, n_rho=80, n_theta=96,
                          binary=False):
  """Generates kinetic.chk and keo_2d_s4s5.chk for the exact 2D (s4,s5) KEO."""
  t0 = time.time()
  s4_grid, s5_grid, w_grid, rho_grid = build_polar_grid(R_max, n_rho, n_theta)
  s4_all = np.concatenate([[0.0], s4_grid])
  s5_all = np.concatenate([[0.0], s5_grid])

  gv, gr, gc, u = eval_exact_4d_keo_vec(s4_all, s5_all)
  print(f"[{time.time() - t0:.2f}s] Evaluated exact 4D+2D KEO on {len(s4_all)} points.")

  # Collect all zero-at-equilibrium 2D functions across gvib, grot, gcor, pseudo (mass-split idx=0,1)
  all_records = {}
  eq_vals = {}
  for l1 in range(6):
    for l2 in range(6):
      for idx in (0, 1):
        for k, val_all in gv[l1][l2][idx].items():
          key = ("gvib", l1, l2, idx, k)
          eq_vals[key] = float(val_all[0])
          all_records[key] = val_all[1:] - val_all[0]
  for l1 in range(3):
    for l2 in range(3):
      for idx in (0, 1):
        for k, val_all in gr[l1][l2][idx].items():
          key = ("grot", l1, l2, idx, k)
          eq_vals[key] = float(val_all[0])
          all_records[key] = val_all[1:] - val_all[0]
  for l1 in range(6):
    for l2 in range(3):
      for idx in (0, 1):
        for k, val_all in gc[l1][l2][idx].items():
          key = ("gcor", l1, l2, idx, k)
          eq_vals[key] = float(val_all[0])
          all_records[key] = val_all[1:] - val_all[0]
  for idx in (0, 1):
    for k, val_all in u[idx].items():
      key = ("u", 0, 0, idx, k)
      eq_vals[key] = float(val_all[0])
      all_records[key] = val_all[1:] - val_all[0]

  env_w = np.exp(-0.25 * (rho_grid / ft) ** 2) * np.sqrt(w_grid)
  keys_nz = [
      k for k in sorted(all_records.keys())
      if np.max(np.abs(all_records[k] * env_w)) > 1e-15
  ]
  mat_w = np.array([all_records[k] * env_w for k in keys_nz])
  u_svd, s_svd, vt_svd = np.linalg.svd(mat_w, full_matrices=False)
  max_w_err = np.max(np.abs(mat_w - (u_svd[:, :rank] * s_svd[:rank]) @ vt_svd[:rank, :]))
  print(
      f"[{time.time() - t0:.2f}s] SVD rank={rank}: s[{rank-1}]={s_svd[rank-1]:.3e}, s[{rank}]={s_svd[rank]:.3e}, "
      f"max Gaussian-weighted error={max_w_err:.3e}"
  )

  basis_r = vt_svd[:rank, :] / env_w[None, :]
  coeffs_r = u_svd[:, :rank] * s_svd[:rank]
  coeffs_map = {keys_nz[i]: coeffs_r[i] for i in range(len(keys_nz))}

  pairs_r = get_i4_i5_pairs(rank=rank, kmax=kmax)
  phi_dict = {(0, 0): np.ones_like(s4_grid)}
  for m_idx, pair in enumerate(pairs_r):
    phi_dict[pair] = basis_r[m_idx]

  # Compute 2D integrals me2d_s4s5 and save keo_2d_s4s5.chk
  me2d = compute_me2d_s4s5(phi_dict, kmax, vmax, ft, s4_grid, s5_grid, w_grid)
  keo2d_path = os.path.join(out_dir, "keo_2d_s4s5.chk")
  info = f"mode=exact ft={ft!r} rmax={R_max!r} n_rho={n_rho} n_theta={n_theta} rank={rank}"
  save_keo_2d_s4s5_chk(keo2d_path, me2d, kmax, vmax, binary=binary, info=info)
  print(f"[{time.time() - t0:.2f}s] Wrote {keo2d_path} (kmax={kmax}, vmax={vmax}, rank={rank})")

  # Write kinetic.chk with exact equilibrium values at (0,0) and `rank` SVD coefficients at pairs_r
  lines = []
  n_written = 0

  def emit_section(sec_tag, k1, k2, keys_list):
    nonlocal n_written
    icount = 0
    for idx in (0, 1):
      i0 = idx + 1
      for k in sorted(keys_list[idx]):
        p1, p2, p3, m_sin, p_cb = k
        i6 = BETA_MODE_MAP[(m_sin, p_cb)]
        rec_key = (sec_tag, k1 - 1, k2 - 1, idx, k)
        # 1. Equilibrium term at (i4, i5) = (0, 0)
        v0 = eq_vals.get(rec_key, 0.0)
        if abs(v0) > 1e-15:
          icount += 1
          n_written += 1
          lines.append(
              f"{k1:5d}{k2:5d}{icount:8d}  0 {v0:23.14e}"
              f" {i0:3d}{p1:3d}{p2:3d}{p3:3d}  0  0{i6:3d}"
          )
        # 2. Non-constant 2D basis terms at pairs_r
        if rec_key in coeffs_map:
          c_vec = coeffs_map[rec_key]
          for m_idx, (i4, i5) in enumerate(pairs_r):
            cm = float(c_vec[m_idx])
            if abs(cm) > 1e-14:
              icount += 1
              n_written += 1
              lines.append(
                  f"{k1:5d}{k2:5d}{icount:8d}  0 {cm:23.14e}"
                  f" {i0:3d}{p1:3d}{p2:3d}{p3:3d}{i4:3d}{i5:3d}{i6:3d}"
              )

  sentinel = "   987654321  0    0 0          0.00000000   0  0  0  0  0  0  0"
  lines.append("  1  1  1   <-  g_vib   ")
  for i in range(6):
    for j in range(6):
      emit_section("gvib", i + 1, j + 1, [gv[i][j][0].keys(), gv[i][j][1].keys()])
  lines.append(sentinel)

  lines.append("  1  1  1   <-  g_rot   ")
  for i in range(3):
    for j in range(3):
      emit_section("grot", i + 1, j + 1, [gr[i][j][0].keys(), gr[i][j][1].keys()])
  lines.append(sentinel)

  lines.append("  1  1  1   <-  g_cor   ")
  for i in range(6):
    for j in range(3):
      emit_section("gcor", i + 1, j + 1, [gc[i][j][0].keys(), gc[i][j][1].keys()])
  lines.append(sentinel)

  lines.append("  1  1  1   <-  pseudo   ")
  emit_section("u", 1, 1, [u[0].keys(), u[1].keys()])
  lines.append(sentinel)

  lines.append(f" {1.0e-6:14.8e} <- sparse threshold used")
  lines.append("End of kinetic")

  kin_path = os.path.join(out_dir, "kinetic.chk")
  with open(kin_path, "w", encoding="utf-8") as f:
    f.write("\n".join(lines) + "\n")
  print(f"[{time.time() - t0:.2f}s] Wrote {kin_path} ({n_written} exact 1D+2D terms)")


def generate_taylor_keo_2d_chk(out_dir, kmax=4, vmax=6, ft=0.2865833376039929, R_max=2.45, n_rho=80, n_theta=96,
                               binary=False):
  """Generates keo_2d_s4s5.chk for Taylor monomials s4^i4 * s5^i5 (0 <= i4+i5 <= kmax)."""
  t0 = time.time()
  s4_grid, s5_grid, w_grid, _ = build_polar_grid(R_max, n_rho, n_theta)
  phi_dict = {}
  for i4 in range(kmax + 1):
    for i5 in range(kmax + 1 - i4):
      phi_dict[(i4, i5)] = (s4_grid**i4) * (s5_grid**i5)
  me2d = compute_me2d_s4s5(phi_dict, kmax, vmax, ft, s4_grid, s5_grid, w_grid)
  keo2d_path = os.path.join(out_dir, "keo_2d_s4s5.chk")
  info = f"mode=taylor ft={ft!r} rmax={R_max!r} n_rho={n_rho} n_theta={n_theta}"
  save_keo_2d_s4s5_chk(keo2d_path, me2d, kmax, vmax, binary=binary, info=info)
  print(f"[{time.time() - t0:.2f}s] Wrote Taylor {keo2d_path} (kmax={kmax}, vmax={vmax})")


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("--mode", choices=["taylor", "exact"], default="exact")
  parser.add_argument("--kmax", type=int, default=4)
  parser.add_argument("--vmax", type=int, default=6)
  parser.add_argument("--ft", type=float, default=0.2865833376039929)
  parser.add_argument("--rmax", type=float, default=2.45)
  parser.add_argument("--out-dir", type=str, default=".")
  parser.add_argument("--binary", action="store_true",
                      help="write keo_2d_s4s5.chk as a binary stream file instead of ASCII")
  args = parser.parse_args()

  os.makedirs(args.out_dir, exist_ok=True)
  if args.mode == "taylor":
    generate_taylor_keo_2d_chk(args.out_dir, kmax=args.kmax, vmax=args.vmax, ft=args.ft, R_max=args.rmax,
                               binary=args.binary)
  else:
    generate_exact_2d_keo(args.out_dir, vmax=args.vmax, ft=args.ft, R_max=args.rmax, binary=args.binary)
