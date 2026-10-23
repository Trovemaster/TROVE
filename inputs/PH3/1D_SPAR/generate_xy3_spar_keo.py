#!/usr/bin/env python3
"""External SPAR KEO Generator for XY3 Pyramidal Molecules (e.g., PH3, NH3).

Generates a TROVE-independent sum-of-products (SPAR) Kinetic Energy Operator
checkpoint (`kinetic.chk`) and basic function definitions (`combination.chk`)
using the t-s matrix formalism in a C3v/D3h symmetry-adapted frame (Trisector
z-axis + Eckart/Mellor symmetric azimuths) with coordinates `R-SYMPHI-BETA`:
  - q_1 = r_1 (X-Y1 bond length, exact analytical in 1/r_1, 1/r_1^2)
  - q_2 = r_2 (X-Y2 bond length, exact analytical in 1/r_2, 1/r_2^2)
  - q_3 = r_3 (X-Y3 bond length, exact analytical in 1/r_3, 1/r_3^2)
  - q_4 = s_4 = (2*delta_1 - delta_2 - delta_3)/sqrt(6) (E_a dihedral, Taylor)
  - q_5 = s_5 = (delta_2 - delta_3)/sqrt(2)             (E_b dihedral, Taylor)
  - q_6 = beta (umbrella angle from trisector z-axis, exact in sin^m(beta)*cos^c(beta))
"""

import argparse
from collections import defaultdict
from fractions import Fraction
import math
import os
import time


# ==============================================================================
# 1. Exact Algebraic Number Field Q(sqrt(2), sqrt(3))
#    Represents x = a + b*sqrt(2) + c*sqrt(3) + d*sqrt(6) with a,b,c,d in Q.
# ==============================================================================
class Q23:
  __slots__ = ("a", "b", "c", "d")

  def __init__(self, a=0, b=0, c=0, d=0):
    self.a = a if isinstance(a, Fraction) else Fraction(a)
    self.b = b if isinstance(b, Fraction) else Fraction(b)
    self.c = c if isinstance(c, Fraction) else Fraction(c)
    self.d = d if isinstance(d, Fraction) else Fraction(d)

  def is_zero(self):
    return self.a == 0 and self.b == 0 and self.c == 0 and self.d == 0

  def __eq__(self, other):
    if isinstance(other, Q23):
      return (
          self.a == other.a
          and self.b == other.b
          and self.c == other.c
          and self.d == other.d
      )
    if isinstance(other, (int, Fraction)):
      return self.a == other and self.b == 0 and self.c == 0 and self.d == 0
    return False

  def __bool__(self):
    return not self.is_zero()

  def __pos__(self):
    return self

  def __neg__(self):
    return Q23(-self.a, -self.b, -self.c, -self.d)

  def __add__(self, other):
    if isinstance(other, (int, Fraction)):
      return Q23(self.a + other, self.b, self.c, self.d)
    return Q23(
        self.a + other.a, self.b + other.b, self.c + other.c, self.d + other.d
    )

  def __radd__(self, other):
    return self.__add__(other)

  def __sub__(self, other):
    if isinstance(other, (int, Fraction)):
      return Q23(self.a - other, self.b, self.c, self.d)
    return Q23(
        self.a - other.a, self.b - other.b, self.c - other.c, self.d - other.d
    )

  def __rsub__(self, other):
    return (-self).__add__(other)

  def __mul__(self, other):
    if isinstance(other, (int, Fraction)):
      f = other if isinstance(other, Fraction) else Fraction(other)
      return Q23(self.a * f, self.b * f, self.c * f, self.d * f)
    # (a1 + b1*s2 + c1*s3 + d1*s6)*(a2 + b2*s2 + c2*s3 + d2*s6)
    a1, b1, c1, d1 = self.a, self.b, self.c, self.d
    a2, b2, c2, d2 = other.a, other.b, other.c, other.d
    return Q23(
        a1 * a2 + 2 * b1 * b2 + 3 * c1 * c2 + 6 * d1 * d2,
        a1 * b2 + b1 * a2 + 3 * c1 * d2 + 3 * d1 * c2,
        a1 * c2 + c1 * a2 + 2 * b1 * d2 + 2 * d1 * b2,
        a1 * d2 + d1 * a2 + b1 * c2 + c1 * b2,
    )

  def __rmul__(self, other):
    return self.__mul__(other)

  def __truediv__(self, other):
    if isinstance(other, (int, Fraction)):
      f = other if isinstance(other, Fraction) else Fraction(other)
      return Q23(self.a / f, self.b / f, self.c / f, self.d / f)
    raise NotImplementedError("Division by general Q23 not needed")

  def to_float(self):
    s2 = math.sqrt(2.0)
    s3 = math.sqrt(3.0)
    s6 = math.sqrt(6.0)
    return (
        float(self.a)
        + float(self.b) * s2
        + float(self.c) * s3
        + float(self.d) * s6
    )

  def __repr__(self):
    parts = []
    for val, sym in (
        (self.a, ""),
        (self.b, "*sqrt(2)"),
        (self.c, "*sqrt(3)"),
        (self.d, "*sqrt(6)"),
    ):
      if val != 0:
        parts.append(f"{val}{sym}")
    return " + ".join(parts) if parts else "0"


ZERO_Q = Q23(0)
ONE_Q = Q23(1)
SQRT2_Q = Q23(0, 1, 0, 0)
SQRT3_Q = Q23(0, 0, 1, 0)
SQRT6_Q = Q23(0, 0, 0, 1)
INV_SQRT2 = Q23(0, Fraction(1, 2), 0, 0)
INV_SQRT6 = Q23(0, 0, 0, Fraction(1, 6))


# ==============================================================================
# 2. Sparse Polynomial-Trig Algebra in (r1^-1, r2^-1, r3^-1, s4, s5, sin(beta), cos(beta))
#    Monomial key: (p1, p2, p3, p4, p5, m_sin, p_cb)
#      - p1, p2, p3 in {0, 1, 2}: powers of 1/r1, 1/r2, 1/r3
#      - p4, p5 >= 0: powers of s4, s5
#      - m_sin in Z: power of sin(beta) (negative means csc(beta)^|m_sin|)
#      - p_cb in {0, 1}: power of cos(beta), reduced via cos^2(beta) = 1 - sin^2(beta)
# ==============================================================================
def poly_clean(p):
  return {k: v for k, v in p.items() if not v.is_zero()}


def poly_add(pA, pB, scaleB=ONE_Q):
  res = dict(pA)
  if isinstance(scaleB, (int, Fraction)):
    scaleB = Q23(scaleB)
  for k, v in pB.items():
    nv = res.get(k, ZERO_Q) + v * scaleB
    if nv.is_zero():
      res.pop(k, None)
    else:
      res[k] = nv
  return res


def poly_scale(pA, c):
  if isinstance(c, (int, Fraction)):
    c = Q23(c)
  if c.is_zero():
    return {}
  return {k: v * c for k, v in pA.items()}


def poly_mul(pA, pB, max_deg=None):
  """Multiplies two sparse polynomials, truncating to p4+p5 <= max_deg and reducing cos^2(beta) -> 1 - sin^2(beta)."""
  if not pA or not pB:
    return {}
  res = defaultdict(lambda: ZERO_Q)
  for (p1a, p2a, p3a, p4a, p5a, ma, ca), va in pA.items():
    for (p1b, p2b, p3b, p4b, p5b, mb, cb), vb in pB.items():
      p4 = p4a + p4b
      p5 = p5a + p5b
      if max_deg is not None and (p4 + p5 > max_deg):
        continue
      prod = va * vb
      if prod.is_zero():
        continue
      p1 = p1a + p1b
      p2 = p2a + p2b
      p3 = p3a + p3b
      m = ma + mb
      if ca + cb < 2:
        key = (p1, p2, p3, p4, p5, m, ca + cb)
        res[key] = res[key] + prod
      else:
        # cos^2(beta) = 1 - sin^2(beta)
        key0 = (p1, p2, p3, p4, p5, m, 0)
        key2 = (p1, p2, p3, p4, p5, m + 2, 0)
        res[key0] = res[key0] + prod
        res[key2] = res[key2] - prod
  return poly_clean(res)


def poly_truncate(pA, max_deg):
  return {k: v for k, v in pA.items() if k[3] + k[4] <= max_deg and not v.is_zero()}


def poly_diff_q(pA, q_idx):
  """Differentiates a sparse polynomial with respect to q_{q_idx} (q_idx in 0..5)."""
  if not pA:
    return {}
  res = defaultdict(lambda: ZERO_Q)
  for (p1, p2, p3, p4, p5, m, c), v in pA.items():
    if q_idx == 0:  # d/dr1 of r1^(-p1) = -p1 * r1^(-(p1+1))
      if p1 > 0:
        key = (p1 + 1, p2, p3, p4, p5, m, c)
        res[key] = res[key] - v * p1
    elif q_idx == 1:  # d/dr2
      if p2 > 0:
        key = (p1, p2 + 1, p3, p4, p5, m, c)
        res[key] = res[key] - v * p2
    elif q_idx == 2:  # d/dr3
      if p3 > 0:
        key = (p1, p2, p3 + 1, p4, p5, m, c)
        res[key] = res[key] - v * p3
    elif q_idx == 3:  # d/ds4
      if p4 > 0:
        key = (p1, p2, p3, p4 - 1, p5, m, c)
        res[key] = res[key] + v * p4
    elif q_idx == 4:  # d/ds5
      if p5 > 0:
        key = (p1, p2, p3, p4, p5 - 1, m, c)
        res[key] = res[key] + v * p5
    elif q_idx == 5:  # d/dbeta of sin^m(b)*cos^c(b)
      # = m * sin^(m-1)(b) * cos^(c+1)(b) - c * sin^(m+1)(b) * cos^(c-1)(b)
      if m != 0:
        if c == 0:
          key = (p1, p2, p3, p4, p5, m - 1, 1)
          res[key] = res[key] + v * m
        else:  # c == 1 -> cos^2(b) = 1 - sin^2(b)
          key0 = (p1, p2, p3, p4, p5, m - 1, 0)
          key2 = (p1, p2, p3, p4, p5, m + 1, 0)
          res[key0] = res[key0] + v * m
          res[key2] = res[key2] - v * m
      if c == 1:
        key_c = (p1, p2, p3, p4, p5, m + 1, 0)
        res[key_c] = res[key_c] - v
  return poly_clean(res)


# ==============================================================================
# 3. Build Unit Vectors u_i(s4, s5, beta) and B_hat, S_hat, G_uv, U
# ==============================================================================
def build_unit_vectors_sparse(max_deg):
  # Dihedral differences opposite to bonds 1, 2, 3:
  #   delta_1 = phi_3 - phi_2 = 2*pi/3 + d3 - d2
  #   delta_2 = phi_1 - phi_3 + 2*pi = 2*pi/3 + d1 - d3
  #   delta_3 = phi_2 - phi_1 = 2*pi/3 + d2 - d1
  # Symmetry-adapted dihedral coordinates (matching TROVE R-SYMPHI-BETA and pot_xy3):
  #   s4 = (2*delta_1 - delta_2 - delta_3)/sqrt(6) = sqrt(3/2)*(d3 - d2)
  #   s5 = (delta_2 - delta_3)/sqrt(2)             = (2*d1 - d2 - d3)/sqrt(2)
  # With trisector/Mellor constraint d1 + d2 + d3 = 0:
  #   d1 = +sqrt(2)/3 * s5
  #   d2 = -1/sqrt(6) * s4 - 1/(3*sqrt(2)) * s5
  #   d3 = +1/sqrt(6) * s4 - 1/(3*sqrt(2)) * s5
  d1 = {(0, 0, 0, 0, 1, 0, 0): Q23(0, Fraction(1, 3), 0, 0)}
  d2 = {
      (0, 0, 0, 1, 0, 0, 0): Q23(0, 0, 0, Fraction(-1, 6)),
      (0, 0, 0, 0, 1, 0, 0): Q23(0, Fraction(-1, 6), 0, 0),
  }
  d3 = {
      (0, 0, 0, 1, 0, 0, 0): Q23(0, 0, 0, Fraction(1, 6)),
      (0, 0, 0, 0, 1, 0, 0): Q23(0, Fraction(-1, 6), 0, 0),
  }

  def cos_sin_series(d_poly):
    # Compute powers d_poly^k up to max_deg
    p_pow = {(0, 0, 0, 0, 0, 0, 0): ONE_Q}
    c_ser = dict(p_pow)
    s_ser = {}
    for k in range(1, max_deg + 1):
      p_pow = poly_mul(p_pow, d_poly, max_deg=max_deg)
      inv_fact = Fraction(1, math.factorial(k))
      if k % 2 == 1:
        sign = Fraction(1 if ((k - 1) // 2) % 2 == 0 else -1)
        s_ser = poly_add(s_ser, p_pow, sign * inv_fact)
      else:
        sign = Fraction(1 if (k // 2) % 2 == 0 else -1)
        c_ser = poly_add(c_ser, p_pow, sign * inv_fact)
    return c_ser, s_ser

  c1, s1 = cos_sin_series(d1)
  c2, s2 = cos_sin_series(d2)
  c3, s3 = cos_sin_series(d3)

  half = Q23(Fraction(-1, 2))
  s3_2 = Q23(0, 0, Fraction(1, 2), 0)

  cp1, sp1 = c1, s1
  cp2 = poly_add(poly_scale(c2, half), poly_scale(s2, -s3_2))
  sp2 = poly_add(poly_scale(c2, s3_2), poly_scale(s2, half))
  cp3 = poly_add(poly_scale(c3, half), poly_scale(s3, s3_2))
  sp3 = poly_add(poly_scale(c3, -s3_2), poly_scale(s3, half))

  sb_poly = {(0, 0, 0, 0, 0, 1, 0): ONE_Q}
  cb_poly = {(0, 0, 0, 0, 0, 0, 1): ONE_Q}

  u1 = [poly_mul(sb_poly, cp1), poly_mul(sb_poly, sp1), dict(cb_poly)]
  u2 = [poly_mul(sb_poly, cp2), poly_mul(sb_poly, sp2), dict(cb_poly)]
  u3 = [poly_mul(sb_poly, cp3), poly_mul(sb_poly, sp3), dict(cb_poly)]
  return [u1, u2, u3]


def derive_xy3_keo(N_max=4):
  t0 = time.time()
  print(f"=== Starting Fast Exact XY3 SPAR KEO Derivation (N_max={N_max}) ===")
  k_max = N_max + 2
  u_vecs = build_unit_vectors_sparse(k_max + 1)

  # Build B_hat_deg[d][row][col] for d = 0..k_max
  B_hat_deg = [[[{} for _ in range(9)] for _ in range(9)] for _ in range(k_max + 1)]

  def add_by_deg(row, col, poly):
    for k, v in poly.items():
      d = k[3] + k[4]
      if d <= k_max and not v.is_zero():
        B_hat_deg[d][row][col][k] = B_hat_deg[d][row][col].get(k, ZERO_Q) + v

  for i in range(3):
    ux, uy, uz = u_vecs[i]
    # e_0 x u_i = (0, -uz, uy)
    # e_1 x u_i = (uz, 0, -ux)
    # e_2 x u_i = (-uy, ux, 0)
    rot_cols = [
        [{}, poly_scale(uz, -1), uy],
        [uz, {}, poly_scale(ux, -1)],
        [poly_scale(uy, -1), ux, {}],
    ]
    for k2 in range(3):
      row = 3 * i + k2
      for k1 in range(3):
        add_by_deg(row, k1, rot_cols[k1][k2])
      add_by_deg(row, 3 + i, u_vecs[i][k2])
      add_by_deg(row, 6, poly_diff_q(u_vecs[i][k2], 3))  # d/ds4
      add_by_deg(row, 7, poly_diff_q(u_vecs[i][k2], 4))  # d/ds5
      add_by_deg(row, 8, poly_diff_q(u_vecs[i][k2], 5))  # d/dbeta

  # Exact analytical inverse S_hat_0 = (B_hat_deg[0])^-1 at s4=s5=0:
  # Using the orthonormal spherical triad R_sph = diag(R_1, R_2, R_3):
  # Let's check: in the triad basis (r_hat_i, beta_hat_i, phi_hat_i):
  #   col 3+j (r_j)   = delta_{ij} r_hat_i
  #   col 8   (beta)  = beta_hat_i  (for all i=1,2,3)
  #   col 0   (rot_x) = -s_i0 * beta_hat_i - cos(b)*c_i0 * phi_hat_i
  #   col 1   (rot_y) = +c_i0 * beta_hat_i - cos(b)*s_i0 * phi_hat_i
  #   col 2   (rot_z) = +sin(b) * phi_hat_i
  #   col 6   (s_4)   = +sin(b) * U_{i, s4} * phi_hat_i
  #   col 7   (s_5)   = +sin(b) * U_{i, s5} * phi_hat_i
  # Notice that the 3x3 meridional block (cols 0, 1, 8 on beta_hat_i) has rows:
  #   [-s_i0, +c_i0, 1] for i=1,2,3.
  # Since sum_i s_i0 = 0, sum_i c_i0 = 0, sum_i s_i0^2 = 3/2, sum_i c_i0^2 = 3/2, sum_i s_i0*c_i0 = 0,
  # its inverse on beta_hat_i is simply:
  #   row 0 (rot_x): -(2/3) * s_i0 * beta_hat_i^T
  #   row 1 (rot_y): +(2/3) * c_i0 * beta_hat_i^T
  #   row 8 (beta):  +(1/3) * beta_hat_i^T
  # And the 3x3 azimuthal block (cols 2, 6, 7 on phi_hat_i) has:
  #   sin(b) * [1, U_{i, s4}, U_{i, s5}] * phi_hat_i
  # plus the coupling from rot_x, rot_y: -cos(b)*c_i0*rot_x - cos(b)*s_i0*rot_y!
  # Since U_{C3v} = [[1/sqrt(3), 2/sqrt(6), 0], [1/sqrt(3), -1/sqrt(6), 1/sqrt(2)], [1/sqrt(3), -1/sqrt(6), -1/sqrt(2)]]
  # is orthogonal, the inverse on phi_hat_i is Immediate!
  c_i0 = [ONE_Q, Q23(Fraction(-1, 2)), Q23(Fraction(-1, 2))]
  s_i0 = [ZERO_Q, Q23(0, 0, Fraction(1, 2), 0), Q23(0, 0, Fraction(-1, 2), 0)]
  # Dual vectors w_s4 = 3 * u_s4, w_s5 = 3 * u_s5 such that w_a . u_b = delta_{ab}:
  w_s4 = [
      ZERO_Q,
      Q23(0, 0, 0, Fraction(-1, 2)),
      Q23(0, 0, 0, Fraction(1, 2)),
  ]
  w_s5 = [
      Q23(0, 1, 0, 0),
      Q23(0, Fraction(-1, 2), 0, 0),
      Q23(0, Fraction(-1, 2), 0, 0),
  ]

  sb_p = {(0, 0, 0, 0, 0, 1, 0): ONE_Q}
  cb_p = {(0, 0, 0, 0, 0, 0, 1): ONE_Q}
  cscb_p = {(0, 0, 0, 0, 0, -1, 0): ONE_Q}
  cotb_p = {(0, 0, 0, 0, 0, -1, 1): ONE_Q}  # cos(b)/sin(b)

  S_hat_0 = [[{} for _ in range(9)] for _ in range(9)]
  for i in range(3):
    r_hat = [
        poly_scale(sb_p, c_i0[i]),
        poly_scale(sb_p, s_i0[i]),
        dict(cb_p),
    ]
    b_hat = [
        poly_scale(cb_p, c_i0[i]),
        poly_scale(cb_p, s_i0[i]),
        poly_scale(sb_p, -1),
    ]
    p_hat = [
        {(0, 0, 0, 0, 0, 0, 0): -s_i0[i]} if not s_i0[i].is_zero() else {},
        {(0, 0, 0, 0, 0, 0, 0): c_i0[i]},
        {},
    ]
    for k2 in range(3):
      col = 3 * i + k2
      # Rows 0, 1, 8 (rot_x, rot_y, beta): purely from b_hat
      S_hat_0[0][col] = poly_scale(b_hat[k2], c_i0[i] * 0 + s_i0[i] * Fraction(-2, 3))
      S_hat_0[1][col] = poly_scale(b_hat[k2], c_i0[i] * Fraction(2, 3))
      S_hat_0[8][col] = poly_scale(b_hat[k2], Fraction(1, 3))
      # Rows 3, 4, 5 (stretches r1, r2, r3): purely from r_hat
      S_hat_0[3 + i][col] = dict(r_hat[k2])
      # Row 2 (rot_z): (1/(3*sin(b))) * p_hat[k2] (since sum_i c_i0 = sum_i s_i0 = 0, rot_z doesn't couple to rot_x, rot_y!)
      S_hat_0[2][col] = poly_mul(cscb_p, poly_scale(p_hat[k2], Fraction(1, 3)))
      # Rows 6, 7 (s4, s5):
      # sum_j w_s4[j] * (p_hat_j . du_j) = sin(b)*s4 - cos(b)*(sum_j w_s4[j]*s_j0)*rot_y
      # Since sum_j w_s4[j]*s_j0 = -3*sqrt(2)/2:
      #   s4 = csc(b)*sum_j w_s4[j]*p_hat_j - cot(b)*(3*sqrt(2)/2)*rot_y
      # And sum_j w_s5[j]*c_j0 = +3*sqrt(2)/2:
      #   s5 = csc(b)*sum_j w_s5[j]*p_hat_j + cot(b)*(3*sqrt(2)/2)*rot_x
      p_s4 = poly_mul(cscb_p, poly_scale(p_hat[k2], w_s4[i]))
      rot_y_part = poly_mul(cotb_p, poly_scale(S_hat_0[1][col], Q23(0, Fraction(-3, 2), 0, 0)))
      S_hat_0[6][col] = poly_add(p_s4, rot_y_part, ONE_Q)

      p_s5 = poly_mul(cscb_p, poly_scale(p_hat[k2], w_s5[i]))
      rot_x_part = poly_mul(cotb_p, poly_scale(S_hat_0[0][col], Q23(0, Fraction(3, 2), 0, 0)))
      S_hat_0[7][col] = poly_add(p_s5, rot_x_part, ONE_Q)

  # Verify B_hat_deg[0] * S_hat_0 == I_9!
  for r in range(9):
    for c in range(9):
      acc = {}
      for m in range(9):
        acc = poly_add(acc, poly_mul(B_hat_deg[0][r][m], S_hat_0[m][c]))
      expected = {(0, 0, 0, 0, 0, 0, 0): ONE_Q} if r == c else {}
      assert acc == expected, f"S_hat_0 check failed at ({r},{c}): got {acc}"
  print(f"[{time.time() - t0:.3f}s] Exact S_hat_0 constructed and verified (B0 * S0 == I_9)!")

  # Recursive power series inversion: S_hat_deg[d] = - S_hat_0 * sum_{m=1}^d B_hat_deg[m] * S_hat_deg[d-m]
  S_hat_deg = [[[{} for _ in range(9)] for _ in range(9)] for _ in range(k_max + 1)]
  S_hat_deg[0] = S_hat_0
  for d in range(1, k_max + 1):
    accum = [[{} for _ in range(9)] for _ in range(9)]
    for m in range(1, d + 1):
      Bm = B_hat_deg[m]
      Sdm = S_hat_deg[d - m]
      for r in range(9):
        for k in range(9):
          if not Bm[r][k]:
            continue
          for c in range(9):
            if not Sdm[k][c]:
              continue
            accum[r][c] = poly_add(accum[r][c], poly_mul(Bm[r][k], Sdm[k][c], max_deg=d))
    for r in range(9):
      for k in range(9):
        if not S_hat_0[r][k]:
          continue
        for c in range(9):
          if not accum[k][c]:
            continue
          S_hat_deg[d][r][c] = poly_add(
              S_hat_deg[d][r][c], poly_mul(S_hat_0[r][k], accum[k][c], max_deg=d), -ONE_Q
          )
    print(f"[{time.time() - t0:.3f}s] Computed S_hat degree {d}/{k_max}")

  # Combine S_hat_deg into full st[l][n][k2]
  shat_row_for_l = [3, 4, 5, 6, 7, 8, 0, 1, 2]
  ir_polys = [
      {(1, 0, 0, 0, 0, 0, 0): ONE_Q},
      {(0, 1, 0, 0, 0, 0, 0): ONE_Q},
      {(0, 0, 1, 0, 0, 0, 0): ONE_Q},
  ]

  st = [[[{} for _ in range(3)] for _ in range(4)] for _ in range(9)]
  for l in range(9):
    s_row = shat_row_for_l[l]
    max_d_l = (N_max + 2) if l < 6 else (N_max + 1)
    for i in range(3):
      for k2 in range(3):
        col = 3 * i + k2
        combined = {}
        for d in range(max_d_l + 1):
          combined = poly_add(combined, S_hat_deg[d][s_row][col])
        if l not in (0, 1, 2):
          combined = poly_mul(combined, ir_polys[i])
        st[l][i + 1][k2] = combined
    for k2 in range(3):
      s_sum = poly_add(poly_add(st[l][1][k2], st[l][2][k2]), st[l][3][k2])
      st[l][0][k2] = poly_scale(s_sum, -ONE_Q)

  st_N = [[[poly_truncate(st[l][n][k2], N_max) for k2 in range(3)] for n in range(4)] for l in range(9)]
  st_N1 = [[[poly_truncate(st[l][n][k2], N_max + 1) for k2 in range(3)] for n in range(4)] for l in range(9)]
  print(f"[{time.time() - t0:.3f}s] Full st[9, 4, 3] constructed. Computing G_uv and U...")

  # Compute g_vib, g_rot, g_cor
  g_vib = [[[{} for _ in range(2)] for _ in range(6)] for _ in range(6)]
  g_rot = [[[{} for _ in range(2)] for _ in range(3)] for _ in range(3)]
  g_cor = [[[{} for _ in range(2)] for _ in range(3)] for _ in range(6)]

  for l1 in range(9):
    for l2 in range(l1, 9):
      sum_X = {}
      for x1 in range(3):
        sum_X = poly_add(sum_X, poly_mul(st_N[l1][0][x1], st_N[l2][0][x1], max_deg=N_max))
      sum_Y = {}
      for n in (1, 2, 3):
        for x1 in range(3):
          sum_Y = poly_add(sum_Y, poly_mul(st_N[l1][n][x1], st_N[l2][n][x1], max_deg=N_max))

      for idx, val in enumerate((sum_X, sum_Y)):
        if l1 < 6 and l2 < 6:
          g_vib[l1][l2][idx] = val
          g_vib[l2][l1][idx] = val
        elif l1 >= 6 and l2 >= 6:
          g_rot[l1 - 6][l2 - 6][idx] = val
          g_rot[l2 - 6][l1 - 6][idx] = val
        elif l1 < 6 and l2 >= 6:
          g_cor[l1][l2 - 6][idx] = val

  print(f"[{time.time() - t0:.3f}s] G_vib, G_rot, G_cor computed. Computing pseudopotential U...")

  def epsil(a, b, c):
    if (a, b, c) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
      return 1
    if (a, b, c) in ((2, 1, 0), (0, 2, 1), (1, 0, 2)):
      return -1
    return 0

  pseudo = [{}, {}]
  for n in range(4):
    idx = 0 if n == 0 else 1
    # Part 1: u1
    u1_n = {}
    for x0 in range(3):
      curl_x0 = {}
      for q1 in range(3):
        for x1 in range(3):
          e = epsil(x0, q1, x1)
          if e != 0:
            curl_x0 = poly_add(curl_x0, st_N[6 + q1][n][x1], Fraction(e))
      u1_n = poly_add(u1_n, poly_mul(curl_x0, curl_x0, max_deg=N_max), Fraction(1, 8))

    # Part 2: u2
    u2_n = {}
    for x0 in range(3):
      for q1 in range(3):
        for x1 in range(3):
          e = epsil(x0, q1, x1)
          if e != 0:
            for q2 in range(6):
              src_rot = st_N1[6 + q1][n][x0] if q2 in (3, 4) else st_N[6 + q1][n][x0]
              ds_rot = poly_truncate(poly_diff_q(src_rot, q2), N_max)
              u2_n = poly_add(
                  u2_n, poly_mul(ds_rot, st_N[q2][n][x1], max_deg=N_max), Fraction(e, 4)
              )

    # Part 3a (u3) and Part 3b (u4)
    u3_n = {}
    u4_n = {}
    for x0 in range(3):
      div_x0 = {}
      for q2 in range(6):
        src_vib = st[q2][n][x0] if q2 in (3, 4) else st_N1[q2][n][x0]
        ds_q2 = poly_truncate(poly_diff_q(src_vib, q2), N_max + 1)
        ds_q2_N = poly_truncate(ds_q2, N_max)
        div_x0 = poly_add(div_x0, ds_q2_N)
        for q1 in range(6):
          src_ds = ds_q2 if q1 in (3, 4) else ds_q2_N
          dds_q1_q2 = poly_truncate(poly_diff_q(src_ds, q1), N_max)
          u3_n = poly_add(
              u3_n, poly_mul(st_N[q1][n][x0], dds_q1_q2, max_deg=N_max), Fraction(-1, 4)
          )
      u4_n = poly_add(u4_n, poly_mul(div_x0, div_x0, max_deg=N_max), Fraction(-1, 8))

    pseudo[idx] = poly_add(pseudo[idx], poly_add(poly_add(u1_n, u2_n), poly_add(u3_n, u4_n)))

  print(f"[{time.time() - t0:.3f}s] Pseudopotential U computed!")
  return g_vib, g_rot, g_cor, pseudo


BETA_MODE_MAP = {
    (0, 0): 0,   # 1
    (-2, 0): 1,  # csc^2(beta)
    (-1, 1): 2,  # cot(beta)
    (1, 1): 3,   # sin(beta)*cos(beta)
    (2, 0): 4,   # sin^2(beta)
}


def write_combination_chk(filepath: str, N_max: int = 4) -> None:
  """Write combination.chk in TROVE moltype.f90 read_basic_function_constructor format."""
  lines = [
      "Start basic function constructor",
      "Mode 1 2",
      "  1 1  -1 I  1.00000000000000  1",
      "  2 1  -2 I  1.00000000000000  1",
      "Mode 2 2",
      "  1 1  -1 I  1.00000000000000  1",
      "  2 1  -2 I  1.00000000000000  1",
      "Mode 3 2",
      "  1 1  -1 I  1.00000000000000  1",
      "  2 1  -2 I  1.00000000000000  1",
      f"Mode 4 {N_max}",
  ]
  for k in range(1, N_max + 1):
    lines.append(f"  {k} 1   {k} I  1.00000000000000  1")
  lines.append(f"Mode 5 {N_max}")
  for k in range(1, N_max + 1):
    lines.append(f"  {k} 1   {k} I  1.00000000000000  1")
  lines.extend([
      "Mode 6 4",
      "  1 1   2 Csc  1.00000000000000  1",
      "  2 1   1 Cot  1.00000000000000  1",
      "  3 2   1 Sin  1.00000000000000  1   1 Cos  1.00000000000000  1",
      "  4 1   2 Sin  1.00000000000000  1",
      "End basic function constructor",
  ])
  with open(filepath, "w", encoding="utf-8") as f:
    f.write("\n".join(lines) + "\n")


def write_kinetic_chk(
    filepath: str,
    g_vib,
    g_rot,
    g_cor,
    pseudo,
    mass_split: bool = True,
    mX: float = 30.9737620,
    mY: float = 1.00782505,
    thresh: float = 1.0e-6,
) -> int:
  """Write kinetic.chk for TROVE (either mass-independent i0..i6 or mass-baked i1..i6)."""
  lines = []
  n_written = 0

  def emit_poly_pair(k1: int, k2: int, poly_pair) -> None:
    nonlocal n_written
    icount = 0
    if mass_split:
      for idx in (0, 1):
        i0 = idx + 1  # i0=1 -> 1/mX (atom 1), i0=2 -> 1/mY (atom 2)
        for key in sorted(poly_pair[idx].keys()):
          val = poly_pair[idx][key].to_float()
          if abs(val) <= 1e-15:
            continue
          p1, p2, p3, p4, p5, m_sin, p_cb = key
          i6 = BETA_MODE_MAP[(m_sin, p_cb)]
          icount += 1
          n_written += 1
          lines.append(
              f"{k1:5d}{k2:5d}{icount:8d}  0 {val:23.14e}"
              f" {i0:3d}{p1:3d}{p2:3d}{p3:3d}{p4:3d}{p5:3d}{i6:3d}"
          )
    else:
      combined = {}
      for idx, inv_m in ((0, 1.0 / mX), (1, 1.0 / mY)):
        for key, q_val in poly_pair[idx].items():
          combined[key] = combined.get(key, 0.0) + inv_m * q_val.to_float()
      for key in sorted(combined.keys()):
        val = combined[key]
        if abs(val) <= 1e-15:
          continue
        p1, p2, p3, p4, p5, m_sin, p_cb = key
        i6 = BETA_MODE_MAP[(m_sin, p_cb)]
        icount += 1
        n_written += 1
        lines.append(
            f"{k1:5d}{k2:5d}{icount:8d}  0 {val:23.14e}"
            f" {p1:3d}{p2:3d}{p3:3d}{p4:3d}{p5:3d}{i6:3d}"
        )

  sentinel = (
      "   987654321  0    0 0          0.00000000   0  0  0  0  0  0  0"
      if mass_split
      else "   987654321  0    0 0          0.00000000   0  0  0  0  0  0"
  )

  lines.append("  1  1  1   <-  g_vib   ")
  for i in range(6):
    for j in range(6):
      emit_poly_pair(i + 1, j + 1, g_vib[i][j])
  lines.append(sentinel)

  lines.append("  1  1  1   <-  g_rot   ")
  for i in range(3):
    for j in range(3):
      emit_poly_pair(i + 1, j + 1, g_rot[i][j])
  lines.append(sentinel)

  lines.append("  1  1  1   <-  g_cor   ")
  for i in range(6):
    for j in range(3):
      emit_poly_pair(i + 1, j + 1, g_cor[i][j])
  lines.append(sentinel)

  lines.append("  1  1  1   <-  pseudo   ")
  emit_poly_pair(1, 1, pseudo)
  lines.append(sentinel)

  lines.append(f" {thresh:14.8e} <- sparse threshold used")
  lines.append("End of kinetic")

  with open(filepath, "w", encoding="utf-8") as f:
    f.write("\n".join(lines) + "\n")
  return n_written


if __name__ == "__main__":
  import os
  out_dir = os.path.dirname(os.path.abspath(__file__))
  g_vib, g_rot, g_cor, pseudo = derive_xy3_keo(N_max=4)
  # Collect all unique (m_sin, p_cb) pairs across the entire KEO
  beta_modes = set()
  total_terms = 0
  for mat in (g_vib, g_rot, g_cor):
    for row in mat:
      for cell in row:
        for p in cell:
          total_terms += len(p)
          for k in p:
            beta_modes.add((k[5], k[6]))
  for p in pseudo:
    total_terms += len(p)
    for k in p:
      beta_modes.add((k[5], k[6]))
  print(f"Total non-zero SPAR terms at N_max=4: {total_terms}")
  print(f"Unique beta basis functions (m_sin, p_cb): {sorted(beta_modes)}")

  comb_path = os.path.join(out_dir, "combination.chk")
  write_combination_chk(comb_path, N_max=4)
  kin_mass_path = os.path.join(out_dir, "kinetic.chk")
  n_mass = write_kinetic_chk(kin_mass_path, g_vib, g_rot, g_cor, pseudo, mass_split=True)
  kin_nomass_path = os.path.join(out_dir, "kinetic_nomass.chk")
  n_nomass = write_kinetic_chk(kin_nomass_path, g_vib, g_rot, g_cor, pseudo, mass_split=False)
  print(f"Wrote {comb_path}")
  print(f"Wrote {kin_mass_path} ({n_mass} mass-split SPAR terms)")
  print(f"Wrote {kin_nomass_path} ({n_nomass} PH3 mass-combined SPAR terms)")

  # Numerical verification against direct finite-difference 12x12 Cartesian t-s matrix inversion
  import numpy as np

  def eval_spar_poly(poly_pair, r1_v, r2_v, r3_v, s4_v, s5_v, beta_v, inv_mX, inv_mY):
    sb_v = math.sin(beta_v)
    cb_v = math.cos(beta_v)
    total = 0.0
    for idx, inv_m in enumerate((inv_mX, inv_mY)):
      p_val = 0.0
      for (p1, p2, p3, p4, p5, m_sin, p_cb), coeff in poly_pair[idx].items():
        term = (
            coeff.to_float()
            * (r1_v ** (-p1))
            * (r2_v ** (-p2))
            * (r3_v ** (-p3))
            * (s4_v**p4)
            * (s5_v**p5)
            * (sb_v**m_sin)
            * (cb_v**p_cb)
        )
        p_val += term
      total += inv_m * p_val
    return total

  def cart_coords(q, mX, mY):
    r1_v, r2_v, r3_v, s4_v, s5_v, beta_v = q
    d1 = math.sqrt(2.0) * s5_v / 3.0
    d2 = -s4_v / math.sqrt(6.0) - s5_v / (3.0 * math.sqrt(2.0))
    d3 = +s4_v / math.sqrt(6.0) - s5_v / (3.0 * math.sqrt(2.0))
    phi = [d1, 2.0 * math.pi / 3.0 + d2, 4.0 * math.pi / 3.0 + d3]
    r_arr = [r1_v, r2_v, r3_v]
    sb_v, cb_v = math.sin(beta_v), math.cos(beta_v)
    pos = np.zeros((4, 3))
    for i in range(3):
      pos[i + 1, 0] = r_arr[i] * sb_v * math.cos(phi[i])
      pos[i + 1, 1] = r_arr[i] * sb_v * math.sin(phi[i])
      pos[i + 1, 2] = r_arr[i] * cb_v
    masses = np.array([mX, mY, mY, mY])
    cm = np.sum(masses[:, None] * pos, axis=0) / np.sum(masses)
    return pos - cm

  def s_matrix_num(q, mX, mY):
    # Exact analytical t_mat at q (no finite-difference noise in t_mat!)
    r1_v, r2_v, r3_v, s4_v, s5_v, beta_v = q
    d1 = math.sqrt(2.0) * s5_v / 3.0
    d2 = -s4_v / math.sqrt(6.0) - s5_v / (3.0 * math.sqrt(2.0))
    d3 = +s4_v / math.sqrt(6.0) - s5_v / (3.0 * math.sqrt(2.0))
    phi = [d1, 2.0 * math.pi / 3.0 + d2, 4.0 * math.pi / 3.0 + d3]
    dphi_ds4 = [0.0, -1.0 / math.sqrt(6.0), 1.0 / math.sqrt(6.0)]
    dphi_ds5 = [math.sqrt(2.0) / 3.0, -1.0 / (3.0 * math.sqrt(2.0)), -1.0 / (3.0 * math.sqrt(2.0))]
    r_arr = [r1_v, r2_v, r3_v]
    sb_v, cb_v = math.sin(beta_v), math.cos(beta_v)

    pos = np.zeros((4, 3))
    u_v = np.zeros((3, 3))
    du_db = np.zeros((3, 3))
    du_dphi = np.zeros((3, 3))
    for i in range(3):
      cp, sp_i = math.cos(phi[i]), math.sin(phi[i])
      u_v[i] = [sb_v * cp, sb_v * sp_i, cb_v]
      du_db[i] = [cb_v * cp, cb_v * sp_i, -sb_v]
      du_dphi[i] = [-sb_v * sp_i, sb_v * cp, 0.0]
      pos[i + 1] = r_arr[i] * u_v[i]

    masses = np.array([mX, mY, mY, mY])
    M_tot = np.sum(masses)
    cm = np.sum(masses[:, None] * pos, axis=0) / M_tot
    pos_cm = pos - cm

    t_mat = np.zeros((12, 12))
    # cols 0..2: stretches r1, r2, r3
    for i in range(3):
      dp = np.zeros((4, 3))
      dp[i + 1] = u_v[i]
      dp -= (mY / M_tot) * u_v[i]
      t_mat[:, i] = dp.reshape(12)
    # col 3: s4
    dp_s4 = np.zeros((4, 3))
    for i in range(3):
      dp_s4[i + 1] = r_arr[i] * dphi_ds4[i] * du_dphi[i]
    dp_s4 -= np.sum(masses[:, None] * dp_s4, axis=0) / M_tot
    t_mat[:, 3] = dp_s4.reshape(12)
    # col 4: s5
    dp_s5 = np.zeros((4, 3))
    for i in range(3):
      dp_s5[i + 1] = r_arr[i] * dphi_ds5[i] * du_dphi[i]
    dp_s5 -= np.sum(masses[:, None] * dp_s5, axis=0) / M_tot
    t_mat[:, 4] = dp_s5.reshape(12)
    # col 5: beta
    dp_b = np.zeros((4, 3))
    for i in range(3):
      dp_b[i + 1] = r_arr[i] * du_db[i]
    dp_b -= np.sum(masses[:, None] * dp_b, axis=0) / M_tot
    t_mat[:, 5] = dp_b.reshape(12)
    # cols 6..8: rotational e_alpha x pos_cm
    e_vecs = [np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, 1.0])]
    for a in range(3):
      t_mat[:, 6 + a] = np.cross(e_vecs[a], pos_cm).reshape(12)
    # cols 9..11: translational e_alpha
    for a in range(3):
      t_mat[:, 9 + a] = np.tile(e_vecs[a], (4, 1)).reshape(12)
    s_full = np.linalg.inv(t_mat)
    return s_full[:9, :].reshape(9, 4, 3)

  mX, mY = 30.973762, 1.00782505
  inv_m = np.array([1.0 / mX, 1.0 / mY, 1.0 / mY, 1.0 / mY])

  for label, q_test in [
      ("s4=s5=0, asymmetric r1,r2,r3", [1.38, 1.42, 1.47, 0.0, 0.0, 0.98]),
      ("s4=0.02, s5=-0.015, asymmetric r", [1.38, 1.42, 1.47, 0.02, -0.015, 0.98]),
  ]:
    s_num = s_matrix_num(q_test, mX, mY)
    g_num = np.einsum("n,lnk,mnk->lm", inv_m, s_num, s_num)
    max_err_gvib = 0.0
    for i in range(6):
      for j in range(6):
        v_spar = eval_spar_poly(g_vib[i][j], *q_test, 1.0 / mX, 1.0 / mY)
        max_err_gvib = max(max_err_gvib, abs(v_spar - g_num[i, j]))
    max_err_grot = 0.0
    for i in range(3):
      for j in range(3):
        v_spar = eval_spar_poly(g_rot[i][j], *q_test, 1.0 / mX, 1.0 / mY)
        max_err_grot = max(max_err_grot, abs(v_spar - g_num[6 + i, 6 + j]))
    max_err_gcor = 0.0
    for i in range(6):
      for j in range(3):
        v_spar = eval_spar_poly(g_cor[i][j], *q_test, 1.0 / mX, 1.0 / mY)
        max_err_gcor = max(max_err_gcor, abs(v_spar - g_num[i, 6 + j]))

    eps = 1e-5
    ds_num = np.zeros((6, 9, 4, 3))
    dds_num = np.zeros((6, 6, 9, 4, 3))
    for k1 in range(6):
      qp = list(q_test)
      qm = list(q_test)
      qp[k1] += eps
      qm[k1] -= eps
      sp_k1 = s_matrix_num(qp, mX, mY)
      sm_k1 = s_matrix_num(qm, mX, mY)
      ds_num[k1] = (sp_k1 - sm_k1) / (2.0 * eps)
      dds_num[k1, k1] = (sp_k1 - 2.0 * s_num + sm_k1) / (eps**2)
      for k2 in range(k1 + 1, 6):
        qpp, qpm, qmp, qmm = list(q_test), list(q_test), list(q_test), list(q_test)
        qpp[k1] += eps; qpp[k2] += eps
        qpm[k1] += eps; qpm[k2] -= eps
        qmp[k1] -= eps; qmp[k2] += eps
        qmm[k1] -= eps; qmm[k2] -= eps
        val_dd = (
            s_matrix_num(qpp, mX, mY)
            - s_matrix_num(qpm, mX, mY)
            - s_matrix_num(qmp, mX, mY)
            + s_matrix_num(qmm, mX, mY)
        ) / (4.0 * eps * eps)
        dds_num[k1, k2] = val_dd
        dds_num[k2, k1] = val_dd

    lev = np.zeros((3, 3, 3))
    lev[0, 1, 2] = lev[1, 2, 0] = lev[2, 0, 1] = 1.0
    lev[2, 1, 0] = lev[0, 2, 1] = lev[1, 0, 2] = -1.0
    s_vib_n = s_num[:6]
    s_rot_n = s_num[6:]
    u1_num = 0.125 * np.einsum("n,xqa,xrb,qna,rnb->", inv_m, lev, lev, s_rot_n, s_rot_n)
    u2_num = 0.25 * np.einsum("n,xqa,kqnx,kna->", inv_m, lev, ds_num[:, 6:, :, :], s_vib_n)
    u3_num = -0.25 * np.einsum("n,knx,kllnx->", inv_m, s_vib_n, dds_num[:, :, :6, :, :])
    u4_num = -0.125 * np.einsum("n,kkna,llna->", inv_m, ds_num[:, :6, :, :], ds_num[:, :6, :, :])
    u_tot_num = u1_num + u2_num + u3_num + u4_num
    u_spar = eval_spar_poly(pseudo, *q_test, 1.0 / mX, 1.0 / mY)
    print(f"Verification [{label}]:")
    print(f"  max |G_vib - G_vib_num| = {max_err_gvib:.3e}")
    print(f"  max |G_rot - G_rot_num| = {max_err_grot:.3e}")
    print(f"  max |G_cor - G_cor_num| = {max_err_gcor:.3e}")
    print(f"  |U_spar - U_num|        = {abs(u_spar - u_tot_num):.3e} (U_spar={u_spar:.8f}, U_num={u_tot_num:.8f})")
