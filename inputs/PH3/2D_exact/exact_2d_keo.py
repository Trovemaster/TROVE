#!/usr/bin/env python3
"""Generates exact algebraic 2D-SPAR KEO (kinetic.chk and BASIC-FUNCTION mass block) for XY3 (R-SYMPHI-BETA)."""

import argparse
from collections import defaultdict
import math
import os
import time

BETA_MODE_MAP = {
    (0, 0): 0,   # 1
    (-2, 0): 1,  # csc^2(beta)
    (-1, 1): 2,  # cot(beta)
    (1, 1): 3,   # sin(beta)*cos(beta)
    (2, 0): 4,   # sin^2(beta)
}

a4 = 1.0 / math.sqrt(6.0)
a5 = 1.0 / (3.0 * math.sqrt(2.0))
s3_2 = math.sqrt(3.0) / 2.0


def trig_product_1d_flt(tA, kA, tB, kB):
  if kA == 0:
    return [] if tA == "S" else [(tB, kB, 1.0)]
  if kB == 0:
    return [] if tB == "S" else [(tA, kA, 1.0)]
  ksum = kA + kB
  kdiff = abs(kA - kB)
  sgn_diff = 1.0 if kA >= kB else -1.0
  out = []
  if tA == "C" and tB == "C":
    out.append(("C", kdiff, 0.5))
    out.append(("C", ksum, 0.5))
  elif tA == "S" and tB == "S":
    out.append(("C", kdiff, 0.5))
    out.append(("C", ksum, -0.5))
  elif tA == "S" and tB == "C":
    out.append(("S", ksum, 0.5))
    if kdiff > 0:
      out.append(("S", kdiff, 0.5 * sgn_diff))
  elif tA == "C" and tB == "S":
    out.append(("S", ksum, 0.5))
    if kdiff > 0:
      out.append(("S", kdiff, -0.5 * sgn_diff))
  return out


def fpoly_clean(p, thresh=1e-13):
  return {k: v for k, v in p.items() if abs(v) > thresh}


def fpoly_add(pA, pB, scaleB=1.0):
  res = dict(pA)
  for k, v in pB.items():
    nv = res.get(k, 0.0) + v * scaleB
    if abs(nv) <= 1e-13:
      res.pop(k, None)
    else:
      res[k] = nv
  return res


def fpoly_scale(pA, c):
  if abs(c) <= 1e-13:
    return {}
  return {k: v * c for k, v in pA.items()}


def fpoly_mul(pA, pB):
  if not pA or not pB:
    return {}
  res = defaultdict(float)
  for (p1a, p2a, p3a, ma, ca, da, t4a, k4a, t5a, k5a), va in pA.items():
    for (p1b, p2b, p3b, mb, cb, db, t4b, k4b, t5b, k5b), vb in pB.items():
      prod0 = va * vb
      p1 = p1a + p1b
      p2 = p2a + p2b
      p3 = p3a + p3b
      m = ma + mb
      d = da + db
      t4_list = trig_product_1d_flt(t4a, k4a, t4b, k4b)
      t5_list = trig_product_1d_flt(t5a, k5a, t5b, k5b)
      for t4, k4, c4 in t4_list:
        for t5, k5, c5 in t5_list:
          prod = prod0 * c4 * c5
          if ca + cb < 2:
            res[(p1, p2, p3, m, ca + cb, d, t4, k4, t5, k5)] += prod
          else:
            res[(p1, p2, p3, m, 0, d, t4, k4, t5, k5)] += prod
            res[(p1, p2, p3, m + 2, 0, d, t4, k4, t5, k5)] -= prod
  return fpoly_clean(res)


def derive_exact_2d_spar_keo():
  c_phi1 = {(0, 0, 0, 0, 0, 0, "C", 0, "C", 2): 1.0}
  s_phi1 = {(0, 0, 0, 0, 0, 0, "C", 0, "S", 2): 1.0}
  c_d2 = {
      (0, 0, 0, 0, 0, 0, "C", 1, "C", 1): 1.0,
      (0, 0, 0, 0, 0, 0, "S", 1, "S", 1): -1.0,
  }
  s_d2 = {
      (0, 0, 0, 0, 0, 0, "S", 1, "C", 1): -1.0,
      (0, 0, 0, 0, 0, 0, "C", 1, "S", 1): -1.0,
  }
  c_d3 = {
      (0, 0, 0, 0, 0, 0, "C", 1, "C", 1): 1.0,
      (0, 0, 0, 0, 0, 0, "S", 1, "S", 1): 1.0,
  }
  s_d3 = {
      (0, 0, 0, 0, 0, 0, "S", 1, "C", 1): 1.0,
      (0, 0, 0, 0, 0, 0, "C", 1, "S", 1): -1.0,
  }

  c_phi2 = fpoly_add(fpoly_scale(c_d2, -0.5), fpoly_scale(s_d2, -s3_2))
  s_phi2 = fpoly_add(fpoly_scale(c_d2, s3_2), fpoly_scale(s_d2, -0.5))
  c_phi3 = fpoly_add(fpoly_scale(c_d3, -0.5), fpoly_scale(s_d3, s3_2))
  s_phi3 = fpoly_add(fpoly_scale(c_d3, -s3_2), fpoly_scale(s_d3, -0.5))

  c_phi = [c_phi1, c_phi2, c_phi3]
  s_phi = [s_phi1, s_phi2, s_phi3]

  Delta_poly = {}
  for i in range(3):
    j, k = (i + 1) % 3, (i + 2) % 3
    term = fpoly_add(
        fpoly_mul(s_phi[k], c_phi[j]), fpoly_mul(c_phi[k], s_phi[j]), -1.0
    )
    Delta_poly = fpoly_add(Delta_poly, term)

  def fpoly_diff_q_no_delta(pA, q_idx):
    res = defaultdict(float)
    for (p1, p2, p3, m, c, d, t4, k4, t5, k5), v in pA.items():
      if q_idx == 0 and p1 > 0:
        res[(p1 + 1, p2, p3, m, c, d, t4, k4, t5, k5)] += -v * p1
      elif q_idx == 1 and p2 > 0:
        res[(p1, p2 + 1, p3, m, c, d, t4, k4, t5, k5)] += -v * p2
      elif q_idx == 2 and p3 > 0:
        res[(p1, p2, p3 + 1, m, c, d, t4, k4, t5, k5)] += -v * p3
      elif q_idx == 3 and k4 > 0:
        if t4 == "C":
          res[(p1, p2, p3, m, c, d, "S", k4, t5, k5)] += -v * (k4 * a4)
        else:
          res[(p1, p2, p3, m, c, d, "C", k4, t5, k5)] += v * (k4 * a4)
      elif q_idx == 4 and k5 > 0:
        if t5 == "C":
          res[(p1, p2, p3, m, c, d, t4, k4, "S", k5)] += -v * (k5 * a5)
        else:
          res[(p1, p2, p3, m, c, d, t4, k4, "C", k5)] += v * (k5 * a5)
      elif q_idx == 5:
        if m != 0:
          if c == 0:
            res[(p1, p2, p3, m - 1, 1, d, t4, k4, t5, k5)] += v * m
          else:
            res[(p1, p2, p3, m - 1, 0, d, t4, k4, t5, k5)] += v * m
            res[(p1, p2, p3, m + 1, 0, d, t4, k4, t5, k5)] += -v * m
        if c == 1:
          res[(p1, p2, p3, m + 1, 0, d, t4, k4, t5, k5)] += -v
    return fpoly_clean(res)

  dDelta_ds4 = fpoly_diff_q_no_delta(Delta_poly, 3)
  dDelta_ds5 = fpoly_diff_q_no_delta(Delta_poly, 4)

  def fpoly_diff_q(pA, q_idx):
    res = fpoly_diff_q_no_delta(pA, q_idx)
    if q_idx in (3, 4):
      dDelta = dDelta_ds4 if q_idx == 3 else dDelta_ds5
      chain_part = defaultdict(float)
      for (p1, p2, p3, m, c, d, t4, k4, t5, k5), v in pA.items():
        if d > 0:
          chain_part[(p1, p2, p3, m, c, d + 1, t4, k4, t5, k5)] += -v * d
      if chain_part:
        res = fpoly_add(res, fpoly_mul(fpoly_clean(chain_part), dDelta))
    return res

  inv_Delta_p = {(0, 0, 0, 0, 0, 1, "C", 0, "C", 0): 1.0}
  M_inv = [[{} for _ in range(3)] for _ in range(3)]
  for i in range(3):
    j, k = (i + 1) % 3, (i + 2) % 3
    M_inv[0][i] = fpoly_mul(fpoly_add(c_phi[j], c_phi[k], -1.0), inv_Delta_p)
    M_inv[1][i] = fpoly_mul(fpoly_add(s_phi[j], s_phi[k], -1.0), inv_Delta_p)
    M_inv[2][i] = fpoly_mul(
        fpoly_add(
            fpoly_mul(s_phi[k], c_phi[j]), fpoly_mul(c_phi[k], s_phi[j]), -1.0
        ),
        inv_Delta_p,
    )

  w_rot_z = [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]
  w_s4 = [0.0, -math.sqrt(6.0) / 2.0, math.sqrt(6.0) / 2.0]
  w_s5 = [math.sqrt(2.0), -math.sqrt(2.0) / 2.0, -math.sqrt(2.0) / 2.0]

  sb_p = {(0, 0, 0, 1, 0, 0, "C", 0, "C", 0): 1.0}
  cb_p = {(0, 0, 0, 0, 1, 0, "C", 0, "C", 0): 1.0}
  cscb_p = {(0, 0, 0, -1, 0, 0, "C", 0, "C", 0): 1.0}
  cotb_p = {(0, 0, 0, -1, 1, 0, "C", 0, "C", 0): 1.0}

  r_hat, b_hat, p_hat = [], [], []
  for i in range(3):
    r_hat.append(
        [fpoly_mul(sb_p, c_phi[i]), fpoly_mul(sb_p, s_phi[i]), dict(cb_p)]
    )
    b_hat.append([
        fpoly_mul(cb_p, c_phi[i]),
        fpoly_mul(cb_p, s_phi[i]),
        fpoly_scale(sb_p, -1.0),
    ])
    p_hat.append([fpoly_scale(s_phi[i], -1.0), dict(c_phi[i]), {}])

  S_hat = [[[{} for _ in range(3)] for _ in range(3)] for _ in range(9)]
  for i in range(3):
    for x in range(3):
      S_hat[i][i][x] = r_hat[i][x]
      S_hat[6][i][x] = fpoly_mul(b_hat[i][x], M_inv[0][i])
      S_hat[7][i][x] = fpoly_mul(b_hat[i][x], M_inv[1][i])
      S_hat[5][i][x] = fpoly_mul(b_hat[i][x], M_inv[2][i])

  for l_row, w_vec in ((3, w_s4), (4, w_s5), (8, w_rot_z)):
    wc, ws = {}, {}
    for j in range(3):
      wc = fpoly_add(wc, c_phi[j], w_vec[j])
      ws = fpoly_add(ws, s_phi[j], w_vec[j])
    for i in range(3):
      for x in range(3):
        term1 = fpoly_mul(fpoly_scale(cscb_p, w_vec[i]), p_hat[i][x])
        comb_rot = fpoly_add(
            fpoly_mul(S_hat[6][i][x], wc), fpoly_mul(S_hat[7][i][x], ws)
        )
        S_hat[l_row][i][x] = fpoly_add(term1, fpoly_mul(cotb_p, comb_rot))

  ir_polys = [
      {(1, 0, 0, 0, 0, 0, "C", 0, "C", 0): 1.0},
      {(0, 1, 0, 0, 0, 0, "C", 0, "C", 0): 1.0},
      {(0, 0, 1, 0, 0, 0, "C", 0, "C", 0): 1.0},
  ]
  st = [[[{} for _ in range(3)] for _ in range(4)] for _ in range(9)]
  for l in range(9):
    for i in range(3):
      for x in range(3):
        st[l][i + 1][x] = (
            fpoly_mul(S_hat[l][i][x], ir_polys[i]) if l >= 3 else S_hat[l][i][x]
        )
    for x in range(3):
      st[l][0][x] = fpoly_scale(
          fpoly_add(fpoly_add(st[l][1][x], st[l][2][x]), st[l][3][x]), -1.0
      )

  g_vib = [[[{} for _ in range(2)] for _ in range(6)] for _ in range(6)]
  g_rot = [[[{} for _ in range(2)] for _ in range(3)] for _ in range(3)]
  g_cor = [[[{} for _ in range(2)] for _ in range(3)] for _ in range(6)]
  for l1 in range(9):
    for l2 in range(l1, 9):
      sum_X, sum_Y = {}, {}
      for x in range(3):
        sum_X = fpoly_add(sum_X, fpoly_mul(st[l1][0][x], st[l2][0][x]))
      for n in (1, 2, 3):
        for x in range(3):
          sum_Y = fpoly_add(sum_Y, fpoly_mul(st[l1][n][x], st[l2][n][x]))
      for idx, val in enumerate((sum_X, sum_Y)):
        if l1 < 6 and l2 < 6:
          g_vib[l1][l2][idx] = val
          g_vib[l2][l1][idx] = val
        elif l1 >= 6 and l2 >= 6:
          g_rot[l1 - 6][l2 - 6][idx] = val
          g_rot[l2 - 6][l1 - 6][idx] = val
        elif l1 < 6 and l2 >= 6:
          g_cor[l1][l2 - 6][idx] = val

  def epsil(a, b, c):
    if (a, b, c) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
      return 1
    if (a, b, c) in ((2, 1, 0), (0, 2, 1), (1, 0, 2)):
      return -1
    return 0

  ds = [
      [
          [[fpoly_diff_q(st[l][n][x], q) for q in range(6)] for x in range(3)]
          for n in range(4)
      ]
      for l in range(9)
  ]
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
            curl_x0 = fpoly_add(curl_x0, st[6 + q1][n][x1], float(e))
            for q2 in range(6):
              u2_n = fpoly_add(
                  u2_n,
                  fpoly_mul(ds[6 + q1][n][x0][q2], st[q2][n][x1]),
                  0.25 * e,
              )
      u1_n = fpoly_add(u1_n, fpoly_mul(curl_x0, curl_x0), 0.125)
      div_x0 = {}
      for q2 in range(6):
        div_x0 = fpoly_add(div_x0, ds[q2][n][x0][q2])
      u4_n = fpoly_add(u4_n, fpoly_mul(div_x0, div_x0), -0.125)
      for q1 in range(6):
        u3_n = fpoly_add(
            u3_n,
            fpoly_mul(st[q1][n][x0], fpoly_diff_q(div_x0, q1)),
            -0.25,
        )
    tot_u = fpoly_add(fpoly_add(u1_n, u2_n), fpoly_add(u3_n, u4_n))
    pseudo[idx] = fpoly_add(pseudo[idx], tot_u)

  return g_vib, g_rot, g_cor, pseudo


def build_mode45_maps(g_vib, g_rot, g_cor, pseudo):
  f4_set = set()
  f5_set = set()
  all_dicts = []
  for i in range(6):
    for j in range(6):
      all_dicts.extend(g_vib[i][j])
  for i in range(3):
    for j in range(3):
      all_dicts.extend(g_rot[i][j])
  for i in range(6):
    for j in range(3):
      all_dicts.extend(g_cor[i][j])
  all_dicts.extend(pseudo)

  for d_map in all_dicts:
    for (p1, p2, p3, m, c, d, t4, k4, t5, k5), v in d_map.items():
      if abs(v) > 1e-13:
        f4_set.add((t4, k4))
        f5_set.add((d, t5, k5))

  f4_list = sorted(f4_set - {("C", 0)}, key=lambda x: (x[1], x[0]))
  f5_list = sorted(f5_set - {(0, "C", 0)}, key=lambda x: (x[0], x[2], x[1]))

  f4_map = {("C", 0): 0}
  for idx, item in enumerate(f4_list, 1):
    f4_map[item] = idx

  f5_map = {(0, "C", 0): 0}
  for idx, item in enumerate(f5_list, 1):
    f5_map[item] = idx

  return f4_list, f4_map, f5_list, f5_map


def format_basic_function_block(f4_list, f5_list):
  lines = ["BASIC-FUNCTION mass"]
  for mode_r in (1, 2, 3):
    lines.append(f"Mode {mode_r} 2")
    lines.append("  1 1  I -1 1.0000000000000000")
    lines.append("  2 1  I -2 1.0000000000000000")

  lines.append(f"Mode 4 {len(f4_list)}")
  for idx, (t4, k4) in enumerate(f4_list, 1):
    fn = "Cos" if t4 == "C" else "Sin"
    freq = k4 * a4
    lines.append(f"  {idx:2d} 1  {fn:3s}  1 {freq:22.16f}")

  lines.append(f"Mode 5 {len(f5_list)}")
  for idx, (d, t5, k5) in enumerate(f5_list, 1):
    if d == 0:
      fn = "Cos" if t5 == "C" else "Sin"
      freq = k5 * a5
      lines.append(f"  {idx:2d} 1  {fn:3s}  1 {freq:22.16f}")
    elif k5 == 0 and t5 == "C":
      lines.append(f"  {idx:2d} 1  DELTA_XY3 {-d:2d} {1.0:22.16f}")
    else:
      fn = "Cos" if t5 == "C" else "Sin"
      freq = k5 * a5
      lines.append(
          f"  {idx:2d} 2  {fn:3s}  1 {freq:22.16f}   DELTA_XY3 {-d:2d} {1.0:22.16f}"
      )

  lines.append("Mode 6 4")
  lines.append("  1 1  Csc  2 1.0000000000000000")
  lines.append("  2 1  Cot  1 1.0000000000000000")
  lines.append("  3 2  Sin  1 1.0000000000000000   Cos  1 1.0000000000000000")
  lines.append("  4 1  Sin  2 1.0000000000000000")
  lines.append("END")
  return "\n".join(lines)


def write_kinetic_chk(out_dir, g_vib, g_rot, g_cor, pseudo, f4_map, f5_map):
  lines = []
  n_written = 0

  def emit_section(k1, k2, dicts_pair):
    nonlocal n_written
    icount = 0
    for idx in (0, 1):
      i0 = idx + 1
      d_map = dicts_pair[idx]
      # Map to (p1, p2, p3, i4, i5, i6) and sort
      mapped = []
      for (p1, p2, p3, m, c, d, t4, k4, t5, k5), val in d_map.items():
        if abs(val) > 1e-13:
          i4 = f4_map[(t4, k4)]
          i5 = f5_map[(d, t5, k5)]
          i6 = BETA_MODE_MAP[(m, c)]
          mapped.append(((p1, p2, p3, i4, i5, i6), val))
      mapped.sort(key=lambda x: x[0])
      for (p1, p2, p3, i4, i5, i6), val in mapped:
        icount += 1
        n_written += 1
        lines.append(
            f"{k1:5d}{k2:5d}{icount:8d}  0 {val:23.14e}"
            f" {i0:3d}{p1:3d}{p2:3d}{p3:3d}{i4:3d}{i5:3d}{i6:3d}"
        )

  sentinel = "   987654321  0    0 0          0.00000000   0  0  0  0  0  0  0"
  lines.append("  1  1  1   <-  g_vib   ")
  for i in range(6):
    for j in range(6):
      emit_section(i + 1, j + 1, g_vib[i][j])
  lines.append(sentinel)

  lines.append("  1  1  1   <-  g_rot   ")
  for i in range(3):
    for j in range(3):
      emit_section(i + 1, j + 1, g_rot[i][j])
  lines.append(sentinel)

  lines.append("  1  1  1   <-  g_cor   ")
  for i in range(6):
    for j in range(3):
      emit_section(i + 1, j + 1, g_cor[i][j])
  lines.append(sentinel)

  lines.append("  1  1  1   <-  pseudo   ")
  emit_section(1, 1, pseudo)
  lines.append(sentinel)

  lines.append(f" {1.0e-6:14.8e} <- sparse threshold used")
  lines.append("End of kinetic")

  kin_path = os.path.join(out_dir, "kinetic.chk")
  with open(kin_path, "w", encoding="utf-8") as f:
    f.write("\n".join(lines) + "\n")
  return kin_path, n_written


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("--out-dir", type=str, default=".")
  args = parser.parse_args()

  t0 = time.time()
  os.makedirs(args.out_dir, exist_ok=True)
  g_vib, g_rot, g_cor, pseudo = derive_exact_2d_spar_keo()
  f4_list, f4_map, f5_list, f5_map = build_mode45_maps(
      g_vib, g_rot, g_cor, pseudo
  )
  bf_block = format_basic_function_block(f4_list, f5_list)
  bf_path = os.path.join(args.out_dir, "basic_function_mass.inp")
  with open(bf_path, "w", encoding="utf-8") as f:
    f.write(bf_block + "\n")

  kin_path, n_written = write_kinetic_chk(
      args.out_dir, g_vib, g_rot, g_cor, pseudo, f4_map, f5_map
  )
  print(
      f"[{time.time() - t0:.2f}s] Wrote {kin_path} ({n_written} exact 2D-SPAR terms, "
      f"Mode 4: {len(f4_list)} funcs, Mode 5: {len(f5_list)} funcs) and {bf_path}"
  )


if __name__ == "__main__":
  main()
