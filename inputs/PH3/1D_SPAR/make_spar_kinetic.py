#!/usr/bin/env python3
"""Writes the 1D SPAR (Taylor in s4,s5 to total degree N_max) kinetic.chk for PH3 (R-SYMPHI-BETA).

Usage: python3 make_spar_kinetic.py N_max [output_file]
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from generate_xy3_spar_keo import derive_xy3_keo, write_kinetic_chk

n_max = int(sys.argv[1])
out = sys.argv[2] if len(sys.argv) > 2 else "kinetic.chk"
g_vib, g_rot, g_cor, u = derive_xy3_keo(N_max=n_max)
n_terms = write_kinetic_chk(out, g_vib, g_rot, g_cor, u, mass_split=True)
print(f"N_max={n_max}: {n_terms} terms written to {out}")
