#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 22 21:01:55 2025

@author: boyde
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
quaia_inspect_dipoles.py  (0V1)

Read quaia_dipole_summary.txt and add separation from CMB dipole.
"""

import numpy        as     np
from   pathlib      import Path
from   quaia_config import OUT_DIR

# CMB dipole in Galactic coords
L_CMB_DEG = 264.0
B_CMB_DEG = 48.0

def ang_sep(l1_deg, b1_deg, l2_deg, b2_deg):
    l1 = np.deg2rad(l1_deg)
    b1 = np.deg2rad(b1_deg)
    l2 = np.deg2rad(l2_deg)
    b2 = np.deg2rad(b2_deg)

    x1 = np.cos(b1) * np.cos(l1)
    y1 = np.cos(b1) * np.sin(l1)
    z1 = np.sin(b1)

    x2 = np.cos(b2) * np.cos(l2)
    y2 = np.cos(b2) * np.sin(l2)
    z2 = np.sin(b2)

    dot = x1 * x2 + y1 * y2 + z1 * z2
    dot = np.clip(dot, -1.0, 1.0)
    return np.rad2deg(np.arccos(dot))

def main():
    txt_path = OUT_DIR / "quaia_dipole_summary.txt"
    if not txt_path.exists():
        raise SystemExit(f"Missing {txt_path}")

    lines = txt_path.read_text().strip().splitlines()
    print(lines[0])  # header

    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.split()
        tag   = parts[0]
        # last 4 columns are: N_good, amp, l_deg, b_deg, a0
        N_good = int(parts[-5])
        amp    = float(parts[-4])
        l_deg  = float(parts[-3])
        b_deg  = float(parts[-2])
        a0     = float(parts[-1])

        if np.isfinite(l_deg) and np.isfinite(b_deg):
            sep = ang_sep(l_deg, b_deg, L_CMB_DEG, B_CMB_DEG)
        else:
            sep = np.nan

        print(f"{tag:12s} N={N_good:6d}  amp={amp:8.2e}  "
              f"(l,b)=({l_deg:7.2f},{b_deg:7.2f})  "
              f"a0={a0:9.2e}  sep_CMB={sep:6.1f}°")

if __name__ == "__main__":
    main()
