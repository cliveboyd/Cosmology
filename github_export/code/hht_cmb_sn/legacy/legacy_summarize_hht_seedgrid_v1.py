#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 11 22:41:22 2026

@author: boyde



"""

#!/usr/bin/env python3
# summarize_hht_seedgrid_v1.py

import pandas  as     pd
from   pathlib import Path

CSV = Path("/Users/boyde/.spyder-py3/plamb_runs/hht_bridge/FR_log1pz_hht_export_v4/hht_seasonal_robustness_summary_v4.csv")

df  = pd.read_csv(CSV)

# If your CSV contains multiple kinds of rows, filter by what run_hht_seedgrid writes.
# Common patterns: a 'seed' column, or a 'tag', etc. Adjust as needed:
# df = df[df["seed"].notna()].copy()

print("\n[loaded]", CSV)
print("[rows]", len(df), " [cols]", len(df.columns))

# Choose the columns you care about (adjust names to match your schema-4 header)
cols = [
        "PLV_p_raw",
        "PLV_p_det",
        "Rayleigh_p_raw",
        "Rayleigh_p_det",
        "wald_p_det",
        "runL_best_p",
        "runningL_spin_global_p_det",
       ]

present = [c for c in cols if c in df.columns]
missing = [c for c in cols if c not in df.columns]

print("\n[present cols]", present)

if missing:
    print("[missing cols]", missing)

# Convert to numeric safely
for c in present:
    df[c] = pd.to_numeric(df[c], errors="coerce")

def stability(col):
    s = df[col].dropna()
    if len(s) == 0:
        return {"n": 0}
    return {
            "n"              : int(s.shape[0]),
            "min"            : float(s.min()),
            "q25"            : float(s.quantile(0.25)),
            "median"         : float(s.median()),
            "q75"            : float(s.quantile(0.75)),
            "max"            : float(s.max()),
            "frac_p_lt_0p05" : float((s < 0.05).mean()),
            "frac_p_lt_0p01" : float((s < 0.01).mean()),
        }

print("\n[stability]")
for c in present:
    print("\n", c)
    print(stability(c))

# Quick NA check
if "runningL_spin_global_p_det" in df.columns:
    print("\n[runningL_spin_global_p_det NA fraction]",
          df["runningL_spin_global_p_det"].isna().mean())
