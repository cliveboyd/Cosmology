#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  1 08:30:35 2026

@author: boyde
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

BASE = Path("/Users/boyde/.spyder-py3/quaia_imf_grid")
METR = BASE / "grid_metrics_with_shiftnull.csv"  # or grid_metrics.csv if you prefer

# ---- choose a group you want to inspect (from stacked_summary.csv top rows) ----
BCUT = 10
NBINS = 1080
CE = 2000

TARGET_LO = 140.0
TARGET_HI = 240.0

# ---- helpers ----
def est_period_zero_cross(y, x):
    y = np.asarray(y)
    zc = np.where(np.signbit(y[:-1]) != np.signbit(y[1:]))[0]
    if len(zc) < 4:
        return np.nan
    dx = np.diff(x[zc])
    return 2.0 * np.nanmean(dx)

def best_shift_corr(y, yref):
    """Return integer roll shift that maximizes correlation with yref."""
    y = np.asarray(y)
    yref = np.asarray(yref)
    y0 = y - y.mean()
    r0 = yref - yref.mean()
    # FFT-based circular cross-correlation
    fy = np.fft.rfft(y0)
    fr = np.fft.rfft(r0)
    cc = np.fft.irfft(fy * np.conj(fr), n=len(y0))
    # cc[k] corresponds to y rolled by k aligning to yref
    return int(np.argmax(cc))

def zero_crossings_interp(y, x):
    """Return interpolated x positions where y crosses 0."""
    y = np.asarray(y)
    x = np.asarray(x)
    idx = np.where(np.signbit(y[:-1]) != np.signbit(y[1:]))[0]
    xs = []
    for i in idx:
        y0, y1 = y[i], y[i+1]
        x0, x1 = x[i], x[i+1]
        if y1 == y0:
            continue
        t = -y0 / (y1 - y0)
        xs.append(x0 + t*(x1 - x0))
    return np.array(xs)

# ---- load metrics & select run tags ----
dfm = pd.read_csv(METR)
sel = dfm[(dfm["bcut"]==BCUT) & (dfm["nbins"]==NBINS) & (dfm["ce_trials"]==CE)].copy()
tags = sel["tag"].tolist()

print(f"Selected runs: {len(tags)} tags for bcut={BCUT}, nbins={NBINS}, ce={CE}")
print(tags)

# ---- load each run, pick IMF in target period band ----
runs = []
for tag in tags:
    d = BASE / tag
    sigp = d / "signal.csv"
    imfp = d / "imfs.npy"

    if not sigp.exists() or not imfp.exists():
        print("Missing files for", tag, "->", d)
        continue

    sig = pd.read_csv(sigp)
    x = sig["x_deg"].to_numpy()
    imfs = np.load(imfp)  # (n_imf, n)

    periods = np.array([est_period_zero_cross(imfs[k], x) for k in range(imfs.shape[0])])
    ok = np.where((periods >= TARGET_LO) & (periods <= TARGET_HI))[0]

    if len(ok) == 0:
        # fallback: use IMF6 if it exists else last IMF
        k = 5 if imfs.shape[0] >= 6 else (imfs.shape[0]-1)
    else:
        mid = 0.5*(TARGET_LO + TARGET_HI)
        k = ok[np.argmin(np.abs(periods[ok] - mid))]

    runs.append((tag, x, imfs[k], periods[k], k+1))

print("\nLoaded runs:", len(runs))
for tag, _, _, per, imf_idx in runs:
    print(f"  {tag} -> picked IMF{imf_idx} period~{per:.1f} deg")

# ---- align and stack ----
ref_tag, xref, yref, _, _ = runs[0]
aligned = []
shifts = []
for tag, x, y, per, imf_idx in runs:
    # x should match xref
    if len(x) != len(xref) or np.max(np.abs(x-xref)) > 1e-6:
        print("x grid mismatch for", tag, "skipping")
        continue
    s = best_shift_corr(y, yref)
    shifts.append(s)
    aligned.append(np.roll(y, s))

aligned = np.array(aligned)  # (n_runs, n)
stack_mean = aligned.mean(axis=0)
stack_std  = aligned.std(axis=0)

# ---- extract feature longitudes ----
i_max = int(np.argmax(stack_mean))
i_min = int(np.argmin(stack_mean))
l_peak   = float(xref[i_max])
l_trough = float(xref[i_min])

zc = zero_crossings_interp(stack_mean, xref)

def circ_dist(a, b):
    d = (a-b+180) % 360 - 180
    return abs(d)

print("\nSTACK FEATURES")
print("  ref:", ref_tag)
print("  mean shift applied (bins):", np.mean(shifts), "std:", np.std(shifts))
print("  l_peak  :", l_peak)
print("  l_trough:", l_trough)
print("  zero crossings (deg):", np.round(zc, 2)[:12], "..." if len(zc)>12 else "")

print("\nDistances to candidate directions:")
for L in [0.0, 180.0, 264.0]:
    print(f"  |l_peak-{L:>5.1f}| (circular) = {circ_dist(l_peak, L):.2f} deg")

# ---- plot ----
plt.figure(figsize=(10,4))
plt.plot(xref, stack_mean, label="stack mean")
plt.fill_between(xref, stack_mean-stack_std, stack_mean+stack_std, alpha=0.25, label="±1σ")
plt.axvline(l_peak, linestyle="--", linewidth=1, label=f"peak {l_peak:.1f}°")
plt.axvline(l_trough, linestyle="--", linewidth=1, label=f"trough {l_trough:.1f}°")
plt.title(f"Stacked IMF (bcut={BCUT}, nbins={NBINS}, ce={CE})")
plt.xlabel("galactic longitude (deg)")
plt.ylabel("IMF amplitude")
plt.legend()
plt.tight_layout()
plt.show()
