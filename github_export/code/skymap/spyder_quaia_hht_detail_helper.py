#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  8 07:11:05 2025

@author: boyde

quaia_hht_detail_helper.py

"""

import os
import numpy as np
import matplotlib.pyplot as plt

# --- config ---
bcut     = 25
base_out = "/Users/boyde/.spyder-py3/plamb_runs/hht_quaia"
base_in  = "/Users/boyde/.spyder-py3/skymap2/quaia_outputs"

detail_npz = os.path.join(base_out, f"bcut{bcut}_allz_apar_signed", "hht_quaia_detail.npz")
stats_npz  = os.path.join(base_out, f"bcut{bcut}_allz_apar_signed", "hht_quaia_imf_stats.npz")
series_txt = os.path.join(base_in,  f"quaia_hht_bcut{bcut}_allz_apar_signed.txt")

# --- load HHT detail ---
d         = np.load(detail_npz)
Ls        = d["Ls"]        # log(1+z) grid
IMFs      = d["IMFs"]      # shape (n_imf, n_grid)
z_raw     = d["z_raw"]     # original z_mid used for the series
y_raw     = d["y_raw"]     # original signed A_par(z_mid)

# convert L -> z grid
z_grid    = np.expm1(Ls)

# --- load IMF stats and pick low-frequency IMF ---
s         = np.load(stats_npz, allow_pickle=True)
imf_stats = s["imf_stats"]   # dtype: ('IMF','energy','f_med','rank')

names     = imf_stats.dtype.names

if ("rank" in names) and ("IMF" in names):
    idx = np.where(imf_stats["rank"] == 1)[0]
    
    if idx.size > 0:
        imf_idx    = int(imf_stats["IMF"][idx[0]])  # 1-based index
        
    else:
        # fall back to smallest f_med
        imf_idx    = int(imf_stats["IMF"][np.argmin(imf_stats["f_med"])])
else:
    # very defensive: smallest f_med if no rank field
    imf_idx        = int(imf_stats["IMF"][np.argmin(imf_stats["f_med"])])

print(f"[select] using IMF {imf_idx} as low-frequency mode")

lowfreq             = IMFs[imf_idx - 1, :]   # 0-based index in the array

# --- load the signed A_par(z_mid) series explicitly (same as z_raw/y_raw) ---
z_mid, A_par_signed = np.loadtxt(series_txt, unpack=True)

# --- fit scale + offset to match IMF shape to the binned A_par(z) ---
lf_at_zmid          = np.interp(z_mid, z_grid, lowfreq)

A_mat               = np.vstack([lf_at_zmid, np.ones_like(lf_at_zmid)]).T
coeffs, *_          = np.linalg.lstsq(A_mat, A_par_signed, rcond=None)
scale, offset       = coeffs
lowfreq_rescaled    = scale * lowfreq + offset

print(f"[fit] scale={scale:.3g}, offset={offset:.3g}")

# --- make diagnostic plot ---
plt.figure(figsize      = (9,5))
plt.axhline(0.0, color  = "grey", linestyle="--", linewidth=0.8)

plt.scatter(z_mid, A_par_signed,
            s=40, alpha = 0.85,
            label       = r"$A_{\parallel}(z)$ signed (z-bin centres)")

plt.plot(z_grid, lowfreq_rescaled,
         linewidth      = 2.0,
         label          = f"Low-freq IMF {imf_idx} (scaled to match)")

plt.xlabel("z")
plt.ylabel(r"$A_{\parallel}$ (CMB-aligned, signed)")
plt.title(f"|b| > {bcut}° : signed $A_{{\\parallel}}(z)$ and low-frequency IMF")
plt.legend()
plt.tight_layout()

out_fig = os.path.join(base_out,
                       f"bcut{bcut}_allz_apar_signed",
                       "hht_quaia_lowfreq_signed_overlay.png")
plt.savefig(out_fig, dpi=150)
print(f"[plot] saved {out_fig}")
plt.show()
