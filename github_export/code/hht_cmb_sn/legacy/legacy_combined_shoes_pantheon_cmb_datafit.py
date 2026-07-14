#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Combined Cosmology Analysis:
- SH0ES+Pantheon+ Segment Chebyshev Fit
- Planck 2018 CMB Power Spectra (TT, TE, EE)
Author: Clive Stewart Boyd
"""

import numpy as np
import matplotlib.pyplot as plt
from   astropy.io import fits
from   numpy.polynomial.chebyshev import Chebyshev

# Load SH0ES+Pantheon+ Data
y_path = "/Users/boyde/Downloads/ally_shoes_ceph_topantheonwt6.0_112221.fits"
L_path = "/Users/boyde/Downloads/alll_shoes_ceph_topantheonwt6.0_112221.fits"
C_path = "/Users/boyde/Downloads/allc_shoes_ceph_topantheonwt6.0_112221.fits"

with fits.open(y_path) as hdul_y, fits.open(L_path) as hdul_L, fits.open(C_path) as hdul_C:
    y  = hdul_y[0].data.astype(np.float64)
    L  = hdul_L[0].data.astype(np.float64)
    C  = hdul_C[0].data.astype(np.float64)

W      = np.linalg.inv(C)
LW     = L @ W
x_best = np.linalg.inv(LW @ L.T) @ (LW @ y)
mu_fit = L.T @ x_best
mu_obs = y
mu_err = np.sqrt(np.diag(C))
z      = np.linspace(0.01, 2.3, len(mu_obs))

segments = [
    (1.414, 1.420), (1.420, 1.706), (1.706, 1.7105), (1.712, 1.747), (1.747, 2.059),
    (2.059, 2.063), (2.063, 2.113), (2.114, 2.119), (2.119, 2.124)
]

# Load Planck 2018 TT, TE, EE
tt_file = "/Users/boyde/Downloads/COM_PowerSpect_CMB-TT-full_R3.01.txt"
te_file = "/Users/boyde/Downloads/COM_PowerSpect_CMB-TE-full_R3.01.txt"
ee_file = "/Users/boyde/Downloads/COM_PowerSpect_CMB-EE-full_R3.01.txt"

def load_planck_data(filepath):
    data     = np.loadtxt(filepath, comments="#")
    ell      = data[:, 0]
    Dl       = data[:, 1]
    err_low  = data[:, 2]
    err_high = data[:, 3]
    err      = [np.abs(err_low), np.abs(err_high)]
    return ell, Dl, err

ell_tt, Dl_tt, err_tt = load_planck_data(tt_file)
ell_te, Dl_te, err_te = load_planck_data(te_file)
ell_ee, Dl_ee, err_ee = load_planck_data(ee_file)

# Plot Combined
fig, ax = plt.subplots(2, 1, figsize=(14, 10), gridspec_kw={'height_ratios': [2, 3]})

# Top: SH0ES+Pantheon+
ax[0].errorbar(z, mu_obs, yerr=mu_err, fmt='o', markersize=2, alpha=0.4, label="Observed μ")
colors = plt.cm.tab10.colors
for i, (z_min, z_max) in enumerate(segments):
    mask   = (z >= z_min) & (z <= z_max)
    z_seg  = z[mask]
    mu_seg = mu_obs[mask]
    if len(z_seg) < 6:
        continue
    z_norm   = 2 * (z_seg - z_seg.min()) / (z_seg.max() - z_seg.min()) - 1
    model    = Chebyshev.fit(z_norm, mu_seg, deg=6)
    mu_model = model(z_norm)
    ax[0].plot(z_seg, mu_model, linewidth=2, label=f"Segment {i+1}: z={z_min:.3f}-{z_max:.3f}", color=colors[i % len(colors)])
ax[0].set_xlabel("Redshift (z)")
ax[0].set_ylabel("Distance Modulus (μ)")
ax[0].set_title("Segmented Chebyshev Fits Overlayed on SH0ES+Pantheon+ Data")
ax[0].legend()
ax[0].grid(True)

# Bottom: Planck TT, TE, EE
ax[1].errorbar(ell_tt, Dl_tt, yerr=err_tt, fmt='.', alpha=0.6, label="Planck 2018 TT")
ax[1].errorbar(ell_te, Dl_te, yerr=err_te, fmt='.', alpha=0.6, label="Planck 2018 TE")
ax[1].errorbar(ell_ee, Dl_ee, yerr=err_ee, fmt='.', alpha=0.6, label="Planck 2018 EE")
ax[1].set_xlabel("Multipole ℓ")
ax[1].set_ylabel(r"$D_\ell = \ell(\ell+1)C_\ell / 2\pi$ [$\mu$K²]")
ax[1].set_title("Planck 2018 CMB Power Spectra (TT, TE, EE)")
ax[1].legend()
ax[1].grid(True)

plt.tight_layout()
plt.savefig("/Users/boyde/Downloads/combined_shoes_planck_overlay.png")
plt.show()
