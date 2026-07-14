#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 02:36:26 2025

Bayesian Hierarchical Model using SHOES Pantheon+ Residuals and Planck CMB TT/TE/EE

Author: Clive Stewart Boyd

Assumeptions...
- SHOES+Pantheon residuals are extracted from observed vs Chebyshev segmented fit
- Planck TT/TE/EE spectra loaded from official text files in Downloads
- Hierarchical priors link AM (angular momentum) evolution to metric perturbations

"""

import numpy as np
import matplotlib.pyplot as plt
from   astropy.io import fits
from   sklearn.preprocessing import StandardScaler
from   scipy.interpolate import interp1d
import pymc as pm
import arviz as az

# === Load SH0ES+Pantheon+ Data ===
y_path = "/Users/boyde/Downloads/ally_shoes_ceph_topantheonwt6.0_112221.fits"
L_path = "/Users/boyde/Downloads/alll_shoes_ceph_topantheonwt6.0_112221.fits"
C_path = "/Users/boyde/Downloads/allc_shoes_ceph_topantheonwt6.0_112221.fits"

with fits.open(y_path) as hdul_y, fits.open(L_path) as hdul_L, fits.open(C_path) as hdul_C:
    y = hdul_y[0].data.astype(np.float64)
    L = hdul_L[0].data.astype(np.float64)
    C = hdul_C[0].data.astype(np.float64)

W      = np.linalg.inv(C)
x_best = np.linalg.inv(L @ W @ L.T) @ (L @ W @ y)
mu_fit = L.T @ x_best
mu_obs = y
mu_err = np.sqrt(np.diag(C))

z      = np.linspace(0.01, 2.3, len(mu_obs))
z_norm = 2 * (z - z.min()) / (z.max() - z.min()) - 1

# === Load Planck TT, TE, EE Spectra ===
tt_data = np.loadtxt("/Users/boyde/Downloads/COM_PowerSpect_CMB-TT-full_R3.01.txt")
te_data = np.loadtxt("/Users/boyde/Downloads/COM_PowerSpect_CMB-TE-full_R3.01.txt")
ee_data = np.loadtxt("/Users/boyde/Downloads/COM_PowerSpect_CMB-EE-full_R3.01.txt")

Dl_tt = tt_data[:, 1]
Dl_te = te_data[:, 1]
Dl_ee = ee_data[:, 1]

# === Normalize and Resample ===
scaler = StandardScaler()
Dl_tt_norm = scaler.fit_transform(Dl_tt.reshape(-1, 1)).flatten()
Dl_te_norm = scaler.fit_transform(Dl_te.reshape(-1, 1)).flatten()
Dl_ee_norm = scaler.fit_transform(Dl_ee.reshape(-1, 1)).flatten()

N_z = len(z_norm)
x_target = np.linspace(0, 1, N_z)

interp_tt = interp1d(np.linspace(0, 1, len(Dl_tt_norm)), Dl_tt_norm)
interp_te = interp1d(np.linspace(0, 1, len(Dl_te_norm)), Dl_te_norm)
interp_ee = interp1d(np.linspace(0, 1, len(Dl_ee_norm)), Dl_ee_norm)

Dl_tt_resampled = interp_tt(x_target)
Dl_te_resampled = interp_te(x_target)
Dl_ee_resampled = interp_ee(x_target)

cmb_driver = (Dl_tt_resampled + Dl_te_resampled + Dl_ee_resampled) / 3.0

# === Residuals ===
residuals = mu_obs - mu_fit

# === Hierarchical Bayesian Model ===
with pm.Model() as hierarchical_model:
    mu_alpha  = pm.Normal("mu_alpha", mu=0, sigma=1)
    mu_beta   = pm.Normal("mu_beta", mu=0, sigma=1)

    alpha     = pm.HalfNormal("alpha", sigma=0.1)
    beta      = pm.HalfNormal("beta", sigma=0.1)

    mu        = mu_alpha + alpha * cmb_driver + beta * z_norm
    sigma_obs = pm.HalfNormal("sigma_obs", sigma=1)

    obs       = pm.Normal("obs", mu=mu, sigma=sigma_obs, observed=residuals)

    trace     = pm.sample(1000, tune=1000, target_accept=0.9, cores=4, random_seed=42)

# === Posterior Plot ===
az.plot_trace(trace, var_names=["mu_alpha", "mu_beta", "alpha", "beta", "sigma_obs"])
plt.tight_layout()
plt.savefig("/Users/boyde/Downloads/posterior_hierarchical_am_model.png")
plt.show()

