#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Bayesian Fit for Kerr–de Sitter Cosmology with Rotation and Accretion
Fits to: Pantheon+SH0ES SN Ia, BAO distance proxies
Author: Clive Stewart Boyd
"""

import numpy as np
import matplotlib.pyplot as plt
import pymc as pm
import arviz as az
from astropy.io import fits

# === Load SN Ia Data ===
y_path = "/Users/boyde/Downloads/ally_shoes_ceph_topantheonwt6.0_112221.fits"
L_path = "/Users/boyde/Downloads/alll_shoes_ceph_topantheonwt6.0_112221.fits"
C_path = "/Users/boyde/Downloads/allc_shoes_ceph_topantheonwt6.0_112221.fits"

with fits.open(y_path) as hdul_y, fits.open(L_path) as hdul_L, fits.open(C_path) as hdul_C:
    y = hdul_y[0].data.astype(np.float64)
    L = hdul_L[0].data.astype(np.float64)
    C = hdul_C[0].data.astype(np.float64)

W = np.linalg.inv(C)
LW = L @ W
x_best = np.linalg.inv(LW @ L.T) @ (LW @ y)
mu_obs = y
mu_err = np.sqrt(np.diag(C))
z_data = np.linspace(0.01, 2.3, len(mu_obs))

# === BAO proxy data ===
bao_z = np.array([0.38, 0.51, 0.61])
bao_dm = np.array([1512.39, 1975.22, 2306.68])
bao_err = np.array([24.0, 29.0, 37.0])

c = 299792.458  # speed of light in km/s

# === Theoretical Functions ===
def kerr_de_sitter_mu(z, H0, Om, OL, A_acc):
    def E(z):
        acc = A_acc * (1 - np.exp(-z / 0.5))
        return np.sqrt(Om * (1 + z) ** 3 + OL + acc)

    z = np.atleast_1d(z)
    integral = []
    for zi in z:
        if zi <= 0:
            integral.append(0.0)
        else:
            zz = np.linspace(0.001, zi, 200)
            Ez = E(zz)
            if Ez.shape != zz.shape:
                if Ez.size == 1:
                    Ez = np.full_like(zz, Ez.item())
                else:
                    raise ValueError(f"Shape mismatch: E(z)={Ez.shape}, z={zz.shape}")
            integral.append(np.trapz(1.0 / Ez, zz))
    dL = (c / H0) * (1 + z) * np.array(integral)
    return 5 * np.log10(dL) + 25

def dM_bao(z, H0, Om, OL, A_acc):
    def E(z):
        acc = A_acc * (1 - np.exp(-z / 0.5))
        return np.sqrt(Om * (1 + z) ** 3 + OL + acc)

    z = np.atleast_1d(z)
    integral = []
    for zi in z:
        if zi <= 0:
            integral.append(0.0)
        else:
            zz = np.linspace(0.001, zi, 200)
            Ez = E(zz)
            if Ez.shape != zz.shape:
                if Ez.size == 1:
                    Ez = np.full_like(zz, Ez.item())
                else:
                    raise ValueError(f"Shape mismatch in BAO E(z): {Ez.shape} vs {zz.shape}")
            integral.append(np.trapz(1.0 / Ez, zz))
    return (c / H0) * np.array(integral)

# === Bayesian Model ===
with pm.Model() as model:
    H0 = pm.Uniform("H0", lower=60, upper=80)
    Om = pm.Uniform("Omega_m", lower=0.1, upper=0.5)
    A_acc = pm.HalfNormal("A_acc", sigma=0.3)
    OL = pm.Deterministic("Omega_L", 1 - Om)

    def loglike_fn(H0_val, Om_val, OL_val, A_acc_val):
        try:
            mu_theory = kerr_de_sitter_mu(z_data, H0_val, Om_val, OL_val, A_acc_val)
            bao_theory = dM_bao(bao_z, H0_val, Om_val, OL_val, A_acc_val)
            sn_term = -0.5 * np.sum(((mu_obs - mu_theory) / mu_err) ** 2)
            bao_term = -0.5 * np.sum(((bao_dm - bao_theory) / bao_err) ** 2)
            return sn_term + bao_term
        except Exception as e:
            raise RuntimeError(f"Error in log-likelihood calculation: {e}")

    pm.DensityDist("likelihood", logp=lambda H0=H0, Om=Om, OL=OL, A_acc=A_acc: loglike_fn(H0, Om, OL, A_acc))
    trace = pm.sample(1000, tune=1000, chains=2, target_accept=0.95, return_inferencedata=True)

# === Posterior Summary ===
az.plot_trace(trace)
plt.tight_layout()
plt.show()

summary = az.summary(trace, hdi_prob=0.68)
print(summary)