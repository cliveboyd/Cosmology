#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Joint Bayesian MCMC Fit: SN Ia + BAO + CMB with Accretion Decay
Model: Kerr-de Sitter Cosmology with Decaying Accretion Term
Author: Clive Stewart Boyd
Updated: 2025-07-25
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import emcee
import corner
import arviz as az
from scipy.integrate import quad

# === Load Data ===
sn_data = pd.read_csv("/Users/boyde/Downloads/lcparam_full_long (2).txt", sep=r'\s+', comment='#',
                      names=['CID', 'zCMB', 'zHD', 'zHDERR', 'MU', 'MUERR'])
z_sn, mu_obs, mu_err = sn_data['zCMB'].values, sn_data['MU'].values, sn_data['MUERR'].values
mask = z_sn > 0.005
z_sn, mu_obs, mu_err = z_sn[mask], mu_obs[mask], mu_err[mask]

bao = pd.read_csv("/Users/boyde/Downloads/eboss_dr16_bao_summary_corrected.csv")
z_bao  = bao["z"].values
DM_obs = bao["DM_over_rd"].values * 147.0
DM_err = bao["DM_err"].values * 147.0
H_obs  = bao["Hz_rd"].values / 147.0
H_err  = bao["Hz_err"].values / 147.0

# === Constants ===
c = 299792.458  # km/s
z_star = 1089.92
DA_star_obs = 14700.0  # Mpc
DA_star_err = 100.0
decay_power = 6  # Controls high-z suppression of A_acc

# === Model Functions ===
def A_acc_z(z, A_acc):
    return A_acc / (1 + z)**(decay_power - 2)

def E_z(z, Om, Ol, A_acc):
    return np.sqrt(Om*(1 + z)**3 + Ol + A_acc_z(z, A_acc)*(1 + z)**2)

def H_z(z, H0, Om, Ol, A_acc):
    return H0 * E_z(z, Om, Ol, A_acc)

def DM_z(z, H0, Om, Ol, A_acc):
    return (c / H0) * np.array([quad(lambda zp: 1 / E_z(zp, Om, Ol, A_acc), 0, zi)[0] for zi in z])

def mu_theory(z, H0, Om, Ol, A_acc):
    DL = (1 + z) * DM_z(z, H0, Om, Ol, A_acc)
    return 5 * np.log10(DL) + 25

def DA_z(z, H0, Om, Ol, A_acc):
    integral = quad(lambda zp: 1 / E_z(zp, Om, Ol, A_acc), 0, z)[0]
    return (c / H0) * integral / (1 + z)

# === Likelihood ===
def log_likelihood(theta):
    H0, Om, Ol, A_acc = theta
    mu_model = mu_theory(z_sn, H0, Om, Ol, A_acc)
    chi2_sn = np.sum(((mu_obs - mu_model) / mu_err) ** 2)

    DM_model = DM_z(z_bao, H0, Om, Ol, A_acc)
    H_model = H_z(z_bao, H0, Om, Ol, A_acc)
    chi2_DM = np.sum(((DM_obs - DM_model) / DM_err) ** 2)
    chi2_H  = np.sum(((H_obs - H_model) / H_err) ** 2)

    DA_model = DA_z(z_star, H0, Om, Ol, A_acc)
    chi2_CMB = ((DA_model - DA_star_obs) / DA_star_err) ** 2

    return -0.5 * (chi2_sn + chi2_DM + chi2_H + chi2_CMB)

# === Priors and Posterior ===
def log_prior(theta):
    H0, Om, Ol, A_acc = theta
    if 60 < H0 < 80 and 0.05 < Om < 0.5 and 0.6 < Ol < 0.95 and 0.0 <= A_acc < 0.2:
        return 0.0
    return -np.inf

def log_posterior(theta):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta)

# === MCMC Sampling ===
ndim, nwalkers = 4, 50
initial = np.array([72.0, 0.28, 0.72, 0.02])
pos = initial + 1e-3 * np.random.randn(nwalkers, ndim)
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior)

print("Running MCMC...")
sampler.run_mcmc(pos, 5000, progress=True)

# === Posterior Analysis ===
samples = sampler.get_chain(discard=1000, flat=True)
param_names = ["H0", "Omega_m", "Omega_L", "A_acc"]
df = pd.DataFrame(samples, columns=param_names)
df.to_csv("/Users/boyde/Downloads/posteriors_acc_decay.csv", index=False)

H0, Om, Ol, A_acc = np.median(samples, axis=0)
summary           = az.summary(samples, hdi_prob=0.68)
summary.to_csv("/Users/boyde/Downloads/summary_stats_acc_decay.csv")

# === Plots ===
corner.corner(samples, labels=param_names, quantiles=[0.16, 0.5, 0.84], show_titles=True)
plt.savefig("/Users/boyde/Downloads/corner_acc_decay.png")
plt.show()

# Residuals
mu_model = mu_theory(z_sn, H0, Om, Ol, A_acc)
DM_model = DM_z(z_bao, H0, Om, Ol, A_acc)
H_model  = H_z(z_bao, H0, Om, Ol, A_acc)

plt.figure(); plt.errorbar(z_sn, mu_obs - mu_model, yerr=mu_err, fmt='o'); plt.axhline(0, ls='--'); plt.title("μ(z) Residuals")
plt.savefig("/Users/boyde/Downloads/mu_residuals_decay.png")

plt.figure(); plt.errorbar(z_bao, DM_obs - DM_model, yerr=DM_err, fmt='o'); plt.axhline(0, ls='--'); plt.title("D_M(z) Residuals")
plt.savefig("/Users/boyde/Downloads/DM_residuals_decay.png")

plt.figure(); plt.errorbar(z_bao, H_obs - H_model, yerr=H_err, fmt='o'); plt.axhline(0, ls='--'); plt.title("H(z) Residuals")
plt.savefig("/Users/boyde/Downloads/H_residuals_decay.png")

# Energy density evolution
z_vals = np.linspace(0, 3, 300)
E2 = E_z(z_vals, Om, Ol, A_acc)**2
rho_m = Om * (1 + z_vals)**3 / E2
rho_L = Ol / E2
rho_acc = A_acc_z(z_vals, A_acc) * (1 + z_vals)**2 / E2

plt.figure(figsize=(8, 5))
plt.plot(z_vals, rho_m, label=r'$\Omega_m(z)$')
plt.plot(z_vals, rho_L, label=r'$\Omega_\Lambda(z)$')
plt.plot(z_vals, rho_acc, label=r'$\Omega_{\rm acc}(z)$')
plt.xlabel("Redshift z"); plt.ylabel("Fractional Energy Density")
plt.title("Relative Energy Contributions (with Accretion Decay)")
plt.legend(); plt.grid(); plt.tight_layout()
plt.savefig("/Users/boyde/Downloads/energy_budget_acc_decay.png")
plt.show()

# Metrics
def compute_metrics(logL, k, N):
    chi2 = -2 * logL
    dof = N - k
    return chi2 / dof, chi2 + 2*k, chi2 + k*np.log(N)

N_total = len(mu_obs) + len(DM_obs) + len(H_obs) + 1
logL = log_likelihood([H0, Om, Ol, A_acc])
red_chi2, aic, bic = compute_metrics(logL, ndim, N_total)
print(f"Reduced χ²: {red_chi2:.2f}, AIC: {aic:.2f}, BIC: {bic:.2f}")
print(f"DA(z=1089.92) = {DA_z(z_star, H0, Om, Ol, A_acc):.2f} Mpc (Planck: 14700 ± 100 Mpc)")

