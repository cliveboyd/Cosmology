#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PLamb â€“ Pantheon(+SH0ES) model fitting â€” V10V6c_plus
Adds BAO, CC (cosmic chronometers), and Planck distance-prior (l_A, R) options.

SU(2) TEST-BRANCH NOTE
---------------------
This Windows-migration version has been extended with a first-pass,
phenomenological chiral SU(2)-inspired dark-energy test branch.

This is NOT yet a derived Yang-Mills field-equation solver. It is the safe
first diagnostic described in the SU2 proposal: add an effective density term
Omega_chi(z) to the Friedmann expansion and test whether the existing
Pantheon+SH0ES / BAO / Planck likelihood machinery can tolerate or prefer it.

New model:
  --model SU2
  --model SU2R

New SU(2)-effective parameters:
  --Omega_chi0 FLOAT   Present-day density fraction for the SU(2)-like component.
  --su2-fraction FLOAT SU2R fraction of the dark sector assigned to Omega_chi0.
  --Omega_dark0 FLOAT  SU2R total dark-sector density before the Lambda/SU2 split.
  --w0-chi FLOAT       Present-day equation of state for Omega_chi(z).
  --wa-chi FLOAT       CPL-style redshift evolution parameter for Omega_chi(z).

Implemented form:
  w_chi(z) = w0_chi + wa_chi * z / (1 + z)
  Omega_chi(z) = Omega_chi0 * a^[-3(1+w0_chi+wa_chi)] * exp[-3 wa_chi (1-a)]
  where a = 1 / (1 + z).

Flat-mode closure:
  When --flat is used, the code enforces present-day closure as:
      Omega_Lambda = 1 - Omega_m - Omega_chi0
  so the SU(2)-like term replaces part of the ordinary Lambda budget rather
  than being added on top of a flat LCDM universe.
  In --model SU2R, flat closure instead samples su2_fraction directly:
      Omega_dark0 = 1 - Omega_m
      Omega_chi0  = su2_fraction * Omega_dark0
      Omega_Lambda = (1 - su2_fraction) * Omega_dark0

Model hygiene:
  For LCDM, KdS, and FR runs, Omega_chi0 is automatically frozen to zero and
  w0_chi / wa_chi are frozen to their neutral values. This preserves the old
  baseline behavior.

First validation script:
  C:/Users/clive/Documents/Cosmology/test_run_plamb_su2_spyder.py

Recommended first SU2 test:
  Use --sampler none with --model SU2, --flat, --Omega_chi0 0.05,
  --w0-chi -1.0, --wa-chi 0.0. This should behave like a Lambda split and
  therefore reproduce the LCDM chi2 while showing the AIC/BIC penalty for
  the extra active parameters.

MODEL GUIDE
-----------
Available model branches:

  --model LCDM
      Baseline Lambda-CDM distance model.
      Automatically freezes:
        A_acc = 0, n_acc = 0, gamma_c = 0, epsilon_M = 0,
        Omega_chi0 = 0, w0_chi = -1, wa_chi = 0.
      Use this as the control run for every new data or code change.

  --model KdS
      Kerr-de Sitter / accretion-style branch.
      Allows A_acc and n_acc to act as phenomenological acceleration terms.
      Automatically freezes propagation/magnitude-evolution terms:
        gamma_c = 0, epsilon_M = 0.
      Use this to compare the older rotating/accretion hypothesis against
      LCDM and SU2 on the same datasets.

  --model FR
      Full Relativity exploratory wrapper.
      Allows A_acc, n_acc, gamma_c, and epsilon_M to be active unless fixed.
      This is the broadest existing non-SU2 model and is useful for checking
      whether SU2 is really adding evidence or merely duplicating freedom that
      FR already had.

  --model SU2
      First-pass chiral SU(2)-inspired effective dark-energy branch.
      Automatically freezes the older FR/KdS operator terms:
        A_acc = 0, n_acc = 0, gamma_c = 0, epsilon_M = 0.
      Allows Omega_chi0, w0_chi, and wa_chi to be active unless fixed.
      In --flat mode, Omega_chi0 replaces part of Omega_Lambda.

  --model SU2R
      Reparameterized SU2 branch with the same physical likelihood as SU2.
      In flat mode it samples H0, Om, su2_fraction, and w0_chi rather than
      H0, Om, Omega_chi0, and w0_chi. This directly tests the dark-sector
      split exposed by the SU2 chain audit.

OPERATOR / PARAMETER GUIDE
--------------------------
Core expansion parameters:
  H0          Hubble constant in km/s/Mpc.
  Om          Present-day matter density fraction Omega_m.
  Ol          Present-day Lambda density fraction Omega_Lambda.
  Ok          Derived curvature term, Ok = 1 - Om - Ol - Omega_chi0.

KdS / FR phenomenological operators:
  A_acc       Amplitude of the accretion/extra-acceleration term in E(z)^2.
  n_acc       Redshift power of that term: A_acc * (1+z)^n_acc.
  gamma_c     Linear variable-c coefficient used by c_of_z(z).
  epsilon_M   Supernova magnitude-evolution coefficient.

SU2 effective-sector operators:
  Omega_chi0  Present-day density fraction of the SU2-like component.
  Omega_dark0 Total late-time dark-sector density used by SU2R.
  su2_fraction Fraction of Omega_dark0 assigned to the SU2-like component.
  w0_chi      Present-day equation of state for the SU2-like component.
  wa_chi      CPL redshift-evolution term for w_chi(z).

Probe operators / likelihood blocks:
  SN          Pantheon+SH0ES distance-modulus likelihood.
  BAO         BAO distance-scale likelihood using DM, DH, or DV rows.
  CC          Cosmic chronometer H(z) likelihood if --cc is supplied.
  Planck      Distance-prior likelihood using l_A and R.

Sampler modes:
  --sampler none
      Fast single-point evaluation. Use first for migration and smoke tests.
  --sampler grid
      One-parameter scan. Freeze all but one active parameter.
  --sampler emcee
      Full MCMC. Requires emcee; use only after sampler-none tests pass.

SPYDER SAMPLE SCRIPTS
---------------------
Open and run these directly in Spyder:
  test_run_plamb_lcdm_spyder.py
  test_run_plamb_kds_spyder.py
  test_run_plamb_fr_spyder.py
  test_run_plamb_su2_spyder.py

New flags (all optional):
  --bao PATH            BAO table (CSV or whitespace). Columns: z,kind,value[,sigma]
  --bao-cov PATH        Covariance matrix for BAO (NxN text/CSV). If provided, 'sigma' is optional in table.
  --rd FLOAT            Sound horizon at drag epoch r_d [Mpc] for BAO (default 147.09)

  --cc PATH             Cosmic-chronometer H(z) table (CSV/whitespace). Columns: z,H,sigma

  --planck PATH         Distance prior file for (l_A, R). JSON or NPZ with 'means' and 'cov' (2x2).
  --planck-zstar FLOAT  Redshift of decoupling z_* (default 1089.0)
  --rs FLOAT            Sound horizon r_s used in l_A (defaults to --rd)
  --planck-nint INT     Fixed Simpson panels for z~z_* integral (default 4000)

Per-probe Ï‡Â² contributions are added to the SN likelihood and reported.

PLAMB Noether/conservation scaffold:
  --plamb-conservation-mode free
      Preserve the current phenomenological PLAMB behavior.
  --plamb-conservation-mode Gc_const
      Apply a small-redshift diagnostic constraint motivated by G/c constant:
      epsilon_M = (2.5 / ln(10)) * gamma_c.
  --plamb-conservation-mode Gc2_const
      Apply a small-redshift diagnostic constraint motivated by G/c^2 constant:
      epsilon_M = (5.0 / ln(10)) * gamma_c.

These modes are not a proof of Noether consistency. They are a conservative
way to turn Peter Lamb's conservation hypothesis into an explicit fit option.
"""

from __future__ import annotations

import argparse
import os
import time
import warnings
import json
import numpy as np

try:
    import emcee
    
except Exception:
    emcee       = None

try:
    import matplotlib.pyplot as plt
    
except Exception:
    plt         = None

try:
    import corner
    
except Exception:
    corner      = None

# ---------------- constants ----------------
C_KMS           = 299_792.458                      # km/s
LOG10_E         = np.log10(np.e)
LN10            = np.log(10.0)
PLAMB_MAG_PER_LOG10 = 2.5


def plamb_conservation_power(mode: str) -> float:
    m = (mode or 'free').lower()
    if m == 'gc_const':
        return 1.0
    if m == 'gc2_const':
        return 2.0
    return 0.0


def plamb_constrained_epsilon(gamma_c: float, mode: str) -> float:
    power = plamb_conservation_power(mode)
    if power == 0.0:
        return float(gamma_c) * 0.0
    return float(PLAMB_MAG_PER_LOG10 * power * float(gamma_c) / LN10)


def apply_plamb_conservation(full: dict, mode: str) -> dict:
    if plamb_conservation_power(mode) == 0.0:
        return full
    full['epsilon_M'] = plamb_constrained_epsilon(full.get('gamma_c', 0.0), mode)
    return full


def plamb_conservation_summary(mode: str, gamma_c: float, epsilon_M: float, z_ref: float = 1.0) -> dict:
    power = plamb_conservation_power(mode)
    c_ratio = float(max(1e-12, 1.0 + float(gamma_c) * float(z_ref)))
    out = {
        'plamb_conservation_mode': mode,
        'plamb_conservation_power_p': float(power),
        'plamb_constrained_epsilon_M': float(plamb_constrained_epsilon(gamma_c, mode)) if power else float(epsilon_M),
        'plamb_epsilon_M_relation': 'free' if power == 0.0 else f'epsilon_M = {PLAMB_MAG_PER_LOG10 * power / LN10:.8g} * gamma_c',
        'c_ratio_z1': c_ratio,
    }
    if power:
        out.update({
            'G_ratio_z1_for_mode': float(c_ratio ** power),
            'inertia_proxy_ratio_z1': float(c_ratio ** (-power)),
            'conservation_hypothesis': f'Noether diagnostic scaffold: assumes G/c^{power:g} constant and ties epsilon_M to gamma_c; not a full action-level derivation',
        })
    else:
        out.update({
            'G_ratio_z1_for_mode': float('nan'),
            'inertia_proxy_ratio_z1': float('nan'),
            'conservation_hypothesis': 'free PLAMB fit; Noether consistency not constrained',
        })
    return out


# ---------------- utils --------------------
def log10(x):
    return np.log(x) * LOG10_E


def ensure_dir(path: str):
    if not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)


def format_model_tag(model, best, k_active, chi2, dof, AIC, BIC, frozen_items):
    H0, Om, Ol, A_acc, n_acc, gamma_c, epsilon_M, Omega_chi0, w0_chi, wa_chi, *extra = best
    plamb_alpha = extra[0] if extra else 0.5
    
    frz_parts = [f"{k}={v:g}" for k, v in frozen_items.items()] if frozen_items else []
    frz_txt   = (" | frozen: " + ", ".join(frz_parts)) if frz_parts else ""
    metrics   = f"Ï‡Â²/dof={chi2/dof:.3f}  (k={k_active}, AIC={AIC:.1f}, BIC={BIC:.1f})"
    if model.upper() == 'PLAMB':
        params = (f"H0={H0:.2f}, Omega_inertia={Om:.3f}, "
                  f"clump_amp={A_acc:.3f}, clump_index={n_acc:.3f}, "
                  f"c_drift={gamma_c:.4f}, inertia_drift={epsilon_M:.3f}, "
                  f"alpha={plamb_alpha:.3f}")
    else:
        params = (f"H0={H0:.2f}, Î©m={Om:.3f}, Î©Î›={Ol:.3f}, "
                  f"A_acc={A_acc:.3f}, n_acc={n_acc:.3f}, Î³_c={gamma_c:.3f}, Îµ_M={epsilon_M:.3f}, "
                  f"Omega_chi0={Omega_chi0:.3f}, w0_chi={w0_chi:.3f}, wa_chi={wa_chi:.3f}")
    
    return f"{model} | {params}\n{metrics}{frz_txt}"


# --------- BH interpretation helper --------
def w_eff_from_nacc(n_acc: float) -> float:
    return (n_acc / 3.0) - 1.0


# --------------- data load (SNe) -----------
def _find_col_case_insensitive(names, target):
    t               = target.lower()
    
    for n in names:
        if n.lower() == t:
            return n
        
    return None


def load_pantheon(data_path: str):
    arr             = np.genfromtxt(data_path, names=True, dtype=None, encoding=None)
    names           = arr.dtype.names

    z_name          = _find_col_case_insensitive(names, 'zCMB') or _find_col_case_insensitive(names, 'zcmb')
    mu_name         = _find_col_case_insensitive(names, 'MU_SH0ES') or _find_col_case_insensitive(names, 'mu_sh0es')
    er_name         = _find_col_case_insensitive(names, 'MU_SH0ES_ERR_DIAG') or _find_col_case_insensitive(names, 'mu_sh0es_err_diag')

    if not (z_name and mu_name and er_name):
        raise ValueError(f"Could not find required columns in {data_path}. Found: {names}")

    z               = np.asarray(arr[z_name],  dtype=float)
    mu              = np.asarray(arr[mu_name], dtype=float)
    sigma           = np.asarray(arr[er_name], dtype=float)

    return {'z': z, 'mu': mu, 'sigma': sigma, 'N': z.size}


def load_covariance(cov_path: str, N: int) -> np.ndarray | None:
    if not cov_path:
        return None
    
    if not os.path.exists(cov_path):
        warnings.warn(f"Covariance file not found: {cov_path}")
    
        return None

    try:
        C = np.loadtxt(cov_path)
        
    except Exception as e:
        warnings.warn(f"Failed to read covariance file as ASCII text: {e}")
        return None

    if C.ndim == 2 and C.shape == (N, N):
        return C

    flat        = C.ravel()
    NN          = N * N
    if flat.size == NN:
        return flat.reshape(N, N)
    
    if flat.size > NN:
        return flat[-NN:].reshape(N, N)

    warnings.warn(f"Covariance has size {flat.size}, cannot reshape to ({N},{N}); ignoring.")
    return None


# -------------- cosmology ------------------
def c_of_z(z, gamma_c):
    cz          = C_KMS * (1.0 + gamma_c * z)
    
    return np.maximum(cz, 1.0)


def deltaM_of_z(z, epsilon_M, mode='linear_z'):
    if epsilon_M == 0.0:
        return 0.0 * z
    
    if mode == 'log1pz':
        return epsilon_M * np.log1p(z)
    
    return epsilon_M * z


def omega_chi_z(z, Omega_chi0, w0_chi, wa_chi):
    z           = np.asarray(z, dtype=float)
    
    if Omega_chi0 == 0.0:
        return np.zeros_like(z, dtype=float)
    
    a           = 1.0 / (1.0 + z)
    exponent    = -3.0 * wa_chi * (1.0 - a)
    
    return Omega_chi0 * a ** (-3.0 * (1.0 + w0_chi + wa_chi)) * np.exp(exponent)


def E_of_z(z, Om, Ol, A_acc, n_acc, Omega_chi0=0.0, w0_chi=-1.0, wa_chi=0.0):
    zp1         = 1.0 + z
    Ok          = 1.0 - Om - Ol - Omega_chi0
    E2          = (Om * zp1**3) + (Ol) + (Ok * zp1**2) + (A_acc * zp1**(n_acc))
    E2         += omega_chi_z(z, Omega_chi0, w0_chi, wa_chi)
    
    return np.sqrt(np.clip(E2, 1e-12, None))


def _simpson_integral_0z(zi, nint, H0, Om, Ol, A_acc, n_acc, gamma_c,
                         Omega_chi0=0.0, w0_chi=-1.0, wa_chi=0.0):
    
    grid        = np.linspace(0.0, zi, nint + 1)
    Ez          = E_of_z(grid, Om, Ol, A_acc, n_acc, Omega_chi0, w0_chi, wa_chi)
    cz_ratio    = c_of_z(grid, gamma_c) / C_KMS
    integrand   = cz_ratio / Ez

    h           = (zi - 0.0) / nint
    S           = integrand[0] + integrand[-1] + 4.0 * np.sum(integrand[1:-1:2]) + 2.0 * np.sum(integrand[2:-1:2])
    
    return (h / 3.0) * S


def comoving_distance(z, H0, Om, Ol, A_acc, n_acc, gamma_c, nint_base=200,
                      Omega_chi0=0.0, w0_chi=-1.0, wa_chi=0.0):
    
    if np.ndim(z) == 0:
        z       = np.array([float(z)], dtype=float)
    
    z           = np.asarray(z, dtype=float)

    D_H         = C_KMS / H0
    dc          = np.empty_like(z)

    for i, zi in enumerate(z):
        nint    = int(max(50, int(nint_base * (1.0 + zi))))
        nint   += (nint % 2)
        I0z     = _simpson_integral_0z(zi, nint, H0, Om, Ol, A_acc, n_acc, gamma_c,
                                     Omega_chi0, w0_chi, wa_chi)
        dc[i]   = D_H * I0z

    return dc


def comoving_distance_fixed_panels(zi, H0, Om, Ol, A_acc, n_acc, gamma_c, n_panels=4000,
                                   Omega_chi0=0.0, w0_chi=-1.0, wa_chi=0.0):
    """Comoving distance to scalar z=zi using a fixed Simpson panel count (for z~1100 speed)."""
    D_H         = C_KMS / H0
    nint        = int(max(50, int(n_panels)))
    nint       += (nint % 2)
    I0z         = _simpson_integral_0z(float(zi), nint, H0, Om, Ol, A_acc, n_acc, gamma_c,
                                       Omega_chi0, w0_chi, wa_chi)
    return D_H * I0z


def _safe_sinh(x):
    x_clamp     = np.clip(x, -20.0, 20.0)
    
    return np.sinh(x_clamp)


def transverse_comoving_distance(DC, H0, Om, Ol, Omega_chi0=0.0):
    Ok          = 1.0 - Om - Ol - Omega_chi0
    
    if np.allclose(Ok, 0.0, atol=1e-10):
        return DC

    sqrtOk      = np.sqrt(np.abs(Ok))
    DH          = C_KMS / H0
    x           = DC * sqrtOk / DH

    if Ok > 0:
        return DH / sqrtOk * _safe_sinh(x)
    
    else:
        x_clamp = np.clip(x, -50.0, 50.0)
        
        return DH / sqrtOk * np.sin(x_clamp)


def luminosity_distance(z, H0, Om, Ol, A_acc, n_acc, gamma_c, nint_base=200,
                        Omega_chi0=0.0, w0_chi=-1.0, wa_chi=0.0):
    
    DC          = comoving_distance(z, H0, Om, Ol, A_acc, n_acc, gamma_c, nint_base=nint_base,
                                    Omega_chi0=Omega_chi0, w0_chi=w0_chi, wa_chi=wa_chi)
    
    if not np.all(np.isfinite(DC)):
        return np.full_like(z, np.inf)
    
    DM          = transverse_comoving_distance(DC, H0, Om, Ol, Omega_chi0)
    
    if not np.all(np.isfinite(DM)):
        return np.full_like(z, np.inf)
    
    DL          = (1.0 + np.asarray(z)) * DM
    DL          = np.where(DL > 0, DL, np.inf)
    
    return DL


def mu_model(z, H0, Om, Ol, A_acc, n_acc, gamma_c, epsilon_M, dM_mode='linear_z', nint_base=200,
             Omega_chi0=0.0, w0_chi=-1.0, wa_chi=0.0):
    
    DL          = luminosity_distance(z, H0, Om, Ol, A_acc, n_acc, gamma_c, nint_base=nint_base,
                                      Omega_chi0=Omega_chi0, w0_chi=w0_chi, wa_chi=wa_chi)
    
    if not np.all(np.isfinite(DL)):
        return np.full_like(z, np.inf)
    
    mu_th       = 5.0 * log10(np.maximum(DL, 1e-12)) + 25.0
    mu_th      += deltaM_of_z(np.asarray(z), epsilon_M, mode=dM_mode)
    
    return mu_th


def _plamb_clump_kernel(z, A_clump, n_clump, c_drift):
    """No-expansion pilot kernel for redshift as clumping/inertia accumulation.

    This is a phenomenological test scaffold for Peter Lamb's note, not a
    completed action-level theory. It uses log(1+z) as the static redshift path,
    lets clumping/inertia perturb that path, and lets c vary mildly along it.
    """
    z           = np.asarray(z, dtype=float)
    x           = np.log1p(np.maximum(z, 0.0))
    n_eff       = np.clip(float(n_clump), 0.05, 4.0)
    growth      = np.power(np.maximum(x, 0.0), n_eff)
    clump       = np.exp(np.clip(float(A_clump) * growth, -8.0, 8.0))
    c_ratio     = np.clip(c_of_z(z, c_drift) / C_KMS, 0.05, 20.0)
    return c_ratio * clump / np.maximum(1.0 + z, 1e-12)


def _plamb_integral_0z(zi, nint, A_clump, n_clump, c_drift):
    grid        = np.linspace(0.0, zi, nint + 1)
    integrand   = _plamb_clump_kernel(grid, A_clump, n_clump, c_drift)
    h           = zi / nint
    S           = integrand[0] + integrand[-1] + 4.0 * np.sum(integrand[1:-1:2]) + 2.0 * np.sum(integrand[2:-1:2])
    return (h / 3.0) * S


def plamb_noexp_distances(z, H0, A_clump, n_clump, c_drift, nint_base=200, dimming_alpha=0.5):
    """Return effective (DM, DH, DV, DL) distances for the PLAMB pilot branch."""
    scalar      = np.ndim(z) == 0
    z_arr       = np.array([float(z)], dtype=float) if scalar else np.asarray(z, dtype=float)
    D_H         = C_KMS / H0
    DM          = np.empty_like(z_arr)
    DH_eff      = np.empty_like(z_arr)

    for i, zi in enumerate(z_arr):
        nint    = int(max(50, int(nint_base * (1.0 + min(float(zi), 10.0)))))
        nint   += (nint % 2)
        I0z     = _plamb_integral_0z(float(zi), nint, A_clump, n_clump, c_drift)
        DM[i]   = D_H * I0z
        DH_eff[i] = D_H * _plamb_clump_kernel(np.array([float(zi)]), A_clump, n_clump, c_drift)[0]

    DV          = np.power(np.maximum(z_arr * DM * DM * DH_eff, 1e-30), 1.0 / 3.0)
    DL          = np.power(np.maximum(1.0 + z_arr, 1e-12), float(dimming_alpha)) * DM
    if scalar:
        return DM[0], DH_eff[0], DV[0], DL[0]
    return DM, DH_eff, DV, DL


def plamb_noexp_distance_fixed_panels(zi, H0, A_clump, n_clump, c_drift, n_panels=4000):
    D_H         = C_KMS / H0
    nint        = int(max(50, int(n_panels)))
    nint       += (nint % 2)
    return D_H * _plamb_integral_0z(float(zi), nint, A_clump, n_clump, c_drift)


def mu_model_plamb_noexp(z, H0, A_clump, n_clump, c_drift, inertia_drift,
                         dimming_alpha=0.5,
                         dM_mode='linear_z', nint_base=200):
    _, _, _, DL = plamb_noexp_distances(z, H0, A_clump, n_clump, c_drift,
                                        nint_base=nint_base, dimming_alpha=dimming_alpha)
    if not np.all(np.isfinite(DL)):
        return np.full_like(np.asarray(z, dtype=float), np.inf)
    mu_th       = 5.0 * log10(np.maximum(DL, 1e-12)) + 25.0
    mu_th      += deltaM_of_z(np.asarray(z, dtype=float), inertia_drift, mode=dM_mode)
    return mu_th


# -------------- dispatcher -----------------
def model_mu(z, model_name, H0, Om, Ol, A_acc, n_acc, gamma_c, epsilon_M,
             Omega_chi0, w0_chi, wa_chi, plamb_alpha, dM_mode, nint):
    
    if model_name.upper() == 'LCDM':
        return mu_model(z, H0, Om, Ol, 0.0, 0.0, 0.0, 0.0, dM_mode, nint)
    
    if model_name.upper() == 'KDS':
        return mu_model(z, H0, Om, Ol, A_acc, n_acc, 0.0, 0.0, dM_mode, nint)
    
    if model_name.upper() in {'SU2', 'SU2R'}:
        return mu_model(z, H0, Om, Ol, 0.0, 0.0, 0.0, 0.0, dM_mode, nint,
                        Omega_chi0, w0_chi, wa_chi)

    if model_name.upper() == 'PLAMB':
        return mu_model_plamb_noexp(z, H0, A_acc, n_acc, gamma_c, epsilon_M,
                                    plamb_alpha, dM_mode, nint)
    
    return mu_model(z, H0, Om, Ol, A_acc, n_acc, gamma_c, epsilon_M, dM_mode, nint)


# ------------- priors ----------------------
def _bounds_default():
    return {
            'H0'        : (40.0, 100.0),
            'Om'        : ( 0.0,   1.5),
            'Ol'        : (-0.5,   2.0),
            'A_acc'     : (-2.0,   2.0),
            'n_acc'     : (-6.0,   6.0),
            'gamma_c'   : (-0.6,   0.6),
            'epsilon_M' : (-0.6,   0.6),
            'Omega_chi0': ( 0.0,   1.0),
            'Omega_dark0': ( 0.0,   1.5),
            'su2_fraction': ( 0.0,  1.0),
            'w0_chi'    : (-1.5,  -0.3),
            'wa_chi'    : (-2.0,   2.0),
            'plamb_alpha': (-0.5,   1.5),
            'Ok'        : (-2.0,   2.0),
        }


def _bounds_tight():
    return {
        'H0'            : (50.0,   80.0),
        'Om'            : ( 0.000,  0.40),
        'Ol'            : (-0.10,   1.00),
        'A_acc'         : ( 0.00,   1.50),
        'n_acc'         : ( 0.00,   2.50),
        'gamma_c'       : (-0.70,   0.10),
        'epsilon_M'     : (-0.10,   0.60),
        'Omega_chi0'    : ( 0.0,    0.8),
        'Omega_dark0'   : ( 0.0,    1.0),
        'su2_fraction'  : ( 0.0,    1.0),
        'w0_chi'        : (-1.2,   -0.7),
        'wa_chi'        : (-0.8,    0.8),
        'plamb_alpha'   : ( 0.0,    1.2),
        'Ok'            : (-1.0,    1.0),
    }


def lnprior(H0, Om, Ol, A_acc, n_acc, gamma_c, epsilon_M,
            Omega_chi0=0.0, w0_chi=-1.0, wa_chi=0.0,
            plamb_alpha=0.5, bounds=None):
    
    if isinstance(Omega_chi0, dict) and bounds is None:
        bounds      = Omega_chi0
        Omega_chi0  =  0.0
        w0_chi      = -1.0
        wa_chi      =  0.0
        plamb_alpha =  0.5
    if isinstance(plamb_alpha, dict) and bounds is None:
        bounds       = plamb_alpha
        plamb_alpha  = 0.5
    
    if bounds is None:
        bounds      = _bounds_default()
    
    H0_lo, H0_hi    = bounds['H0']
    Om_lo, Om_hi    = bounds['Om']
    Ol_lo, Ol_hi    = bounds['Ol']
    Aa_lo, Aa_hi    = bounds['A_acc']
    na_lo, na_hi    = bounds['n_acc']
    gc_lo, gc_hi    = bounds['gamma_c']
    eM_lo, eM_hi    = bounds['epsilon_M']
    Oc_lo, Oc_hi    = bounds['Omega_chi0']
    w0_lo, w0_hi    = bounds['w0_chi']
    wa_lo, wa_hi    = bounds['wa_chi']
    pa_lo, pa_hi    = bounds['plamb_alpha']
    Ok_lo, Ok_hi    = bounds['Ok']

    Ok              = 1.0 - Om - Ol - Omega_chi0

    if not (H0_lo <= H0 <= H0_hi):         return -np.inf
    if not (Om_lo <= Om <= Om_hi):         return -np.inf
    if not (Ol_lo <= Ol <= Ol_hi):         return -np.inf
    if not (Aa_lo <= A_acc <= Aa_hi):      return -np.inf
    if not (na_lo <= n_acc <= na_hi):      return -np.inf
    if not (gc_lo <= gamma_c <= gc_hi):    return -np.inf
    if not (eM_lo <= epsilon_M <= eM_hi):  return -np.inf
    if not (Oc_lo <= Omega_chi0 <= Oc_hi): return -np.inf
    if not (w0_lo <= w0_chi <= w0_hi):     return -np.inf
    if not (wa_lo <= wa_chi <= wa_hi):     return -np.inf
    if not (pa_lo <= plamb_alpha <= pa_hi): return -np.inf
    if not (Ok_lo <= Ok <= Ok_hi):         return -np.inf
    
    return 0.0


# -------- likelihood helpers ---------------
def chi2_from_residuals(resid, sigma=None, Cinv=None):
    if np.any(~np.isfinite(resid)):
        return np.inf
    
    if Cinv is not None:
        return float(resid.T @ Cinv @ resid)
    
    if sigma is not None:
        return float(np.sum((resid / sigma)**2))
    
    raise ValueError("Provide either sigma or Cinv.")


# ----------------- Extra probes I/O -----------------
def _smart_load_table(path):
    """Load CSV or whitespace table with named columns; returns (arr, names)."""
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            s       = line.strip()
            
            if not s or s.startswith('#'):
                continue
            
            delim   = ',' if (',' in s) else None
            break
        
        else:
            raise ValueError(f"No data in {path}")
    
    arr             = np.genfromtxt(path, names=True, dtype=None, encoding=None, delimiter=delim, autostrip=True)
    
    return arr, arr.dtype.names


# BAO
def load_bao(path, cov_path=None):
    arr, names      = _smart_load_table(path)
    
    # Column name normalization
    def _find(name_alt):
        for n in names:
            
            if n.lower() == name_alt:
                return n
            
        return None

    z_col           = _find('z')
    kindcol         = _find('kind') or _find('type')
    valcol          = _find('value') or _find('y') or _find('val')
    sigcol          = _find('sigma') or _find('err') or _find('error')

    if not (z_col and kindcol and valcol):
        raise ValueError(f"BAO file must have z, kind, value[, sigma]. Found: {names}")

    z               = np.asarray(arr[z_col], dtype=float)
    kinds           = [str(k).strip() for k in np.asarray(arr[kindcol])]
    vals            = np.asarray(arr[valcol], dtype=float)
    
    sig             = None
    
    if sigcol is not None:
        sig         = np.asarray(arr[sigcol], dtype=float)

    C               = None
    Cinv            = None
    
    if cov_path:
        try:
            C       = np.loadtxt(cov_path, dtype=float)
            C       = 0.5 * (C + C.T)
            
            # Small jitter
            eps     = 1e-12 + 1e-9 * float(np.mean(np.diag(C)))
            
            np.fill_diagonal(C, np.diag(C) + eps)
            Cinv    = np.linalg.pinv(C)
            
        except Exception as e:
            warnings.warn(f"[BAO] covariance load failed ({e}); using diagonal uncertainties.")
            C       = None

    return {'z': z, 'kind': kinds, 'y': vals, 'sigma': sig, 'C': C, 'Cinv': Cinv}


# CC
def load_cc(path):
    arr, names      = _smart_load_table(path)
    z_col           = None; H_col = None; s_col = None
    for n in names:
        ln          = n.lower()
        
        if ln in ('z',):                    z_col = n
        if ln in ('h', 'hz', 'h_z'):        H_col = n
        if ln in ('sigma', 'err', 'error'): s_col = n
        
    if not (z_col and H_col and s_col):
        raise ValueError(f"CC file must have z, H, sigma columns. Found: {names}")
    
    z               = np.asarray(arr[z_col], dtype=float)
    H               = np.asarray(arr[H_col], dtype=float)
    s               = np.asarray(arr[s_col], dtype=float)
    
    return {'z': z, 'H': H, 'sigma': s}


# Planck distance prior
def load_planck_prior(path):
    if path.lower().endswith('.json'):
        with open(path, 'r', encoding='utf-8') as f:
            obj     = json.load(f)
            
        means       = obj.get('means')
        cov         = obj.get('cov')
        
        if isinstance(means, dict):
            lA      = float(means.get('lA'))
            R       = float(means.get('R'))
            m       = np.array([lA, R], dtype=float)
            
        else:
            m       = np.asarray(means, dtype=float).reshape(2)
            
        C           = np.asarray(cov, dtype=float).reshape(2,2)
        
    else:
        # NPZ expected with 'means' and 'cov'
        npz         = np.load(path, allow_pickle=True)
        m           = np.asarray(npz['means'], dtype=float).reshape(2)
        C           = np.asarray(npz['cov'], dtype=float).reshape(2,2)

    C               = 0.5 * (C + C.T)
    eps             = 1e-12 + 1e-9 * float(np.mean(np.diag(C)))
    
    np.fill_diagonal(C, np.diag(C) + eps)
    
    Cinv            = np.linalg.pinv(C)
    
    return {'means': m, 'C': C, 'Cinv': Cinv}


# -------- per-probe predictions/chi2 -------
def predict_bao_vector(bao, H0, Om, Ol, A_acc, n_acc, gamma_c, nint_base, rd,
                       Omega_chi0=0.0, w0_chi=-1.0, wa_chi=0.0, model_name=''):
    
    z               = bao['z']
    
    if model_name.upper() == 'PLAMB':
        DM, DH, DV, _ = plamb_noexp_distances(z, H0, A_acc, n_acc, gamma_c, nint_base=nint_base)
    else:
        # Distances in Mpc
        DC          = comoving_distance(z, H0, Om, Ol, A_acc, n_acc, gamma_c, nint_base=nint_base,
                                        Omega_chi0=Omega_chi0, w0_chi=w0_chi, wa_chi=wa_chi)
        
        DM          = transverse_comoving_distance(DC, H0, Om, Ol, Omega_chi0)
        Ez          = E_of_z(z, Om, Ol, A_acc, n_acc, Omega_chi0, w0_chi, wa_chi)
        H           = H0 * Ez
        DH          = C_KMS / H
        DV          = (z * DH * DM**2)**(1.0/3.0)

    preds           = []
    for i, kind in enumerate(bao['kind']):
        k           = kind.strip().lower().replace(' ', '').replace('_', '')
        over_rd     = False
        
        if 'overrd' in k or '/rd' in k:
            over_rd = True
            k       = k.replace('overrd', '').replace('/rd', '')
            
        if k in ('dm', 'd_m'):
            val     = DM[i]
            
        elif k in ('dh', 'd_h'):
            val     = DH[i]
            
        elif k in ('dv', 'd_v'):
            val     = DV[i]
            
        else:
            raise ValueError(f"[BAO] Unknown kind='{kind}' (row {i}); use DM, DH, DV or *_over_rd")
        
        if over_rd:
            val     = val / rd
            
        preds.append(val)
        
    return np.asarray(preds, dtype=float)


def chi2_bao(bao, H0, Om, Ol, A_acc, n_acc, gamma_c, nint_base, rd,
             Omega_chi0=0.0, w0_chi=-1.0, wa_chi=0.0, model_name=''):
    
    yth             = predict_bao_vector(bao, H0, Om, Ol, A_acc, n_acc, gamma_c, nint_base, rd,
                                         Omega_chi0, w0_chi, wa_chi, model_name)
    resid           = yth - bao['y']
    
    if bao['Cinv'] is not None:
        return float(resid.T @ bao['Cinv'] @ resid)
    
    if bao['sigma'] is None:
        raise ValueError("[BAO] Need either --bao-cov or 'sigma' column.")
        
    return float(np.sum((resid / bao['sigma'])**2))


def chi2_cc(cc, H0, Om, Ol, A_acc, n_acc, gamma_c,
            Omega_chi0=0.0, w0_chi=-1.0, wa_chi=0.0, model_name=''):
    
    z               = cc['z']
    if model_name.upper() == 'PLAMB':
        _, DH, _, _ = plamb_noexp_distances(z, H0, A_acc, n_acc, gamma_c, nint_base=200)
        H           = C_KMS / np.maximum(DH, 1e-12)
    else:
        Ez          = E_of_z(z, Om, Ol, A_acc, n_acc, Omega_chi0, w0_chi, wa_chi)
        H           = H0 * Ez
    resid           = H - cc['H']
    
    return float(np.sum((resid / cc['sigma'])**2))


def predict_planck_vector(planck, H0, Om, Ol, A_acc, n_acc, gamma_c, zstar, rs, nint_fixed,
                          Omega_chi0=0.0, w0_chi=-1.0, wa_chi=0.0, model_name=''):
    
    if model_name.upper() == 'PLAMB':
        DM          = plamb_noexp_distance_fixed_panels(zstar, H0, A_acc, n_acc, gamma_c, n_panels=nint_fixed)
    else:
        # D_M(z*) with fixed Simpson panels for speed
        DCs         = comoving_distance_fixed_panels(zstar, H0, Om, Ol, A_acc, n_acc, gamma_c, n_panels=nint_fixed,
                                                     Omega_chi0=Omega_chi0, w0_chi=w0_chi, wa_chi=wa_chi)
        DM          = transverse_comoving_distance(np.array([DCs]), H0, Om, Ol, Omega_chi0)[0]
    lA              = np.pi * DM / rs
    
    # Standard definition with c0 (not variable-c) in denominator:
    R               = np.sqrt(Om) * H0 * DM / C_KMS
    
    return np.array([lA, R], dtype=float)


def chi2_planck(planck, H0, Om, Ol, A_acc, n_acc, gamma_c, zstar, rs, nint_fixed,
                Omega_chi0=0.0, w0_chi=-1.0, wa_chi=0.0, model_name=''):
    
    yth             = predict_planck_vector(planck, H0, Om, Ol, A_acc, n_acc, gamma_c, zstar, rs, nint_fixed,
                                  Omega_chi0, w0_chi, wa_chi, model_name)
    resid           = yth - planck['means']
    
    return float(resid.T @ planck['Cinv'] @ resid)


# -------- SN likelihood & total -----------
def lnlike_sn(H0, Om, Ol, A_acc, n_acc, gamma_c, epsilon_M,
              Omega_chi0, w0_chi, wa_chi, plamb_alpha,
              data, Cinv, model_name, dM_mode, nint):
    
    mu_th           = model_mu(data['z'], model_name, H0, Om, Ol, A_acc, n_acc, gamma_c, epsilon_M,
                               Omega_chi0, w0_chi, wa_chi, plamb_alpha, dM_mode, nint)
    if not np.all(np.isfinite(mu_th)):
        return -np.inf, np.inf
    
    resid           = mu_th - data['mu']
    chi2            = chi2_from_residuals(resid, sigma=data.get('sigma', None), Cinv=Cinv)
    
    if not np.isfinite(chi2):
        return -np.inf, np.inf
    
    return -0.5 * chi2, chi2


def total_lnlike(H0, Om, Ol, A_acc, n_acc, gamma_c, epsilon_M,
                 Omega_chi0, w0_chi, wa_chi, plamb_alpha,
                 data, Cinv, model_name, dM_mode, nint,
                 extras):
    
    # SN
    if extras.get('use_sn', True):
        ln_sn, chi2_sn = lnlike_sn(H0, Om, Ol, A_acc, n_acc, gamma_c, epsilon_M,
                                   Omega_chi0, w0_chi, wa_chi, plamb_alpha,
                                   data, Cinv, model_name, dM_mode, nint)
        
        if not np.isfinite(ln_sn):
            return -np.inf, {'SN': np.inf, 'BAO': 0.0, 'CC': 0.0, 'Planck': 0.0}

        chi2_map        = {'SN': chi2_sn, 'BAO': 0.0, 'CC': 0.0, 'Planck': 0.0}
        ln_tot          = ln_sn
    else:
        chi2_map        = {'SN': 0.0, 'BAO': 0.0, 'CC': 0.0, 'Planck': 0.0}
        ln_tot          = 0.0

    # BAO
    if extras.get('bao') is not None:
        try:
            c2      = chi2_bao(extras['bao'], H0, Om, Ol, A_acc, n_acc, gamma_c, nint, extras['rd'],
                               Omega_chi0, w0_chi, wa_chi, model_name)
        
        except Exception as e:
            warnings.warn(f"[BAO] skipping due to error: {e}")
            c2      = 0.0
        
        chi2_map['BAO']     = c2
        ln_tot             += -0.5 * c2

    # CC
    if extras.get('cc') is not None:
        try:
            c2              = chi2_cc(extras['cc'], H0, Om, Ol, A_acc, n_acc, gamma_c, Omega_chi0, w0_chi, wa_chi, model_name)
        
        except Exception as e:
            warnings.warn(f"[CC] skipping due to error: {e}")
            c2              = 0.0
        chi2_map['CC']      = c2
        ln_tot             += -0.5 * c2

    # Planck
    if extras.get('planck') is not None:
        try:
            c2              = chi2_planck(extras['planck'], H0, Om, Ol, A_acc, n_acc, gamma_c, extras['zstar'], extras['rs'], extras['planck_nint'],
                                          Omega_chi0, w0_chi, wa_chi, model_name)
        
        except Exception as e:
            warnings.warn(f"[Planck] skipping due to error: {e}")
            c2              = 0.0
            
        chi2_map['Planck']  = c2
        ln_tot             += -0.5 * c2

    return ln_tot, chi2_map


def effective_observation_count(data, extras):
    if extras.get('use_sn', True):
        return int(data['N'])

    n_obs = 0
    if extras.get('bao') is not None:
        n_obs += int(len(extras['bao']['z']))
    if extras.get('cc') is not None:
        n_obs += int(len(extras['cc']['z']))
    if extras.get('planck') is not None:
        n_obs += int(len(extras['planck']['means']))
    return max(n_obs, 1)


# --------- active-parameter machinery ------
def build_parameter_sets(args):
    theta0_dict         = {
                           'H0'        : args.H0,
                           'Om'        : args.Om,
                           'Ol'        : args.Ol,
                           'A_acc'     : args.A_acc,
                           'n_acc'     : args.n_acc,
                           'gamma_c'   : args.gamma_c,
                           'epsilon_M' : args.epsilon_M,
                           'Omega_chi0': args.Omega_chi0,
                           'Omega_dark0': args.Omega_dark0,
                           'su2_fraction': args.su2_fraction,
                           'w0_chi'    : args.w0_chi,
                           'wa_chi'    : args.wa_chi,
                           'plamb_alpha': args.plamb_alpha,
                          }
    frozen_flags        = {
                          'H0'        : args.fix_H0,
                          'Om'        : args.fix_Om,
                          'Ol'        : args.fix_Ol,
                          'A_acc'     : args.fix_A_acc,
                          'n_acc'     : args.fix_n_acc,
                          'gamma_c'   : args.fix_gamma_c,
                          'epsilon_M' : args.fix_epsilon_M,
                          'Omega_chi0': args.fix_Omega_chi0,
                          'Omega_dark0': args.fix_Omega_dark0,
                          'su2_fraction': args.fix_su2_fraction,
                          'w0_chi'    : args.fix_w0_chi,
                          'wa_chi'    : args.fix_wa_chi,
                          'plamb_alpha': not args.free_plamb_alpha,
                         }
    
    active_names        = [k for k in theta0_dict if not frozen_flags[k]]
    frozen_items        = {k: v for k, v in theta0_dict.items() if frozen_flags[k]}
    theta0_active       = np.array([theta0_dict[k] for k in active_names], dtype=float)
    
    return theta0_dict, frozen_flags, active_names, frozen_items, theta0_active


def unpack_active_factory(active_names, frozen_items, plamb_conservation_mode='free', use_su2r_split=False):
    def _unpack(theta_active):
        full            = dict(frozen_items)
        full.update({k: v for k, v in zip(active_names, theta_active)})

        if use_su2r_split and 'su2_fraction' in full:
            dark        = full.get('Omega_dark0', full.get('Ol', 0.0) + full.get('Omega_chi0', 0.0))
            frac        = full['su2_fraction']
            full['Omega_dark0'] = dark
            full['Omega_chi0']  = dark * frac
            full['Ol']          = dark * (1.0 - frac)

        apply_plamb_conservation(full, plamb_conservation_mode)
        
        return (full['H0'],         full['Om'],     full['Ol'],
                full['A_acc'],      full['n_acc'],
                full['gamma_c'],    full['epsilon_M'],
                full['Omega_chi0'], full['w0_chi'], full['wa_chi'],
                full['plamb_alpha'])
    
    return _unpack


def unpack_active_factory_flat(active_names, frozen_items, plamb_conservation_mode='free', use_su2r_split=False):
    def _unpack(theta_active):
        full            = dict(frozen_items)
            
        full.update({k: v for k, v in zip(active_names, theta_active)})
        Om              = full['Om']
        if use_su2r_split and 'su2_fraction' in full:
            dark        = 1.0 - Om
            frac        = full['su2_fraction']
            full['Omega_dark0'] = dark
            full['Omega_chi0']  = dark * frac
            full['Ol']          = dark * (1.0 - frac)
        else:
            full['Ol']  = 1.0 - Om - full.get('Omega_chi0', 0.0)

        apply_plamb_conservation(full, plamb_conservation_mode)
    
        return (full['H0'],         full['Om'],         full['Ol'],
                full['A_acc'],      full['n_acc'],
                full['gamma_c'],    full['epsilon_M'],
                full['Omega_chi0'], full['w0_chi'],     full['wa_chi'],
                full['plamb_alpha'])
    
    return _unpack


def lnprior_active_coords(theta_active, active_names, bounds):
    active             = {k: v for k, v in zip(active_names, theta_active)}
    if 'su2_fraction' in active:
        lo, hi         = bounds['su2_fraction']
        if not (lo <= active['su2_fraction'] <= hi):
            return -np.inf
    if 'Omega_dark0' in active:
        lo, hi         = bounds['Omega_dark0']
        if not (lo <= active['Omega_dark0'] <= hi):
            return -np.inf
    return 0.0


def lnprob_active(theta_active, data, Cinv, model_name, dM_mode, nint, bounds, unpack, active_names, extras):
    lp_active          = lnprior_active_coords(theta_active, active_names, bounds)
    if not np.isfinite(lp_active):
        return -np.inf

    H0, Om, Ol, A_acc, n_acc, gamma_c, epsilon_M, Omega_chi0, w0_chi, wa_chi, plamb_alpha = unpack(theta_active)
    
    lp                  = lnprior(H0, Om, Ol, A_acc, n_acc, gamma_c, epsilon_M,
                                  Omega_chi0, w0_chi, wa_chi, plamb_alpha, bounds)
    
    if not np.isfinite(lp):
        return -np.inf
    
    ll, _               = total_lnlike(H0, Om, Ol, A_acc, n_acc, gamma_c, epsilon_M,
                                       Omega_chi0, w0_chi, wa_chi, plamb_alpha,
                                       data, Cinv, model_name, dM_mode, nint, extras)
    return lp_active + lp + ll


# --------------- diagnostics ----------------
def save_corner_and_traces(outdir, flat_chain, labels, suptitle=""):
    if flat_chain.size == 0:
        return
    
    # Corner
    if corner is not None:
        fig_c           = corner.corner(flat_chain, labels=labels, show_titles=True, title_fmt=".3f")
    
        if suptitle:
            fig_c.suptitle(suptitle, fontsize=9, y=0.995)
        
        fig_c.savefig(os.path.join(outdir, "corner.png"), dpi=160, bbox_inches="tight")
        if plt is not None:
            plt.close(fig_c)   # guard
    
    # Traces
    if plt is None:
        return
    
    npar            = flat_chain.shape[1]
    fig_t, axs      = plt.subplots(npar, 1, figsize=(7.5, 1.8*npar), constrained_layout=True)
    
    if npar == 1:
        axs         = [axs]
        
    for i, ax in enumerate(axs):
        ax.plot(flat_chain[:, i], lw=0.6)
        ax.set_ylabel(labels[i])
    
    axs[-1].set_xlabel("sample index (post-burn)")
    
    if suptitle:
        fig_t.suptitle(suptitle, fontsize=9, y=1.02)
        
    fig_t.savefig(os.path.join(outdir, "traces.png"), dpi=160, bbox_inches="tight")
    plt.close(fig_t)


# Binned RMS of residuals (visual check)
def binned_rms(x, y, nbins=18, min_count=10):
    edges           = np.linspace(x.min(), x.max(), nbins + 1)
    xc, rms         = [], []
    for i in range(nbins):
        m           = (x >= edges[i]) & (x < edges[i+1])
        
        if np.count_nonzero(m) >= min_count:
            xc.append(0.5 * (edges[i] + edges[i+1]))
            rms.append(np.std(y[m], ddof=1))
    
    return np.array(xc), np.array(rms)


# ---------------- main ---------------------
def main():
    parser          = argparse.ArgumentParser(description="PLamb Pantheon+ fitting â€” V10V6c_plus")
    
    parser.add_argument('--model',          default='FR', choices=['LCDM','KdS','FR','SU2','SU2R','PLAMB'])
    parser.add_argument('--data',           required=True)
    parser.add_argument('--cov',            default='')

    parser.add_argument('--dM-mode',        default='linear_z', choices=['linear_z','log1pz'])

    parser.add_argument('--H0',             type=float,     default=70.0)
    parser.add_argument('--Om',             type=float,     default=0.3)
    parser.add_argument('--Ol',             type=float,     default=0.7)
    parser.add_argument('--A_acc',          type=float,     default=0.0)
    parser.add_argument('--n_acc',          type=float,     default=0.0)
    parser.add_argument('--gamma_c',        type=float,     default=0.0)
    parser.add_argument('--epsilon_M',      type=float,     default=0.0)
    parser.add_argument('--Omega_chi0',     type=float,     default=0.0,             help='Present-day density fraction for the phenomenological chiral SU(2)-like component')
    parser.add_argument('--Omega_dark0',    type=float,     default=None,            help='SU2R total dark-sector density before splitting into Lambda and SU2-like pieces')
    parser.add_argument('--su2-fraction',   dest='su2_fraction', type=float, default=None,
                        help='SU2R fraction of Omega_dark0 assigned to the SU2-like component')
    parser.add_argument('--w0-chi',         dest='w0_chi', type=float, default=-1.0, help='Present-day equation of state for Omega_chi(z)')
    parser.add_argument('--wa-chi',         dest='wa_chi', type=float, default=0.0,  help='CPL evolution parameter for Omega_chi(z)')
    parser.add_argument('--plamb-alpha',    dest='plamb_alpha', type=float, default=0.5,
                        help='PLAMB luminosity dimming exponent in DL=(1+z)^alpha D_path')
    parser.add_argument('--free-plamb-alpha', action='store_true',
                        help='Sample PLAMB luminosity dimming exponent alpha; default fixes alpha at --plamb-alpha')
    parser.add_argument('--plamb-conservation-mode',
                        choices=['free', 'Gc_const', 'Gc2_const'],
                        default='free',
                        help='PLAMB Noether diagnostic scaffold: free, G/c constant, or G/c^2 constant. Non-free modes derive epsilon_M from gamma_c.')

    parser.add_argument('--fix_H0',         action='store_true')
    parser.add_argument('--fix_Om',         action='store_true')
    parser.add_argument('--fix_Ol',         action='store_true')
    parser.add_argument('--fix_A_acc',      action='store_true')
    parser.add_argument('--fix_n_acc',      action='store_true')
    parser.add_argument('--fix_gamma_c',    action='store_true')
    parser.add_argument('--fix_epsilon_M',  action='store_true')
    parser.add_argument('--fix_Omega_chi0', action='store_true')
    parser.add_argument('--fix_Omega_dark0', action='store_true')
    parser.add_argument('--fix_su2_fraction', action='store_true')
    parser.add_argument('--fix_w0_chi',     action='store_true')
    parser.add_argument('--fix_wa_chi',     action='store_true')

    parser.add_argument('--sampler',                        default='emcee', choices=['emcee','grid','none'])
    parser.add_argument('--nwalkers',       type=int,       default=32)
    parser.add_argument('--nsteps',         type=int,       default=2000)
    parser.add_argument('--nburn',          type=int,       default=500)
    parser.add_argument('--seed',           type=int,       default=None)
    parser.add_argument('--p0-frac',        type=float,     default=0.05,   help='Default fractional scatter for initial emcee walkers')
    parser.add_argument('--p0-min-scale',   type=float,     default=1e-3,   help='Minimum absolute scatter for initial emcee walkers')
    parser.add_argument('--p0-scales',      default='',                     help='Comma-separated absolute initial walker scatters, in active-parameter order')

    parser.add_argument('--grid-points',    type=int,       default=1000)
    
    parser.add_argument('--grid-min',       type=float,     default=None,   help='Explicit lower bound for one-parameter grid mode')
    
    parser.add_argument('--grid-max',       type=float,     default=None,   help='Explicit upper bound for one-parameter grid mode')
    
    parser.add_argument('--grid-width-frac', type=float,    default=0.2,    help='Fractional half-width around the start value for grid mode when --grid-min/--grid-max are not set')
    
    parser.add_argument('--nint',           type=int,       default=200)
    parser.add_argument('--thin',           type=int,       default=1)
    parser.add_argument('--zmax',           type=float,     default=10.0)
    
    parser.add_argument('--tight-priors',   action='store_true')
    parser.add_argument('--corner',         action='store_true')

    parser.add_argument('--out',            default='plamb_runs/10V6c_plus')
    parser.add_argument('--show',           action='store_true')
    parser.add_argument('--verbose',        action='store_true')
    
    parser.add_argument('--flat',           action='store_true',            help='Enforce Î©k=0 â‡’ Î©Î› = 1 âˆ’ Î©m (remove Î©Î› from the sampled space)')
    
    parser.add_argument('--zmin',           type=float,     default=0.0,    help='Minimum redshift to include (PV robustness check)')
    
    parser.add_argument('--save-chain',     action='store_true',            help='Save flat_chain.npy and flat_ln.npy to --out')
    parser.add_argument('--no-sn',          action='store_true',            help='Disable SN likelihood; use only extra probes such as BAO/CC/Planck')

    # --- New extra-probe flags ---
    parser.add_argument('--bao',                            default='',     help='BAO table path')
    parser.add_argument('--bao-cov',                        default='',     help='BAO covariance path (optional)')
    parser.add_argument('--rd',             type=float,     default=147.09, help='BAO r_d [Mpc]')

    parser.add_argument('--cc',                             default='',     help='Cosmic-chronometer table path')

    parser.add_argument('--planck',                         default='',     help='Planck distance prior (JSON/NPZ with lA,R and 2x2 cov)')
    parser.add_argument('--planck-zstar',   type=float,     default=1089.0, help='Decoupling redshift z_*')
    parser.add_argument('--rs',             type=float,     default=None,   help='r_s used in l_A (defaults to --rd)')
    parser.add_argument('--planck-nint',    type=int,       default=4000,   help='Fixed Simpson panels for z~z_* integral')

    args = parser.parse_args()

    if args.Omega_dark0 is None:
        args.Omega_dark0 = (1.0 - args.Om) if args.flat else (args.Ol + args.Omega_chi0)
    if args.su2_fraction is None:
        dark0 = args.Omega_dark0
        args.su2_fraction = (args.Omega_chi0 / dark0) if abs(dark0) > 1e-12 else 0.5

    # ---- model-implied freezes ----
    if args.model.upper() == 'LCDM':
        args.A_acc        = 0.0 ; args.fix_A_acc      = True
        args.n_acc        = 0.0 ; args.fix_n_acc      = True
        args.gamma_c      = 0.0 ; args.fix_gamma_c    = True
        args.epsilon_M    = 0.0 ; args.fix_epsilon_M  = True
        
    elif args.model.upper() == 'KDS':
        args.gamma_c      = 0.0 ; args.fix_gamma_c    = True
        args.epsilon_M    = 0.0 ; args.fix_epsilon_M  = True
        
    elif args.model.upper() in {'SU2', 'SU2R'}:
        args.A_acc        = 0.0 ; args.fix_A_acc      = True
        args.n_acc        = 0.0 ; args.fix_n_acc      = True
        args.gamma_c      = 0.0 ; args.fix_gamma_c    = True
        args.epsilon_M    = 0.0 ; args.fix_epsilon_M  = True
        if args.model.upper() == 'SU2R':
            args.fix_Ol          = True
            args.fix_Omega_chi0  = True
            if args.flat:
                args.fix_Omega_dark0 = True

    elif args.model.upper() == 'PLAMB':
        # Peter Lamb no-expansion/inertia-clumping pilot:
        #   A_acc     -> clump_amp
        #   n_acc     -> clump_index
        #   gamma_c   -> c_drift
        #   epsilon_M -> inertia_drift
        # It is intentionally separated from SU2/LCDM expansion closure.
        args.Ol           = 0.0 ; args.fix_Ol         = True
        if args.n_acc == 0.0 and not args.fix_n_acc:
            args.n_acc    = 1.0
        if args.flat:
            warnings.warn("[PLAMB] ignoring --flat; this no-expansion pilot does not use Lambda closure.")
            args.flat     = False
        if args.plamb_conservation_mode != 'free':
            if args.fix_epsilon_M is False and args.epsilon_M != 0.0:
                warnings.warn("[PLAMB] --plamb-conservation-mode derives epsilon_M from gamma_c; ignoring supplied --epsilon_M.")
            if args.free_plamb_alpha:
                warnings.warn("[PLAMB] --plamb-conservation-mode keeps plamb_alpha fixed; ignoring --free-plamb-alpha.")
            args.fix_epsilon_M = True
            args.free_plamb_alpha = False

    if args.model.upper() != 'SU2R':
        args.fix_Omega_dark0  = True
        args.fix_su2_fraction = True

    if args.model.upper() not in {'SU2', 'SU2R'}:
        args.Omega_chi0   = 0.0 ; args.fix_Omega_chi0 = True
        args.w0_chi       = -1.0; args.fix_w0_chi     = True
        args.wa_chi       = 0.0 ; args.fix_wa_chi     = True

    if args.model.upper() != 'PLAMB':
        args.free_plamb_alpha = False
        args.plamb_conservation_mode = 'free'

    if args.seed is not None:
        np.random.seed(int(args.seed))

    print("[status] parsing done.")
    if args.verbose:
        print(f"[int] Simpson base panels set to {args.nint}")

    # ---- SNe data ----
    data_all                    = load_pantheon(args.data)
    mask                        = (data_all['z'] >= args.zmin) & (data_all['z'] <= args.zmax)
    if args.thin > 1:
        mask                   &= (np.arange(data_all['N']) % args.thin == 0)

    z                           = data_all['z'][mask]
    mu                          = data_all['mu'][mask]
    sigma                       = data_all['sigma'][mask]
    data                        = {'z': z, 'mu': mu, 'sigma': sigma, 'N': z.size}

    print("[status] data loaded.")
    if args.verbose:
        lo, hi                  = float(np.min(z)), float(np.max(z))
        print(f"[load] named columns OK: zCMB, MU_SH0ES, MU_SH0ES_ERR_DIAG  N={data['N']}")
        print(f"Loaded {data['N']} SNe from {os.path.basename(args.data)}  zâˆˆ[{lo:.4f},{hi:.4f}]")

    # ---- SN covariance ----
    C_full                      = load_covariance(args.cov, data_all['N'])
    C, Cinv                     = None, None
    
    if C_full is not None:
        try:
            mask_idx            = np.where(mask)[0] if mask.dtype == bool else np.asarray(mask, int)
            if C_full.shape == (data_all['N'], data_all['N']):
                C               = C_full[np.ix_(mask_idx, mask_idx)]
    
            elif C_full.shape == (data['N'], data['N']):
                C               = C_full
            
            else:
                warnings.warn(f"[cov] shape {C_full.shape} does not match {data_all['N']} or {data['N']}; using diagonal.")
                C               = None

            if C is not None:
                C               = 0.5 * (C + C.T)
                diag_mean       = float(np.mean(np.diag(C)))
                eps             = 1e-12 + 1e-9 * max(diag_mean, 1.0)
                np.fill_diagonal(C, np.diag(C) + eps)
                
                try:
                    Cinv        = np.linalg.inv(C)
                    if args.verbose:
                        print(f"[cov] using subsetted covariance: {C.shape}")
                
                except Exception as e:
                    warnings.warn(f"[cov] inv failed ({e}); falling back to pinv.")
                    Cinv        = np.linalg.pinv(C)
        
        except Exception as e:
            warnings.warn(f"[cov] subsetting failed ({e}); using diagonal.")
            C                   = None
    else:
        if args.verbose:
            print("[cov] no usable covariance; using diagonal errors.")

    ensure_dir(args.out)

    # ---- Build parameter lists ----
    theta0_dict, frozen_flags, active_names, frozen_items, theta0_active = build_parameter_sets(args)

    if args.flat and 'Ol' in active_names:
        active_names            = [k for k in active_names if k != 'Ol']
        theta0_active           = np.array([theta0_dict[k] for k in active_names], dtype=float)

    if args.model.upper() == 'SU2R':
        frozen_items.pop('Ol', None)
        frozen_items.pop('Omega_chi0', None)
        if args.flat:
            frozen_items.pop('Omega_dark0', None)
        
    use_su2r_split              = args.model.upper() == 'SU2R'
    unpack_active               = (unpack_active_factory_flat(active_names, frozen_items, args.plamb_conservation_mode, use_su2r_split)
                               
                     if   args.flat 
                     else unpack_active_factory(active_names, frozen_items, args.plamb_conservation_mode, use_su2r_split))

    bounds                      = _bounds_tight() if args.tight_priors else _bounds_default()
    if args.model.upper() == 'PLAMB':
        bounds.update({
            'H0'        : (50.0, 85.0),
            'Om'        : (0.0, 1.5),      # interpreted as Omega_inertia for the Planck R stress test
            'Ol'        : (-1e-8, 1e-8),
            'A_acc'     : (-1.5, 1.5),     # clump_amp
            'n_acc'     : (0.05, 4.0),     # clump_index
            'gamma_c'   : (-0.02, 0.02),   # c_drift
            'epsilon_M' : (-1.0, 1.0),     # inertia_drift / mass-luminosity drift
            'Omega_chi0': (0.0, 0.0),
            'w0_chi'    : (-1.0, -1.0),
            'wa_chi'    : (0.0, 0.0),
            'plamb_alpha': (-0.5, 1.5),
            'Ok'        : (-2.0, 2.0),
        })
    
    if args.flat:
        bounds['Ok']            = (-1e-8, 1e-8)

    print(f"[status] active params: {active_names}  frozen: {frozen_items}")

    # ---- Extras: BAO / CC / Planck ----
    extras                      = {'use_sn': not args.no_sn}
    if args.no_sn and args.verbose:
        print("[sn] likelihood disabled by --no-sn")
    
    if args.bao:
        extras['bao']           = load_bao(args.bao, args.bao_cov if args.bao_cov else None)
        extras['rd']            = float(args.rd)
        
        if args.verbose:
            Nbao                = len(extras['bao']['z'])
            print(f"[bao] loaded {Nbao} rows  (r_d={extras['rd']:.3f} Mpc)")

    if args.cc:
        extras['cc']            = load_cc(args.cc)
        
        if args.verbose:
            print(f"[cc]   loaded {extras['cc']['z'].size} points")

    if args.planck:
        extras['planck']        = load_planck_prior(args.planck)
        extras['zstar']         = float(args.planck_zstar)
        extras['rs']            = float(args.rs) if (args.rs is not None) else float(args.rd)
        extras['planck_nint']   = int(args.planck_nint)
        
        if args.verbose:
            print(f"[planck] distance prior active  z*={extras['zstar']}  r_s={extras['rs']}  nint_fixed={extras['planck_nint']}")

    # ---- Sampler / evaluator ----
    bestfit                     = None
    lnpost_best                 = -np.inf
    t0                          = time.time()
    labels_map                  = {'H0':'H0','Om':'Î©m','Ol':'Î©Î›','A_acc':'A_acc','n_acc':'n_acc',
                                   'gamma_c':'Î³_c','epsilon_M':'Îµ_M',
                                   'Omega_chi0':'Omega_chi0','w0_chi':'w0_chi','wa_chi':'wa_chi',
                                   'Omega_dark0':'Omega_dark0','su2_fraction':'su2_fraction',
                                   'plamb_alpha':'plamb_alpha'}
    if args.model.upper() == 'PLAMB':
        labels_map.update({
            'Om'        : 'Omega_inertia',
            'A_acc'     : 'clump_amp',
            'n_acc'     : 'clump_index',
            'gamma_c'   : 'c_drift',
            'epsilon_M' : 'inertia_drift',
            'plamb_alpha': 'dimming_alpha',
        })
    labels_active               = [labels_map[k] for k in active_names]

    flat_chain                  = np.empty((0, len(active_names)))

    if args.sampler == 'emcee':
        if emcee is None:
            raise RuntimeError("emcee not available. Install 'emcee' or use --sampler grid/none.")

        ndim                    = theta0_active.size
        if ndim == 0:
            print("[warn] no active parameters; evaluating single point only.")
            bestfit         = unpack_active(theta0_active)
            lp              = lnprior(*bestfit, bounds)
            lnpost_best, chi2_map = total_lnlike(*bestfit, data, Cinv, args.model, args.dM_mode, args.nint, extras)
            lnpost_best    += lp
        
        else:
            if args.p0_scales.strip():
                base_scale  = np.array([float(x.strip()) for x in args.p0_scales.split(',') if x.strip()], dtype=float)
                
                if base_scale.size != ndim:
                    raise ValueError(f"--p0-scales requires {ndim} values for active parameters {active_names}.")
            else:
                if args.p0_frac <= 0.0 or args.p0_min_scale <= 0.0:
                    raise ValueError("--p0-frac and --p0-min-scale must be positive.")
                
                base_scale  = np.maximum(args.p0_min_scale, np.abs(theta0_active) * args.p0_frac)
            
            if not np.all(np.isfinite(base_scale)) or np.any(base_scale <= 0.0):
                raise ValueError("--p0-scales values must be finite and positive.")
            
            p0              = theta0_active[None, :] + base_scale[None, :] * np.random.randn(args.nwalkers, ndim)
            
            if args.verbose:
                print(f"[emcee] p0 scales: {dict(zip(active_names, base_scale))}")

            sampler         = emcee.EnsembleSampler(args.nwalkers, ndim,
                                                    lnprob_active,
                                                    args=(data, Cinv, args.model, args.dM_mode, args.nint, bounds, unpack_active, active_names, extras))
            if args.verbose:
                print(f"[emcee] starting: walkers={args.nwalkers}  steps={args.nsteps}  burn-in={args.nburn}  ndim={ndim}")

            try:
                sampler.run_mcmc(p0, args.nsteps, progress=args.verbose)
                
                flat_chain      = sampler.get_chain(discard=args.nburn, flat=True)
                flat_ln         = sampler.get_log_prob(discard=args.nburn, flat=True)
                acc_frac_mean   = float(np.mean(sampler.acceptance_fraction))
                
                if args.verbose:
                    print(f"[emcee] mean acceptance fraction = {acc_frac_mean:.3f}")
                
                if flat_chain.size == 0:
                    print("[warn] no valid posterior samples after burn-in; try relaxing priors or reducing nsteps/nburn.")
                    return
                
                imax            = int(np.argmax(flat_ln))
                theta_best      = flat_chain[imax]
                bestfit         = unpack_active(theta_best)
                lnpost_best     = float(flat_ln[imax])
                
                if args.save_chain:
                    np.save(os.path.join(args.out, "flat_chain.npy"), flat_chain)
                    np.save(os.path.join(args.out, "flat_ln.npy"), flat_ln)
            
            except ValueError as e:
                print(f"[error] emcee halted: {e}")
                return

    elif args.sampler == 'grid':
        if theta0_active.size != 1:
            raise ValueError("Grid mode supports exactly one active parameter. Freeze others.")
        
        pname                   = active_names[0]
        p0_val                  = theta0_active[0]
        
        if (args.grid_min is None) != (args.grid_max is None):
            raise ValueError("Set both --grid-min and --grid-max, or neither.")
        
        if args.grid_min is not None:
            if args.grid_min >= args.grid_max:
                raise ValueError("--grid-min must be smaller than --grid-max.")
            grid                = np.linspace(args.grid_min, args.grid_max, args.grid_points)
        
        else:
            width               = args.grid_width_frac * abs(p0_val) if abs(p0_val) > 0 else args.grid_width_frac
            grid                = np.linspace(p0_val - width, p0_val + width, args.grid_points)
        
        if args.verbose:
            print(f"[grid] scanning {pname}: {grid[0]:.6g} .. {grid[-1]:.6g}  N={args.grid_points}")
        
        lnposts                 = []
        
        for val in grid:
            lp_active          = lnprior_active_coords(np.array([val]), [pname], bounds)
            H0, Om, Ol, A_acc, n_acc, gamma_c, epsilon_M, Omega_chi0, w0_chi, wa_chi, plamb_alpha = unpack_active(np.array([val]))
            
            lp                  = lnprior(H0, Om, Ol, A_acc, n_acc, gamma_c, epsilon_M,
                                          Omega_chi0, w0_chi, wa_chi, plamb_alpha, bounds)
            ll, _               = total_lnlike(H0, Om, Ol, A_acc, n_acc, gamma_c, epsilon_M,
                                               Omega_chi0, w0_chi, wa_chi, plamb_alpha,
                                               data, Cinv, args.model, args.dM_mode, args.nint, extras)
            lnposts.append(lp_active + lp + ll)
        
        lnposts                 = np.asarray(lnposts)
        imax                    = int(np.argmax(lnposts))
        theta_best              = np.array([grid[imax]])
        bestfit                 = unpack_active(theta_best)
        lnpost_best             = float(lnposts[imax])

    else:
        lp_active               = lnprior_active_coords(theta0_active, active_names, bounds)
        bestfit                 = unpack_active(theta0_active)
        lp                      = lnprior(*bestfit, bounds)
        ll, _                   = total_lnlike(*bestfit, data, Cinv, args.model, args.dM_mode, args.nint, extras)
        
        if not np.isfinite(lp_active + lp):
            raise ValueError("Initial parameter set violates priors.")
        
        lnpost_best             = lp_active + lp + ll

    t1                          = time.time()
    
    if bestfit is None:
        print("[warn] no bestfit produced; exiting.")
        return

    # ---- Best-fit diagnostics & outputs ----
    H0, Om, Ol, A_acc, n_acc, gamma_c, epsilon_M, Omega_chi0, w0_chi, wa_chi, plamb_alpha = bestfit

    # Recompute SN residuals and total chi2 breakdown
    if extras.get('use_sn', True):
        mu_th                   = model_mu(z, args.model, H0, Om, Ol, A_acc, n_acc, gamma_c, epsilon_M,
                                           Omega_chi0, w0_chi, wa_chi, plamb_alpha, args.dM_mode, args.nint)
        resid                   = mu_th - mu
        chi2_sn                 = chi2_from_residuals(resid, sigma=None if C is not None else data['sigma'],
                                                      Cinv=Cinv if C is not None else None)
    else:
        mu_th                   = np.array([], dtype=float)
        resid                   = np.array([], dtype=float)
        chi2_sn                 = 0.0

    n_info                      = effective_observation_count(data, extras)
    dof                         = n_info - len(active_names)
    k_active                    = len(active_names)

    # Extras chi2
    chi2_map                    = {'SN': chi2_sn, 'BAO': 0.0, 'CC': 0.0, 'Planck': 0.0}
    
    if extras.get('bao') is not None:
        chi2_map['BAO']         = chi2_bao(extras['bao'], H0, Om, Ol, A_acc, n_acc, gamma_c, args.nint, extras['rd'],
                                           Omega_chi0, w0_chi, wa_chi, args.model)
    if extras.get('cc') is not None:
        chi2_map['CC']          = chi2_cc(extras['cc'], H0, Om, Ol, A_acc, n_acc, gamma_c, Omega_chi0, w0_chi, wa_chi, args.model)
    
    if extras.get('planck') is not None:
        chi2_map['Planck']      = chi2_planck(extras['planck'], H0, Om, Ol, A_acc, n_acc, gamma_c,
                                              extras['zstar'], extras['rs'], extras['planck_nint'],
                                              Omega_chi0, w0_chi, wa_chi, args.model)

    chi2_total                  = sum(chi2_map.values())
    AIC                         = chi2_total + 2.0 * k_active
    BIC                         = chi2_total + k_active * np.log(n_info)

    w_eff_best                  = w_eff_from_nacc(n_acc)

    ensure_dir(args.out)
    tag                         = f"{args.model}_{args.sampler}_bestfit.txt"
    outpath                     = os.path.join(args.out, tag)

    best_dict                   = {
                                   'H0'             : float(H0),
                                   'Om'             : float(Om),
                                   'Ol'             : float(Ol),
                                   'Ok'             : float(1.0 - Om - Ol - Omega_chi0),
                                   'A_acc'          : float(A_acc),
                                   'n_acc'          : float(n_acc),
                                   'gamma_c'        : float(gamma_c),
                                   'epsilon_M'      : float(epsilon_M),
                                   'Omega_chi0'     : float(Omega_chi0),
                                   'Omega_dark0'    : float(Ol + Omega_chi0),
                                   'su2_fraction'   : float(Omega_chi0 / (Ol + Omega_chi0)) if abs(Ol + Omega_chi0) > 1e-12 else float('nan'),
                                   'w0_chi'         : float(w0_chi),
                                   'wa_chi'         : float(wa_chi),
                                   'plamb_alpha'    : float(plamb_alpha),
                                   'lnpost'         : float(lnpost_best),
                                   'chi2_SN'        : float(chi2_map['SN']),
                                   'chi2_BAO'       : float(chi2_map['BAO']),
                                   'chi2_CC'        : float(chi2_map['CC']),
                                   'chi2_Planck'    : float(chi2_map['Planck']),
                                   'chi2'           : float(chi2_total),
                                   'dof'            : int(dof),
                                   'chi2/dof'       : float(chi2_total / max(dof, 1)),
                                   'AIC'            : float(AIC),
                                   'BIC'            : float(BIC),
                                   'N'              : int(n_info),
                                   'use_SN'         : bool(extras.get('use_sn', True)),
                                   'active'         : active_names,
                                   'frozen'         : {k: float(v) for k, v in frozen_items.items()},
                                   'w_eff'          : float(w_eff_best),
                                  }
    if args.model.upper() == 'SU2R':
        best_dict.update({
            'model_note'      : 'SU2R reparameterized SU2 branch; samples total dark-sector split instead of Omega_chi0 directly',
            'reparameterized_dark_sector': True,
            'dark_sector_closure': 'flat: Omega_dark0 = 1 - Om; Omega_chi0 = su2_fraction * Omega_dark0; Ol = (1-su2_fraction) * Omega_dark0',
        })
    if args.model.upper() == 'PLAMB':
        conservation = plamb_conservation_summary(args.plamb_conservation_mode, gamma_c, epsilon_M)
        best_dict.update({
            'model_note'      : 'PLAMB no-expansion/inertia-clumping pilot; phenomenological stress test, not full action-level theory',
            'Omega_inertia'   : float(Om),
            'clump_amp'       : float(A_acc),
            'clump_index'     : float(n_acc),
            'c_drift'         : float(gamma_c),
            'inertia_drift'   : float(epsilon_M),
            'dimming_alpha'   : float(plamb_alpha),
            'distance_law'    : 'D_path = c/H0 integral exp(clump_amp log(1+z)^clump_index) c(z)/c0 dz/(1+z); DL = (1+z)^alpha D_path',
            'no_expansion_pilot': True,
        })
        best_dict.update(conservation)
        if args.plamb_conservation_mode != 'free':
            best_dict['frozen']['epsilon_M'] = float(epsilon_M)
            best_dict['derived'] = {
                'epsilon_M': conservation['plamb_epsilon_M_relation'],
                'plamb_alpha': 'fixed input exponent; default 0.5 approximates sqrt(1+z) / 1+z/2 time projection at low z',
            }
    
    with open(outpath, 'w') as f:
        for k, v in best_dict.items():
            f.write(f"{k:12s} = {v}\n")

    print(f"Best-fit (model={args.model}, sampler={args.sampler}): {best_dict}")
    print(f"[fit] chi2_total={chi2_total:.2f}  dof={dof}  chi2/dof={chi2_total/max(dof,1):.3f}   AIC={AIC:.2f}  BIC={BIC:.2f}")
    print(f"[timing] total {t1 - t0:.2f} s")
    print(f"Best-fit saved to {outpath}")

    # Diagnostics plots for SN only (resids + Hubble diagram) if requested
    if args.show and plt is not None and extras.get('use_sn', True):
        tag_text            = format_model_tag(args.model, bestfit, k_active, chi2_total, dof, AIC, BIC, frozen_items)

        fig1, ax1           = plt.subplots(figsize=(7.5,4.8))
        
        ax1.scatter(z, resid, s=14, alpha=0.8, label='residuals')
        ax1.axhline(0.0, ls=':', lw=1.2)
        ax1.set_xlabel("z")
        ax1.set_ylabel("Residual (mag)")
        ax1.set_title("Residuals")
        ax1.legend(loc='best')
        fig1.suptitle(tag_text, fontsize=9, y=0.98)

        order = np.argsort(z)
        fig2, ax2 = plt.subplots(figsize=(7.5,4.8))
        ax2.scatter(z, mu, s=10, label='SN data')
        ax2.plot(z[order], mu_th[order], lw=1.5, label=f"best {args.model}")
        ax2.set_xlabel("z")
        ax2.set_ylabel("Î¼ (mag)")
        ax2.set_title("Hubble diagram")
        ax2.legend(loc='best')
        fig2.suptitle(tag_text, fontsize=9, y=0.98)

        plt.tight_layout()
        plt.show()

    # Simple corner/trace if we sampled
    if args.corner and flat_chain.size:
        corner_label_map    = {'H0':'H0','Om':'Î©m','Ol':'Î©Î›','A_acc':'A_acc','n_acc':'n_acc',
                               'gamma_c':'Î³_c','epsilon_M':'Îµ_M',
                               'Omega_chi0':'Omega_chi0','w0_chi':'w0_chi','wa_chi':'wa_chi',
                               'Omega_dark0':'Omega_dark0','su2_fraction':'su2_fraction',
                               'plamb_alpha':'plamb_alpha'}
        if args.model.upper() == 'PLAMB':
            corner_label_map.update({
                'Om':'Omega_inertia', 'A_acc':'clump_amp', 'n_acc':'clump_index',
                'gamma_c':'c_drift', 'epsilon_M':'inertia_drift',
                'plamb_alpha':'dimming_alpha',
            })
        labels              = [corner_label_map[k] for k in active_names]
        
        suptitle            = format_model_tag(args.model, bestfit, k_active, chi2_total, dof, AIC, BIC, frozen_items)
        
        save_corner_and_traces(args.out, flat_chain, labels, suptitle=suptitle)

        print(f"[files] saved corner: {os.path.join(args.out, 'corner.png')}")
        print(f"[files] saved traces: {os.path.join(args.out, 'traces.png')}")


if __name__ == "__main__":
    main()

