#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  9 19:48:23 2025

@author: boyde
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# === File Paths ===
base_path = '/Users/boyde/Downloads'
tt_file = os.path.join(base_path, 'COM_PowerSpect_CMB-TT-full_R3.01.txt')
te_file = os.path.join(base_path, 'COM_PowerSpect_CMB-TE-full_R3.01.txt')
ee_file = os.path.join(base_path, 'COM_PowerSpect_CMB-EE-full_R3.01.txt')

# === Load Planck CMB Data ===
def load_cmb_data(filepath):
    try:
        data = np.loadtxt(filepath, comments="#")
        ell = data[:, 0]
        Dl = data[:, 1]
        dDl_low = data[:, 2]
        dDl_high = data[:, 3]
        yerr = np.abs([dDl_low, dDl_high])  # Ensure errors are positive
        return ell, Dl, yerr
    except Exception as e:
        print(f"Failed to load {filepath}: {e}")
        return None, None, None

ell_tt, Dl_tt, yerr_tt = load_cmb_data(tt_file)
ell_te, Dl_te, yerr_te = load_cmb_data(te_file)
ell_ee, Dl_ee, yerr_ee = load_cmb_data(ee_file)

# === Plot the Spectra ===
plt.figure(figsize=(12, 7))
if ell_tt is not None:
    plt.errorbar(ell_tt, Dl_tt, yerr=yerr_tt, fmt='o', ms=2, alpha=0.6, label='Planck 2018 TT')
if ell_te is not None:
    plt.errorbar(ell_te, Dl_te, yerr=yerr_te, fmt='^', ms=2, alpha=0.6, label='Planck 2018 TE')
if ell_ee is not None:
    plt.errorbar(ell_ee, Dl_ee, yerr=yerr_ee, fmt='s', ms=2, alpha=0.6, label='Planck 2018 EE')

plt.xlabel('Multipole $\\ell$')
plt.ylabel('$D_\\ell = \\ell(\\ell + 1)C_\\ell / 2\\pi$ [$\\mu K^2$]')
plt.title('Planck 2018 CMB Power Spectra (TT, TE, EE)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

