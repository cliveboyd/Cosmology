#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 05:48:34 2025

@author: boyde

HHT on SN residuals vs log(1+z), with IMF energy spectrum, 
instantaneous-frequency ridges, 
and significance via phase-randomized surrogates.

What to look for (targets)

SN residual microstructure: IMFs with narrow instantaneous-frequency bands 
in log(1+z) → potential quasi-periodic ripples or step-like structure 
consistent with Λ(z) fluctuations or c(z) features.

Scale locking: phase/frequency of leading IMFs vs. BAO-effective scale; 
alignment across redshift bins.

Coupling signals: correlation between IMF energies/phases and your r
ealized δΛ(z) or PBH γ_holo factor.

Integrator fingerprints: differences in IMF content between Simpson vs discrete-event paths 
(mode splitting or extra high-freq IMFs).

PROGRAM TESTED NO OUTPUT ON CONSOLE UNTIL PLOT

DONE
%run '/Users/boyde/.spyder-py3/hht_sn_residuals.py' \
     '/Users/boyde/.spyder-py3/PLamb_Test_10V6c_plus.py' \
     '/Users/boyde/Downloads/Pantheon+SH0ES.dat' \
     'plamb_runs/hht_sn_FRPBH'

DONE 1) integrator fingerprint test (Simpson vs discrete-events)
# Simpson baseline
%run '/Users/boyde/.spyder-py3/PLamb_Test_10V6c_plus.py' \
     --model FR --data '/Users/boyde/Downloads/Pantheon+SH0ES.dat' \
     --sampler none --out 'plamb_runs/hht_compare/simpson'


DONE
%run '/Users/boyde/.spyder-py3/hht_sn_residuals.py' \
     '/Users/boyde/.spyder-py3/PLamb_Test_10V6c_plus.py' \
     '/Users/boyde/Downloads/Pantheon+SH0ES.dat' \
     'plamb_runs/hht_compare/simpson'

DONE # discrete-events (Δz = 0.01)
%run '/Users/boyde/.spyder-py3/PLamb_Test_10V6c_plus.py' \
     --model FR --data '/Users/boyde/Downloads/Pantheon+SH0ES.dat' \
     --sampler none --discrete-events --event-size 0.01 \
     --out 'plamb_runs/hht_compare/discrete_001'

DONE
%run '/Users/boyde/.spyder-py3/hht_sn_residuals.py' \
     '/Users/boyde/.spyder-py3/PLamb_Test_10V6c_plus.py' \
     '/Users/boyde/Downloads/Pantheon+SH0ES.dat' \
     'plamb_runs/hht_compare/discrete_001'
   
DONE     # realize Λ(z) once (OU) and keep params fixed (sampler none)
%run '/Users/boyde/.spyder-py3/PLamb_Test_10V6c_plus.py' \
     --model FR --data '/Users/boyde/Downloads/Pantheon+SH0ES.dat' \
     --sampler none --fluct-Lambda --Lambda-model OU --Lambda-amp 1e-120 \
     --Lambda-corrz 0.3 --seed-Lambda 777 \
     --out 'plamb_runs/hht_couple/ou_seed777'

DONE 2) Λ(z) coupling / phase-locking probe
# run HHT on those residuals (the script will pick up the model via import)
%run '/Users/boyde/.spyder-py3/hht_sn_residuals.py' \
     '/Users/boyde/.spyder-py3/PLamb_Test_10V6c_plus.py' \
     '/Users/boyde/Downloads/Pantheon+SH0ES.dat' \
     'plamb_runs/hht_couple/ou_seed777'
     

TODO ERROR 1) Turn Λ(z) on (actually couple it) and test phase-locking in HHT
%run '/Users/boyde/.spyder-py3/PLamb_Test_10V6c_plus.py' \
  --model FR \
  --data '/Users/boyde/Downloads/Pantheon+SH0ES.dat' \
  --sampler none \
  --fluct-Lambda --Lambda-model OU --Lambda-amp 1e-120 \
  --Lambda-corrz 0.3 --seed-Lambda 777 \
  --Lambda-coupling --Lambda-gamma 1.0 \
  --out 'plamb_runs/hht_couple/ou_seed777_ON'

DONE
%run '/Users/boyde/.spyder-py3/hht_sn_residuals.py' \
  '/Users/boyde/.spyder-py3/PLamb_Test_10V6c_plus.py' \
  '/Users/boyde/Downloads/Pantheon+SH0ES.dat' \
  'plamb_runs/hht_couple/ou_seed777_ON'


DONE 2) Sensitivity to event discretization (confirm robustness)
%run '/Users/boyde/.spyder-py3/PLamb_Test_10V6c_plus.py' \
  --model FR \
  --data '/Users/boyde/Downloads/Pantheon+SH0ES.dat' \
  --sampler none --discrete-events --event-size 0.005 \
  --out 'plamb_runs/hht_compare/discrete_0005'

DONE
%run '/Users/boyde/.spyder-py3/hht_integrator_ab.py' \
  '/Users/boyde/.spyder-py3/PLamb_Test_10V6c_plus.py' \
  '/Users/boyde/Downloads/Pantheon+SH0ES.dat'

DONE 3) Low-z leverage test (are high-freq IMFs driven by nearby SNe?)
%run '/Users/boyde/.spyder-py3/hht_sn_residuals.py' \
  '/Users/boyde/.spyder-py3/PLamb_Test_10V6c_plus.py' \
  '/Users/boyde/Downloads/Pantheon+SH0ES.dat' \
  'plamb_runs/hht_sn_zmin002' --zmin 0.02

DONE
%run '/Users/boyde/.spyder-py3/hht_sn_residuals.py' \
  '/Users/boyde/.spyder-py3/PLamb_Test_10V6c_plus.py' \
  '/Users/boyde/Downloads/Pantheon+SH0ES.dat' \
  'plamb_runs/hht_sn_zmin010' --zmin 0.10

DONE 4) Noether toggle to see if symmetry suppresses certain IMFs
# Hard symmetry: accretion off
%run '/Users/boyde/.spyder-py3/hht_sn_residuals.py' \
  '/Users/boyde/.spyder-py3/PLamb_Test_10V6c_plus.py' \
  '/Users/boyde/Downloads/Pantheon+SH0ES.dat' \
  'plamb_runs/hht_sn_noether/conserve_time' \
  --noether conserve_time --zmin 0.02

DONE # Soft symmetry combo (gentle pull)
%run '/Users/boyde/.spyder-py3/hht_sn_residuals.py' \
  '/Users/boyde/.spyder-py3/PLamb_Test_10V6c_plus.py' \
  '/Users/boyde/Downloads/Pantheon+SH0ES.dat' \
  'plamb_runs/hht_sn_noether/combo_soft6' \
  --noether combo_soft --noether-strength 6 --zmin 0.02

DONE 5) Significance bands via surrogates (are IMFs > noise?)
%run '/Users/boyde/.spyder-py3/hht_surr_significance.py' \
     '/Users/boyde/.spyder-py3/PLamb_Test_10V6c_plus.py' \
     '/Users/boyde/Downloads/Pantheon+SH0ES.dat' \
     --nsurr 500 --method simpson --out 'plamb_runs/hht_surr'

DONE 6) BAO “locking” probe (does an IMF track a BAO-scale trend?)
%run '/Users/boyde/.spyder-py3/hht_ridge_bao_lock.py' \
     '/Users/boyde/Downloads/bao_long.csv' \
     --rd 147.09 --out 'plamb_runs/hht_bao_lock'




"""

#!/usr/bin/env python3
# HHT on SN residuals (Pantheon+); compares FR/PBH fits; optional surrogates.

import os, sys, json, glob, math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from   sklearn.preprocessing import PolynomialFeatures
from   sklearn.linear_model  import LinearRegression
from   datetime              import datetime

import matplotlib.pyplot     as     plt
plt.ioff()

try:
    from PyEMD import EMD as _EMD
    HAS_EMD                                   = True
except Exception:
    HAS_EMD                                   = False
    
# ---------- optional deps ----------
try:
    from PyEMD import EEMD as _EEMD
    HAS_EEMD                                   = True
except Exception:
    HAS_EEMD                                   = False


try:
    from scipy.signal import hilbert as _hilbert
    HAS_HILBERT                               = True
except Exception:
    HAS_HILBERT                               = False


path = "/Users/boyde/.spyder-py3/skymap2/quaia_outputs/" \
       "quaia_hht_bcut25_z0p30_2p30_apar_simple.txt"

data   = np.loadtxt(path)
z      = data[:, 0]
signal = data[:, 1]   # A_par(z) aligned with CMB

# Optional: model-agnostic de-trend

# ---------- helpers (align =) ----------
def _load_bestfit_txt(path):
    out                                         = {}
    with open(path, 'r') as f:
        for ln in f:
            if '=' not in ln:                   continue
            k, v                                 = ln.split('=', 1)
            k                                    = k.strip()
            v                                    = v.strip()
            try:  out[k]                         = float(v)
            except: out[k]                       = v
    return out

# ---------- add near the top (helpers) ----------
def _get_col_case_insensitive(arr, key_candidates):
    """
    Return array column using a case-insensitive match over dtype.names.
    key_candidates: str or list[str] of possible logical names (e.g., ['zcmb','z'])
    """
    if isinstance(key_candidates, str):
        key_candidates = [key_candidates]
    name_map                                 = {n.lower(): n for n in arr.dtype.names}
    for k in key_candidates:
        if k.lower() in name_map:
            return np.asarray(arr[name_map[k.lower()]], float)
    raise KeyError(f"No matching column for {key_candidates}; available={arr.dtype.names}")


def _instantaneous_freq_phase(signal, x):
    """
    Instantaneous freq over x (monotone, here x = log(1+z)).
    Prefers analytic-signal Hilbert; falls back to finite differences.
    """
    y                                           = np.asarray(signal, float)
    x                                           = np.asarray(x, float)
    dx                                          = np.gradient(x)
    dx                                          = np.where(np.abs(dx) < 1e-9, 1e-9, dx)   # clamp to avoid spikes

    if HAS_HILBERT:
        a                                      = _hilbert(y)
        phase                                  = np.unwrap(np.angle(a))
        inst_freq                              = np.gradient(phase) / np.maximum(dx, 1e-12)
        amp_envelope                           = np.abs(a)
        return inst_freq, phase, amp_envelope
    # fallback: monotone phase via normalized cumulative integral
    y_norm                                      = (y - y.mean()) / (y.std() + 1e-12)
    phase                                       = np.unwrap(np.cumsum(y_norm) * dx / np.max(dx))
    inst_freq                                   = np.gradient(phase) / np.maximum(dx, 1e-12)
    amp_envelope                                = np.abs(y_norm)
    return inst_freq, phase, amp_envelope

def _binned_stat(x, y, nbins=20, min_count=8, fn=np.mean):
    edges                                       = np.linspace(x.min(), x.max(), nbins + 1)
    xc, val                                     = [], []
    for i in range(nbins):
        m                                       = (x >= edges[i]) & (x < edges[i+1])
        if np.count_nonzero(m) >= min_count:
            xc.append(0.5*(edges[i] + edges[i+1]))
            val.append(fn(y[m]))
    return np.array(xc), np.array(val)

def _phase_randomize(y):
    """Fourier-phase randomization surrogate (preserves power)."""
    Y                                           = np.fft.rfft(y)
    mag                                         = np.abs(Y)
    ph                                          = np.angle(Y)
    rnd                                         = np.random.uniform(-np.pi, np.pi, size=ph.shape)
    Ys                                          = mag * np.exp(1j * rnd)
    ys                                          = np.fft.irfft(Ys, n=y.size)
    return ys


def _cbb_surrogate(y, block=48):
    """
    Circular Block Bootstrap surrogate: preserves local correlation structure.
    """
    y                                           = np.asarray(y, float)
    n                                           = y.size
    k                                           = int(np.ceil(n / block))
    starts                                      = np.random.randint(0, n, size=k)
    segs                                        = []
    for s in starts:
        e                                       = s + block
        if e <= n:
            segs.append(y[s:e])
        else:
            wrap                                = (s + block) % n
            segs.append(np.r_[y[s:], y[:wrap]])
    ys                                          = np.concatenate(segs)[:n]
    return ys


# ---------- main ----------
def main(main_py, pantheon_path, cov_path=None, model='PBH', dM_mode='linear_z',
         outdir='plamb_runs/hht_sn', nint=220, n_surr=200, seed=123):

    np.random.seed(int(seed))
    outdir                                      = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)
    print(f"[out] writing to: {outdir}")


    # import your main once
    import importlib.util
    spec                                        = importlib.util.spec_from_file_location("plambmain", main_py)
    mod                                         = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    # load data
    arr                                         = np.genfromtxt(pantheon_path, names=True, dtype=None, encoding=None)
    z                                           = _get_col_case_insensitive(arr, ['zcmb','z'])
    mu                                          = _get_col_case_insensitive(arr, ['mu_sh0es','mu'])

    # pick params (use defaults or bestfit if exists in outdir)
    bf_list                                     = glob.glob(os.path.join(outdir, '*_bestfit.txt'))
    if bf_list:
        best                                    = _load_bestfit_txt(bf_list[0])
        H0                                      = best.get('H0', 70.0)
        Om                                      = best.get('Om', 0.3)
        Ol                                      = best.get('Ol', 0.7)
        A_acc                                   = best.get('A_acc', 0.0)
        n_acc                                   = best.get('n_acc', 0.0)
        gamma_c                                 = best.get('gamma_c', 0.0)
        epsilon_M                               = best.get('epsilon_M', 0.0)
    else:
        H0                                      = 70.0
        Om                                      = 0.3
        Ol                                      = 0.7
        A_acc                                   = 0.0
        n_acc                                   = 0.0
        gamma_c                                 = 0.0
        epsilon_M                               = 0.0

    # compute residuals with your model_mu (Λ(z) etc is already threaded via args_global)
    mu_th                                       = mod.model_mu(z, model, H0, Om, Ol, A_acc, n_acc, gamma_c, epsilon_M, dM_mode, nint)
    resid                                       = mu_th - mu
    
    # --- OPTION A: Polynomial de-trend in redshift z (recommended for model-agnostic trend removal) ---
    # This replaces the Savitzky–Golay block below. Do ONE de-trend method only.
    DO_POLY_DETREND = True     # set False to skip
    
    if DO_POLY_DETREND:
        # Use z or log(1+z) as the regressor. Pick ONE. Using log(1+z) aligns with your HHT axis.
        Z = np.log1p(z).reshape(-1, 1)   # <-- alternative: z.reshape(-1,1)
        y = resid.astype(float)
    
        poly = PolynomialFeatures(5, include_bias=True)  # degree=5 is a good starting point
        X = poly.fit_transform(Z)
        lr = LinearRegression().fit(X, y)
        yhat = lr.predict(X)
        resid = y - yhat
    
        # (optional) keep meta to save into outputs later if you want traceability
        detrend_info = {
            "method": "poly_in_log1pz",
            "degree": 5,
            "coef": lr.coef_.tolist(),
            "intercept": float(lr.intercept_),
        }
    else:
        detrend_info = None

    # gentle low-order trend removal so EMD focuses on microstructure
    if not DO_POLY_DETREND:
        try:
            from scipy.signal import savgol_filter as _savgol
            win = 81 if resid.size > 200 else max(9, (resid.size // 5) * 2 + 1)
            resid_trend = _savgol(resid, win, 3)
            resid = resid - resid_trend
        except Exception:
            pass

    # sort by log(1+z) for quasi-stationarity along the HHT axis
    L                                           = np.log1p(z)
    order                                       = np.argsort(L)
    Ls, rs                                      = L[order], resid[order]

    # ---- robust clip to avoid mode-splitting on outliers ----
    m_resid                                     = np.nanmedian(rs)
    s_resid                                     = np.nanmedian(np.abs(rs - m_resid)) * 1.4826
    good                                        = np.abs(rs - m_resid) <= 3.5 * max(s_resid, 1e-6)
    keep_frac                                   = float(np.count_nonzero(good)) / float(rs.size)
    if keep_frac < 0.6:   # fail-safe: don’t decimate the data
        print(f"[clip] skipped (keep_frac={keep_frac:.2f} < 0.60).")
    else:
        Ls                                      = Ls[good]
        rs                                      = rs[good]
        print(f"[clip] kept {keep_frac:.2f} of points.")

    # ---- resample to a uniform grid in L for Hilbert/EMD stability ----
    N_uni                                       = max(512, int(Ls.size // 1.5))
    Lmin, Lmax                                  = Ls.min(), Ls.max()
    L_uni                                       = np.linspace(Lmin, Lmax, N_uni)

    try:
        from scipy.interpolate import CubicSpline as _CS
        rs_spline                               = _CS(Ls, rs)(L_uni)
    except Exception:
        rs_spline                               = np.interp(L_uni, Ls, rs)

    # ---- cosine taper at edges to reduce boundary artifacts ----
    M_tap                                       = max(12, N_uni // 50)
    taper                                       = np.ones_like(rs_spline)
    w                                           = 0.5 * (1 - np.cos(np.linspace(0, np.pi, M_tap)))
    taper[:M_tap]                                = w
    taper[-M_tap:]                               = w[::-1]
    rs_uni                                      = rs_spline * taper

    # replace working arrays with uniform/tapered versions
    Ls                                          = L_uni
    rs                                          = rs_uni
    print(f"[grid] uniform L bins = {Ls.size}, L ∈ [{Ls.min():.4f},{Ls.max():.4f}]")


    # EMD / IMF decomposition (prefer EEMD to reduce mode mixing)
    if (not HAS_EMD) and (not globals().get('HAS_EEMD', False)):
        print("[warn] PyEMD not found; proceeding with a single-IMF proxy.")
        IMFs                                    = rs[None, :]
    else:
        if globals().get('HAS_EEMD', False):
            try:
                eemd                            = _EEMD(trials=100, noise_width=0.2)
                IMFs                            = eemd.eemd(rs, Ls)
            except Exception:
                emd                             = _EMD()
                IMFs                            = emd.emd(rs)
        else:
            emd                                 = _EMD()
            IMFs                                = emd.emd(rs)

        if IMFs is None or IMFs.size == 0:
            IMFs                                = rs[None, :]

    # instantaneous frequency and amplitude for each IMF
    inst_freqs, envelopes                       = [], []
    for k in range(IMFs.shape[0]):
        f, ph, A                                = _instantaneous_freq_phase(IMFs[k], Ls)
        inst_freqs.append(f)
        envelopes.append(A)
    inst_freqs                                   = np.array(inst_freqs, dtype=float)
    envelopes                                    = np.array(envelopes,  dtype=float)

    # IMF energy spectrum
    energies                                     = np.sum(IMFs**2, axis=1)
    ranks                                        = np.argsort(energies)[::-1]

    # surrogate significance on total IMF energy at peak frequency
    peak_w                                       = []
    for k in ranks[: min(3, IMFs.shape[0])]:
        w_med                                   = np.nanmedian(inst_freqs[k][np.isfinite(inst_freqs[k])])
        peak_w.append(w_med)

    surr_peaks                                   = []
    for s in range(n_surr):
        ys                                      = _phase_randomize(rs) if (s % 2 == 0) else _cbb_surrogate(rs, block=48)

        if globals().get('HAS_EEMD', False):
            try:
                eemd_s                          = _EEMD(trials=50, noise_width=0.2)
                imfs_s                          = eemd_s.eemd(ys, Ls)
            except Exception:
                imfs_s                          = _EMD().emd(ys) if HAS_EMD else ys[None, :]
        elif HAS_EMD:
            imfs_s                              = _EMD().emd(ys)
        else:
            imfs_s                              = ys[None, :]

        if imfs_s is None or imfs_s.size == 0:
            imfs_s                              = ys[None, :]

        e_s                                     = np.sum(imfs_s**2, axis=1).max()
        surr_peaks.append(e_s)

    surr_peaks                                   = np.array(surr_peaks)
    e_max                                        = energies.max()
    p_val                                        = float((np.sum(surr_peaks >= e_max) + 1) / (n_surr + 1))


    # ----- Export IMF stats (with ranks, p_val, model and params) -----
    # ranks[] is indices sorted by energy desc; build inverse map -> rank_of[k]
    rank_of                                     = np.empty(IMFs.shape[0], dtype=int)
    for r, k in enumerate(ranks, 1):
        rank_of[k]                              = r

    imf_stats                                   = []
    for k in range(IMFs.shape[0]):
        f                                       = inst_freqs[k]
        f_med                                   = np.nanmedian(f[np.isfinite(f)]) if np.any(np.isfinite(f)) else np.nan
        imf_stats.append((k+1, energies[k], f_med, int(rank_of[k])))

    imf_stats                                   = np.array(imf_stats,
                                                          dtype=[('IMF','i4'),
                                                                 ('energy','f8'),
                                                                 ('f_med','f8'),
                                                                 ('rank','i4')])

    # --- save detailed HHT arrays for downstream tests (BAO locking, etc.) ---
    detail_npz_path        = os.path.join(outdir, 'hht_sn_detail.npz')
    
    # when saving the detailed NPZ:
    np.savez_compressed(detail_npz_path,
        Ls=Ls, IMFs=IMFs, inst_freqs=inst_freqs, envelopes=envelopes,
        energies=energies, ranks=ranks,
        detrend_info=json.dumps(detrend_info) if detrend_info is not None else "none"
    )

    print(f"[save] HHT detail -> {detail_npz_path}")


    csv_path                                    = os.path.join(outdir, 'hht_sn_imf_stats.csv')
    npz_path                                    = os.path.join(outdir, 'hht_sn_imf_stats.npz')

    import csv
    with open(csv_path, 'w', newline='') as f:
        writer                                  = csv.writer(f)
        writer.writerow(['IMF','energy','median_freq','rank','p_val',
                         'model','H0','Om','Ol','A_acc','n_acc','gamma_c','epsilon_M'])
        for row in imf_stats:
            writer.writerow([row['IMF'], f"{row['energy']:.10g}", f"{row['f_med']:.10g}", row['rank'],
                             f"{p_val:.6g}", model, f"{H0:.10g}", f"{Om:.10g}", f"{Ol:.10g}",
                             f"{A_acc:.10g}", f"{n_acc:.10g}", f"{gamma_c:.10g}", f"{epsilon_M:.10g}"])

    np.savez_compressed(npz_path,
        imf_stats = imf_stats,
        p_val     = p_val,
        model     = model,
        H0        = H0, Om = Om, Ol = Ol,
        A_acc     = A_acc, n_acc = n_acc,
        gamma_c   = gamma_c, epsilon_M = epsilon_M
    )

    print(f"[save] IMF stats -> {csv_path}, {npz_path}")


    run_id                                     = datetime.now().strftime("%Y%m%d-%H%M%S")
    np.savez_compressed(npz_path,
        imf_stats = imf_stats,
        p_val     = p_val,
        model     = model,
        H0        = H0, Om = Om, Ol = Ol,
        A_acc     = A_acc, n_acc = n_acc,
        gamma_c   = gamma_c, epsilon_M = epsilon_M,
        run_id    = run_id
    )

    # ------- plots -------
    fig, axs                                     = plt.subplots(2, 2, figsize=(11.5, 7.8))
    axs[0,0].scatter(Ls, rs, s=8, alpha=0.7)
    axs[0,0].axhline(0, ls=':', lw=1)
    axs[0,0].set_xlabel('log(1+z)'); axs[0,0].set_ylabel('μ residual (mag)')
    axs[0,0].set_title('SN residuals (sorted)')

    # IMF overlay
    off                                          = 0.0
    for k in ranks:
        axs[0,1].plot(Ls, IMFs[k] + off, lw=1.0, label=f'IMF{k+1}')
        off                                     += 0.07
    axs[0,1].set_xlabel('log(1+z)'); axs[0,1].set_ylabel('IMFs (offset)')
    axs[0,1].set_title('EMD decomposition')
    axs[0,1].legend(fontsize=8, ncol=2)

    # Energy spectrum + surrogate
    axs[1,0].bar(np.arange(1, IMFs.shape[0]+1), energies, width=0.7)
    axs[1,0].axhline(np.median(surr_peaks), ls='--', lw=1.0, label='surrogate median')
    axs[1,0].set_xlabel('IMF index'); axs[1,0].set_ylabel('energy')
    axs[1,0].set_title(f'IMF energy spectrum (p≈{p_val:.3f})')
    axs[1,0].legend()

    # Instantaneous frequency ridge: show the top IMF
    # Instantaneous frequency ridges: overlay top K IMFs
    if IMFs.shape[0] > 0:
        K                                       = min(3, IMFs.shape[0])
    
        for j, k in enumerate(ranks[:K], 1):
            f                                   = inst_freqs[k]
            m                                   = np.isfinite(f)
    
            # optional light smoothing for readability (median filter)
            try:
                from scipy.signal import medfilt as _medfilt
                ks                              = 9 if np.count_nonzero(m) >= 9 else max(3, (np.count_nonzero(m) | 1))
                f_plot                          = _medfilt(f[m], kernel_size=ks)
            except Exception:
                f_plot                          = f[m]
    
            axs[1,1].plot(Ls[m], f_plot, lw=1.2, label=f'IMF{k+1}')
    
        # robust autoscale to avoid extreme spikes dominating the view
        try:
            y_all                               = np.concatenate([inst_freqs[k][np.isfinite(inst_freqs[k])] for k in ranks[:K]])
            if y_all.size:
                lo, hi                          = np.percentile(y_all, [1, 99])
                pad                              = 0.05 * (hi - lo) if hi > lo else 1.0
                axs[1,1].set_ylim(lo - pad, hi + pad)
        except Exception:
            pass
    
        axs[1,1].set_xlabel('log(1+z)'); axs[1,1].set_ylabel('ω_inst')
        axs[1,1].set_title(f'Instantaneous frequency (top {K} IMFs)')
        axs[1,1].legend(fontsize=8, ncol=K)
    
    
    fig.tight_layout()
    png_path                                     = os.path.join(outdir, 'hht_sn.png')
    fig.savefig(png_path, dpi=160, bbox_inches='tight')
    print(f"[files] {png_path}")
    plt.close(fig)

    # dump summary
    np.savez_compressed(os.path.join(outdir, 'hht_sn_summary.npz'),
        Ls = Ls, resid = rs, energies = energies, ranks = ranks,
        peak_w = np.array(peak_w), p_val = p_val
    )

    plt.close('all')


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python hht_sn_residuals.py <path_to_10V6_plus.py> <Pantheon+ dat> [outdir]")
        sys.exit(1)
    main_py       = sys.argv[1]
    pantheon_path = sys.argv[2]
    outdir        = sys.argv[3] if len(sys.argv) > 3 else 'plamb_runs/hht_sn'
    main(main_py, pantheon_path, None, outdir=outdir)
