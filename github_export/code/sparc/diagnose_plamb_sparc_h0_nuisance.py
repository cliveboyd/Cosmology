#!/usr/bin/env python3
r"""
Extended SPARC rotation audit for PLAMB development.

This script targets three investigations:

1. Test PLAMB-like functional forms that can approach the empirical RAR.
2. Scan SPARC Hubble-flow distance rescaling over an H0 range.
3. Add a controlled per-galaxy stellar mass-to-light nuisance check.

The distance rescaling is deliberately explicit. For rows with SPARC distance
method f_D=1, the sample documentation says Hubble-flow distances used H0=73.
For a trial H0, the script rescales those distances as

    D_new = D_old * (H0_ref / H0_trial)

and correspondingly rescales radius and baryonic model velocities:

    R_new = R_old * scale
    V_baryon_component,new = V_old * sqrt(scale)

Observed Doppler velocities Vobs are not rescaled.

Run:
    python diagnose_plamb_sparc_h0_nuisance.py
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import replace
from pathlib import Path
from typing import Any, Callable

import numpy as np

try:
    from scipy.optimize import differential_evolution, minimize, minimize_scalar
except Exception as exc:  # pragma: no cover - reported at runtime
    differential_evolution = None
    minimize = None
    minimize_scalar = None
    SCIPY_IMPORT_ERROR = exc
else:
    SCIPY_IMPORT_ERROR = None

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as exc:  # pragma: no cover - optional plotting
    plt = None
    MATPLOTLIB_IMPORT_ERROR = exc
else:
    MATPLOTLIB_IMPORT_ERROR = None

from diagnose_plamb_sparc_rotation import (
    BASE_OUTDIR as SPARC_BASE_OUTDIR,
    DEFAULT_G0,
    DEFAULT_POINTS,
    DEFAULT_SAMPLE,
KPC_M,
    Dataset,
    acceleration_from_v2,
    baryon_v2,
    load_dataset,
    velocity_from_acceleration,
)


ROOT = Path(__file__).resolve().parents[3]
OUTDIR = ROOT / "plamb_runs" / "diagnostics" / "sparc_rotation_h0_nuisance"
SPARC_H0_REF = 73.0
C_LIGHT_M_S = 299_792_458.0
MPC_M = 3.085677581491367e22


class ModelSpec:
    def __init__(
        self,
        name: str,
        param_names: list[str],
        bounds: list[tuple[float, float]],
        initial: list[float],
        predictor: Callable[[np.ndarray, Dataset], np.ndarray],
        note: str,
        has_stellar_ml: bool = True,
    ) -> None:
        self.name = name
        self.param_names = param_names
        self.bounds = bounds
        self.initial = initial
        self.predictor = predictor
        self.note = note
        self.has_stellar_ml = has_stellar_ml


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_h0_values(value: str) -> list[float]:
    values: list[float] = []
    for token in value.split(","):
        token = token.strip()
        if not token:
            continue
        values.append(float(token))
    if not values:
        raise ValueError("No H0 values were supplied.")
    return values


def filter_model_specs(specs: list[ModelSpec], model_filter: str) -> list[ModelSpec]:
    requested = [token.strip() for token in model_filter.split(",") if token.strip()]
    if not requested or any(token.lower() == "all" for token in requested):
        return specs

    by_name = {spec.name: spec for spec in specs}
    missing = [name for name in requested if name not in by_name]
    if missing:
        available = ", ".join(sorted(by_name))
        raise ValueError(f"Unknown model(s): {', '.join(missing)}. Available models: {available}")
    return [by_name[name] for name in requested]


def rescale_for_h0(data: Dataset, h0: float, h0_ref: float, mode: str) -> Dataset:
    if mode == "none":
        scale = np.ones_like(data.rad_kpc)
    elif mode == "all":
        scale = np.full_like(data.rad_kpc, h0_ref / h0, dtype=float)
    elif mode == "hubble_flow_only":
        scale = np.where(data.distance_method == 1, h0_ref / h0, 1.0)
    else:
        raise ValueError(f"Unknown distance scaling mode: {mode}")

    velocity_scale = np.sqrt(scale)
    return replace(
        data,
        rad_kpc=data.rad_kpc * scale,
        vgas=data.vgas * velocity_scale,
        vdisk=data.vdisk * velocity_scale,
        vbul=data.vbul * velocity_scale,
        d_mpc=data.d_mpc * scale,
    )


def fixed_baryon(theta: np.ndarray, data: Dataset) -> np.ndarray:
    del theta
    return np.sqrt(np.maximum(baryon_v2(data, 0.5, 0.7), 0.0))


def free_baryon(theta: np.ndarray, data: Dataset) -> np.ndarray:
    ydisk, ybul = theta
    return np.sqrt(np.maximum(baryon_v2(data, ydisk, ybul), 0.0))


def rar_free(theta: np.ndarray, data: Dataset) -> np.ndarray:
    ydisk, ybul, log10_gdagger = theta
    gbar = acceleration_from_v2(baryon_v2(data, ydisk, ybul), data.rad_kpc)
    gdagger = 10.0**log10_gdagger
    x = np.sqrt(np.maximum(gbar / gdagger, 1e-30))
    gpred = gbar / np.maximum(1.0 - np.exp(-x), 1e-12)
    return velocity_from_acceleration(gpred, data.rad_kpc)


def plamb_power_free_g0(theta: np.ndarray, data: Dataset) -> np.ndarray:
    ydisk, ybul, amp, index, log10_g0 = theta
    gbar = acceleration_from_v2(baryon_v2(data, ydisk, ybul), data.rad_kpc)
    g0 = 10.0**log10_g0
    ratio = np.clip(g0 / np.maximum(gbar, 1e-30), 1e-8, 1e8)
    gpred = gbar * (1.0 + amp * ratio**index)
    return velocity_from_acceleration(gpred, data.rad_kpc)


def plamb_exp_bridge(theta: np.ndarray, data: Dataset) -> np.ndarray:
    ydisk, ybul, log10_g0, exponent = theta
    gbar = acceleration_from_v2(baryon_v2(data, ydisk, ybul), data.rad_kpc)
    g0 = 10.0**log10_g0
    x = np.maximum(gbar / g0, 1e-30) ** exponent
    gpred = gbar / np.maximum(1.0 - np.exp(-x), 1e-12)
    return velocity_from_acceleration(gpred, data.rad_kpc)


def plamb_optical_depth_cosmo_g0(theta: np.ndarray, data: Dataset, g0: float) -> np.ndarray:
    ydisk, ybul = theta
    gbar = acceleration_from_v2(baryon_v2(data, ydisk, ybul), data.rad_kpc)
    tau = np.sqrt(np.maximum(gbar / g0, 1e-30))
    gpred = gbar / np.maximum(1.0 - np.exp(-tau), 1e-12)
    return velocity_from_acceleration(gpred, data.rad_kpc)


def plamb_optical_depth_kappa(theta: np.ndarray, data: Dataset, g0_cosmo: float) -> np.ndarray:
    ydisk, ybul, log10_kappa = theta
    gbar = acceleration_from_v2(baryon_v2(data, ydisk, ybul), data.rad_kpc)
    g0 = (10.0**log10_kappa) * g0_cosmo
    tau = np.sqrt(np.maximum(gbar / g0, 1e-30))
    gpred = gbar / np.maximum(1.0 - np.exp(-tau), 1e-12)
    return velocity_from_acceleration(gpred, data.rad_kpc)


def plamb_optical_depth_kappa_p(theta: np.ndarray, data: Dataset, g0_cosmo: float) -> np.ndarray:
    ydisk, ybul, log10_kappa, exponent = theta
    gbar = acceleration_from_v2(baryon_v2(data, ydisk, ybul), data.rad_kpc)
    g0 = (10.0**log10_kappa) * g0_cosmo
    tau = np.maximum(gbar / g0, 1e-30) ** exponent
    gpred = gbar / np.maximum(1.0 - np.exp(-tau), 1e-12)
    return velocity_from_acceleration(gpred, data.rad_kpc)


def plamb_screened_power(theta: np.ndarray, data: Dataset) -> np.ndarray:
    ydisk, ybul, amp, index, screen, log10_g0 = theta
    gbar = acceleration_from_v2(baryon_v2(data, ydisk, ybul), data.rad_kpc)
    g0 = 10.0**log10_g0
    ratio = np.clip(g0 / np.maximum(gbar, 1e-30), 1e-8, 1e8)
    screening = 1.0 + np.clip(1.0 / ratio, 1e-8, 1e8) ** screen
    gpred = gbar * (1.0 + amp * ratio**index / screening)
    return velocity_from_acceleration(gpred, data.rad_kpc)


def acceleration_cH0_over_2pi(h0_km_s_mpc: float) -> float:
    h0_si = h0_km_s_mpc * 1000.0 / MPC_M
    return C_LIGHT_M_S * h0_si / (2.0 * math.pi)


def plamb_noether_screened_fixed(theta: np.ndarray, data: Dataset, g0: float) -> np.ndarray:
    ydisk, ybul = theta
    gbar = acceleration_from_v2(baryon_v2(data, ydisk, ybul), data.rad_kpc)
    ratio = np.clip(g0 / np.maximum(gbar, 1e-30), 1e-8, 1e8)
    # Fixed physics test: sqrt low-acceleration inertia term, screened in the
    # high-acceleration limit, with g0 tied to cH0/2pi and no free amplitude.
    gpred = gbar * (1.0 + ratio**0.5 / (1.0 + (1.0 / ratio)))
    return velocity_from_acceleration(gpred, data.rad_kpc)


def plamb_noether_screened_free_screen(theta: np.ndarray, data: Dataset, g0: float) -> np.ndarray:
    ydisk, ybul, screen = theta
    gbar = acceleration_from_v2(baryon_v2(data, ydisk, ybul), data.rad_kpc)
    ratio = np.clip(g0 / np.maximum(gbar, 1e-30), 1e-8, 1e8)
    screening = 1.0 + np.clip(1.0 / ratio, 1e-8, 1e8) ** screen
    gpred = gbar * (1.0 + ratio**0.5 / screening)
    return velocity_from_acceleration(gpred, data.rad_kpc)


def plamb_noether_screened_free_shape(theta: np.ndarray, data: Dataset, g0: float) -> np.ndarray:
    ydisk, ybul, index, screen = theta
    gbar = acceleration_from_v2(baryon_v2(data, ydisk, ybul), data.rad_kpc)
    ratio = np.clip(g0 / np.maximum(gbar, 1e-30), 1e-8, 1e8)
    screening = 1.0 + np.clip(1.0 / ratio, 1e-8, 1e8) ** screen
    gpred = gbar * (1.0 + ratio**index / screening)
    return velocity_from_acceleration(gpred, data.rad_kpc)


def plamb_log_inertia(theta: np.ndarray, data: Dataset) -> np.ndarray:
    ydisk, ybul, amp, index, log10_g0 = theta
    gbar = acceleration_from_v2(baryon_v2(data, ydisk, ybul), data.rad_kpc)
    g0 = 10.0**log10_g0
    ratio = np.clip(g0 / np.maximum(gbar, 1e-30), 1e-8, 1e8)
    gpred = gbar * (1.0 + amp * np.log1p(ratio**index))
    return velocity_from_acceleration(gpred, data.rad_kpc)


def model_specs(h0_trial: float) -> list[ModelSpec]:
    g0_cosmo = acceleration_cH0_over_2pi(h0_trial)

    def noether_fixed(theta: np.ndarray, data: Dataset) -> np.ndarray:
        return plamb_noether_screened_fixed(theta, data, g0_cosmo)

    def optical_depth_cosmo_g0(theta: np.ndarray, data: Dataset) -> np.ndarray:
        return plamb_optical_depth_cosmo_g0(theta, data, g0_cosmo)

    def optical_depth_kappa(theta: np.ndarray, data: Dataset) -> np.ndarray:
        return plamb_optical_depth_kappa(theta, data, g0_cosmo)

    def optical_depth_kappa_p(theta: np.ndarray, data: Dataset) -> np.ndarray:
        return plamb_optical_depth_kappa_p(theta, data, g0_cosmo)

    def noether_free_screen(theta: np.ndarray, data: Dataset) -> np.ndarray:
        return plamb_noether_screened_free_screen(theta, data, g0_cosmo)

    def noether_free_shape(theta: np.ndarray, data: Dataset) -> np.ndarray:
        return plamb_noether_screened_free_shape(theta, data, g0_cosmo)

    return [
        ModelSpec(
            "BARYON_FREE_ML",
            ["ydisk", "ybul"],
            [(0.05, 1.25), (0.05, 1.60)],
            [0.5, 0.7],
            free_baryon,
            "Newtonian baryons with global disk/bulge M/L.",
        ),
        ModelSpec(
            "RAR_FREE_GDAGGER",
            ["ydisk", "ybul", "log10_gdagger"],
            [(0.05, 1.25), (0.05, 1.60), (-11.6, -9.2)],
            [0.5, 0.7, math.log10(DEFAULT_G0)],
            rar_free,
            "Empirical RAR benchmark with fitted acceleration scale.",
        ),
        ModelSpec(
            "PLAMB_POWER_FREE_G0",
            ["ydisk", "ybul", "inertia_amp", "inertia_index", "log10_g0"],
            [(0.05, 1.25), (0.05, 1.60), (0.0, 4.0), (0.0, 1.4), (-11.8, -9.1)],
            [0.45, 0.65, 0.9, 0.5, math.log10(DEFAULT_G0)],
            plamb_power_free_g0,
            "PLAMB low-acceleration inertia power law with fitted pivot.",
        ),
        ModelSpec(
            "PLAMB_EXP_BRIDGE",
            ["ydisk", "ybul", "log10_g0", "bridge_exponent"],
            [(0.05, 1.25), (0.05, 1.60), (-11.8, -9.1), (0.15, 1.20)],
            [0.5, 0.7, math.log10(DEFAULT_G0), 0.5],
            plamb_exp_bridge,
            "PLAMB/RAR bridge: g_bar/[1-exp(-(g_bar/g0)^p)]. RAR is p=0.5.",
        ),
        ModelSpec(
            "PLAMB_OPTICAL_DEPTH_COSMO_G0",
            ["ydisk", "ybul"],
            [(0.05, 1.25), (0.05, 1.60)],
            [0.5, 0.7],
            optical_depth_cosmo_g0,
            "Derived PLAMB/RAR test: tau=sqrt(g_bar/g0), g0=cH0/2pi, fit only global stellar M/L.",
        ),
        ModelSpec(
            "PLAMB_OPTICAL_DEPTH_KAPPA",
            ["ydisk", "ybul", "log10_kappa"],
            [(0.05, 1.25), (0.05, 1.60), (-0.50, 0.50)],
            [0.5, 0.7, 0.0],
            optical_depth_kappa,
            "Minimal deformation: tau=sqrt(g_bar/g0), g0=kappa*cH0/2pi, fit global stellar M/L and kappa.",
        ),
        ModelSpec(
            "PLAMB_OPTICAL_DEPTH_KAPPA_P",
            ["ydisk", "ybul", "log10_kappa", "bridge_exponent"],
            [(0.05, 1.25), (0.05, 1.60), (-0.50, 0.50), (0.15, 1.20)],
            [0.5, 0.7, 0.0, 0.5],
            optical_depth_kappa_p,
            "Minimal deformation: tau=(g_bar/g0)^p, g0=kappa*cH0/2pi, fit global stellar M/L, kappa, and p.",
        ),
        ModelSpec(
            "PLAMB_SCREENED_POWER",
            ["ydisk", "ybul", "inertia_amp", "inertia_index", "screen_index", "log10_g0"],
            [(0.05, 1.25), (0.05, 1.60), (0.0, 6.0), (0.0, 1.6), (0.0, 3.0), (-11.8, -9.1)],
            [0.45, 0.65, 1.0, 0.55, 1.0, math.log10(DEFAULT_G0)],
            plamb_screened_power,
            "PLAMB power boost with high-acceleration screening.",
        ),
        ModelSpec(
            "PLAMB_NOETHER_SCREENED_FIXED",
            ["ydisk", "ybul"],
            [(0.05, 1.25), (0.05, 1.60)],
            [0.5, 0.7],
            noether_fixed,
            "Constrained PLAMB screen: g0=cH0/2pi, amp=1, inertia_index=0.5, screen_index=1.",
        ),
        ModelSpec(
            "PLAMB_NOETHER_SCREENED_FREE_SCREEN",
            ["ydisk", "ybul", "screen_index"],
            [(0.05, 1.25), (0.05, 1.60), (0.0, 3.0)],
            [0.5, 0.7, 1.0],
            noether_free_screen,
            "Constrained PLAMB screen: g0=cH0/2pi, amp=1, inertia_index=0.5, screen fitted.",
        ),
        ModelSpec(
            "PLAMB_NOETHER_SCREENED_FREE_SHAPE",
            ["ydisk", "ybul", "inertia_index", "screen_index"],
            [(0.05, 1.25), (0.05, 1.60), (0.0, 1.4), (0.0, 3.0)],
            [0.5, 0.7, 0.5, 1.0],
            noether_free_shape,
            "Constrained PLAMB screen: g0=cH0/2pi, amp=1, inertia/screen exponents fitted.",
        ),
        ModelSpec(
            "PLAMB_LOG_INERTIA",
            ["ydisk", "ybul", "inertia_amp", "inertia_index", "log10_g0"],
            [(0.05, 1.25), (0.05, 1.60), (0.0, 8.0), (0.0, 2.0), (-11.8, -9.1)],
            [0.45, 0.65, 1.0, 0.7, math.log10(DEFAULT_G0)],
            plamb_log_inertia,
            "PLAMB log-response branch: g_bar*[1 + A log(1+(g0/g_bar)^n)].",
        ),
    ]


def chi2_for(theta: np.ndarray, spec: ModelSpec, data: Dataset) -> float:
    try:
        pred = spec.predictor(theta, data)
    except (FloatingPointError, ValueError, OverflowError):
        return 1.0e100
    if pred.shape != data.vobs.shape or not np.all(np.isfinite(pred)):
        return 1.0e100
    residual = (data.vobs - pred) / data.sigma_v
    value = float(np.sum(residual**2))
    return value if math.isfinite(value) else 1.0e100


def fit_model(spec: ModelSpec, data: Dataset, seed: int, maxiter: int) -> tuple[np.ndarray, float]:
    if differential_evolution is None or minimize is None:
        raise RuntimeError(f"SciPy optimizer import failed: {SCIPY_IMPORT_ERROR}")
    objective = lambda theta: chi2_for(np.asarray(theta, dtype=float), spec, data)
    result = differential_evolution(
        objective,
        spec.bounds,
        seed=seed,
        tol=1e-6,
        maxiter=maxiter,
        popsize=10,
        polish=False,
        updating="immediate",
        workers=1,
    )
    local = minimize(
        objective,
        result.x,
        method="L-BFGS-B",
        bounds=spec.bounds,
        options={"maxiter": 1600, "ftol": 1e-9},
    )
    if local.success and float(local.fun) <= float(result.fun):
        return np.asarray(local.x, dtype=float), float(local.fun)
    return np.asarray(result.x, dtype=float), float(result.fun)


def summarize_fit(
    h0: float,
    mode: str,
    spec: ModelSpec,
    theta: np.ndarray,
    chi2: float,
    n: int,
    n_galaxies: int,
    row_type: str = "global",
    chi2_data: float | None = None,
    chi2_prior: float = 0.0,
    extra_k: int = 0,
) -> dict[str, Any]:
    k = len(theta) + extra_k
    aic = chi2 + 2.0 * k
    bic = chi2 + k * math.log(max(n, 1))
    row: dict[str, Any] = {
        "row_type": row_type,
        "H0_trial": h0,
        "distance_scaling_mode": mode,
        "model": spec.name,
        "N_points": n,
        "N_galaxies": n_galaxies,
        "k": k,
        "chi2": chi2,
        "chi2_data": chi2 if chi2_data is None else chi2_data,
        "chi2_prior": chi2_prior,
        "dof": n - k,
        "chi2_dof": chi2 / max(n - k, 1),
        "AIC": aic,
        "BIC": bic,
        "log10_g0_cosmo": math.log10(acceleration_cH0_over_2pi(h0)),
        "note": spec.note,
    }
    for name, value in zip(spec.param_names, theta):
        row[name] = float(value)
    return row


def subset_dataset(data: Dataset, mask: np.ndarray) -> Dataset:
    return replace(
        data,
        galaxy=data.galaxy[mask],
        rad_kpc=data.rad_kpc[mask],
        vobs=data.vobs[mask],
        sigma_v=data.sigma_v[mask],
        vgas=data.vgas[mask],
        vdisk=data.vdisk[mask],
        vbul=data.vbul[mask],
        quality=data.quality[mask],
        distance_method=data.distance_method[mask],
        inc_deg=data.inc_deg[mask],
        d_mpc=data.d_mpc[mask],
    )


def apply_stellar_ml_eta(theta: np.ndarray, eta: float) -> np.ndarray:
    shifted = np.array(theta, dtype=float, copy=True)
    shifted[0] *= eta
    shifted[1] *= eta
    return shifted


def galaxy_ml_nuisance(
    spec: ModelSpec,
    theta: np.ndarray,
    data: Dataset,
    sigma_ln_eta: float,
    eta_min: float,
    eta_max: float,
) -> tuple[float, float, list[dict[str, Any]]]:
    if minimize_scalar is None:
        raise RuntimeError(f"SciPy optimizer import failed: {SCIPY_IMPORT_ERROR}")
    if not spec.has_stellar_ml or len(theta) < 2:
        return chi2_for(theta, spec, data), 0.0, []

    total_data = 0.0
    total_prior = 0.0
    detail_rows: list[dict[str, Any]] = []
    low = math.log(eta_min)
    high = math.log(eta_max)
    galaxies = sorted(set(data.galaxy.tolist()))
    for galaxy in galaxies:
        mask = data.galaxy == galaxy
        sub = subset_dataset(data, mask)

        def objective(log_eta: float) -> float:
            eta = math.exp(log_eta)
            theta_eta = apply_stellar_ml_eta(theta, eta)
            data_chi2 = chi2_for(theta_eta, spec, sub)
            prior = (log_eta / sigma_ln_eta) ** 2 if sigma_ln_eta > 0.0 else 0.0
            return data_chi2 + prior

        result = minimize_scalar(objective, bounds=(low, high), method="bounded", options={"xatol": 1e-4})
        log_eta = float(result.x)
        eta = math.exp(log_eta)
        theta_eta = apply_stellar_ml_eta(theta, eta)
        data_chi2 = chi2_for(theta_eta, spec, sub)
        prior = (log_eta / sigma_ln_eta) ** 2 if sigma_ln_eta > 0.0 else 0.0
        total_data += data_chi2
        total_prior += prior
        detail_rows.append(
            {
                "model": spec.name,
                "Galaxy": galaxy,
                "eta_star_ml": eta,
                "log_eta": log_eta,
                "chi2_data": data_chi2,
                "chi2_prior": prior,
                "N_points": int(np.sum(mask)),
                "Q": int(sub.quality[0]),
                "f_D": int(sub.distance_method[0]),
            }
        )
    return total_data, total_prior, detail_rows


def make_h0_plot(rows: list[dict[str, Any]], path: Path) -> None:
    if plt is None:
        return
    global_rows = [row for row in rows if row["row_type"] == "global"]
    if not global_rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9.2, 6.2))
    for model in sorted({str(row["model"]) for row in global_rows}):
        model_rows = sorted([row for row in global_rows if row["model"] == model], key=lambda row: float(row["H0_trial"]))
        h0 = [float(row["H0_trial"]) for row in model_rows]
        aic = [float(row["AIC"]) for row in model_rows]
        ax.plot(h0, aic, marker="o", linewidth=1.6, markersize=3.8, label=model)
    ax.set_xlabel("Trial H0 used to rescale SPARC Hubble-flow distances")
    ax.set_ylabel("AIC")
    ax.set_title("SPARC H0 Distance-Rescaling Scan")
    ax.grid(True, alpha=0.22)
    ax.legend(fontsize=7, loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def write_report(
    fit_rows: list[dict[str, Any]],
    nuisance_rows: list[dict[str, Any]],
    args: argparse.Namespace,
    n_points: int,
    n_galaxies: int,
) -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    global_rows = [row for row in fit_rows if row["row_type"] == "global"]
    nuisance_summary = [row for row in fit_rows if row["row_type"] == "galaxy_ml_nuisance"]
    best_global_aic = min(global_rows, key=lambda row: float(row["AIC"]))
    best_global_bic = min(global_rows, key=lambda row: float(row["BIC"]))
    best_nuisance_data = min(nuisance_summary, key=lambda row: float(row["chi2_data"])) if nuisance_summary else None

    lines = [
        "# Extended SPARC PLAMB/RAR/H0/Nuisance Audit",
        "",
        "## Executive Summary",
        "",
        f"- Points fitted per scan: {n_points}",
        f"- Galaxies fitted per scan: {n_galaxies}",
        f"- H0 grid: {args.h0_values}",
        f"- Distance scaling mode: {args.distance_scaling_mode}",
        f"- Best global AIC row: {best_global_aic['model']} at H0={float(best_global_aic['H0_trial']):.3f}, AIC={float(best_global_aic['AIC']):.3f}, chi2/dof={float(best_global_aic['chi2_dof']):.3f}",
        f"- Best global BIC row: {best_global_bic['model']} at H0={float(best_global_bic['H0_trial']):.3f}, BIC={float(best_global_bic['BIC']):.3f}, chi2/dof={float(best_global_bic['chi2_dof']):.3f}",
    ]
    if best_nuisance_data:
        lines.append(
            f"- Lowest data chi2 after per-galaxy M/L nuisance: {best_nuisance_data['model']} "
            f"at H0={float(best_nuisance_data['H0_trial']):.3f}, "
            f"chi2_data={float(best_nuisance_data['chi2_data']):.3f}"
        )
    lines.extend(
        [
            "",
            "## PLAMB Functional Tests",
            "",
            "- `PLAMB_POWER_FREE_G0`: low-acceleration power-law inertia boost with fitted pivot.",
            "- `PLAMB_EXP_BRIDGE`: smooth exponential bridge. If its exponent is near `0.5`, it is effectively reproducing the empirical RAR form.",
            "- `PLAMB_OPTICAL_DEPTH_COSMO_G0`: derived optical-depth branch with `tau=sqrt(g_bar/g0)` and `g0=cH0/2pi`; only global stellar M/L is fitted.",
            "- `PLAMB_OPTICAL_DEPTH_KAPPA`: same square-root optical depth, but allows `g0=kappa*cH0/2pi`.",
            "- `PLAMB_OPTICAL_DEPTH_KAPPA_P`: allows both `g0=kappa*cH0/2pi` and `tau=(g_bar/g0)^p`.",
            "- `PLAMB_SCREENED_POWER`: power-law inertia boost with high-acceleration screening.",
            "- `PLAMB_NOETHER_SCREENED_FIXED`: fixed physics check with `g0=cH0/2pi`, unit amplitude, square-root inertia index, and fixed screening.",
            "- `PLAMB_NOETHER_SCREENED_FREE_SCREEN`: same cosmological `g0`, but lets the high-acceleration screen vary.",
            "- `PLAMB_NOETHER_SCREENED_FREE_SHAPE`: same cosmological `g0` and unit amplitude, but lets the inertia/screen exponents vary.",
            "- `PLAMB_LOG_INERTIA`: logarithmic inertia response.",
            "",
            "## Distance Rescaling",
            "",
            "For SPARC rows with `f_D=1`, the original sample table notes Hubble-flow distances using H0=73. This scan applies:",
            "",
            "`D_new = D_old * (73 / H0_trial)`",
            "",
            "and rescales radii and baryonic model velocities while leaving observed Doppler velocities fixed.",
            "",
            "## Top Global Rows By AIC",
            "",
            "| H0 | Model | chi2/dof | AIC | BIC | Key parameters |",
            "|---:|---|---:|---:|---:|---|",
        ]
    )
    for row in sorted(global_rows, key=lambda item: float(item["AIC"]))[:15]:
        pieces = []
        for key in [
            "ydisk",
            "ybul",
            "log10_gdagger",
            "log10_g0_cosmo",
            "log10_g0",
            "log10_kappa",
            "bridge_exponent",
            "inertia_amp",
            "inertia_index",
            "screen_index",
        ]:
            if key in row and row[key] != "":
                try:
                    pieces.append(f"{key}={float(row[key]):.5g}")
                except (TypeError, ValueError):
                    pieces.append(f"{key}={row[key]}")
        lines.append(
            f"| {float(row['H0_trial']):.3f} | {row['model']} | {float(row['chi2_dof']):.3f} | "
            f"{float(row['AIC']):.3f} | {float(row['BIC']):.3f} | {', '.join(pieces)} |"
        )

    lines.extend(
        [
            "",
            "## Nuisance Interpretation",
            "",
            "The galaxy-level nuisance pass applies one stellar M/L scale `eta` per galaxy after the global fit, with a Gaussian prior in `ln eta`.",
            "This is not yet a full hierarchical posterior, but it tells us whether the residuals are mostly stellar-population nuisance or functional-form mismatch.",
            "",
            "## Outputs",
            "",
            f"- Fit rows: `{OUTDIR / 'sparc_h0_nuisance_fit_summary.csv'}`",
            f"- Galaxy nuisance rows: `{OUTDIR / 'sparc_h0_nuisance_galaxy_ml.csv'}`",
            f"- H0 plot: `{OUTDIR / 'sparc_h0_nuisance_aic_vs_h0.png'}`",
            "",
            "## Next Research Options",
            "",
            "1. Turn the best PLAMB bridge branch into a proper likelihood with galaxy-level hierarchical priors.",
            "2. Test whether the bridge acceleration scale links to `c H0`, `c H0/2pi`, or another Noether/conservation scale.",
            "3. Repeat with external SPARC quality/distance subsets and compare residuals with surface brightness and HI mass.",
            "4. Add MaNGA/Galaxy Zoo angular-spin diagnostics separately for SU2 chirality/handedness.",
        ]
    )
    (OUTDIR / "sparc_h0_nuisance_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    meta = {
        "args": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
        "n_fit_rows": len(fit_rows),
        "n_galaxy_nuisance_rows": len(nuisance_rows),
        "best_global_aic": best_global_aic,
        "best_global_bic": best_global_bic,
    }
    (OUTDIR / "sparc_h0_nuisance_metadata.json").write_text(json.dumps(meta, indent=2, sort_keys=True), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Extended SPARC PLAMB H0/nuisance scan.")
    parser.add_argument("--sample", type=Path, default=DEFAULT_SAMPLE)
    parser.add_argument("--points", type=Path, default=DEFAULT_POINTS)
    parser.add_argument("--quality-max", type=int, default=2)
    parser.add_argument("--distance-method", choices=["all", "non_hubble", "hubble_flow"], default="all")
    parser.add_argument("--err-floor-kms", type=float, default=3.0)
    parser.add_argument("--h0-values", default="60,62.5,65,67.5,70,72.5,75")
    parser.add_argument("--h0-ref", type=float, default=SPARC_H0_REF)
    parser.add_argument(
        "--distance-scaling-mode",
        choices=["hubble_flow_only", "all", "none"],
        default="hubble_flow_only",
    )
    parser.add_argument("--maxiter", type=int, default=120)
    parser.add_argument("--seed", type=int, default=260713)
    parser.add_argument(
        "--models",
        default="all",
        help="Comma-separated model names to run, or 'all'.",
    )
    parser.add_argument("--skip-galaxy-ml", action="store_true")
    parser.add_argument("--galaxy-ml-sigma", type=float, default=0.35)
    parser.add_argument("--eta-min", type=float, default=0.25)
    parser.add_argument("--eta-max", type=float, default=4.0)
    parser.add_argument("--output-label", default="")
    args = parser.parse_args()

    outdir = OUTDIR
    if args.output_label:
        safe = "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in args.output_label)
        globals()["OUTDIR"] = outdir / safe

    base_data = load_dataset(args.sample, args.points, args.quality_max, args.distance_method, args.err_floor_kms)
    h0_values = parse_h0_values(args.h0_values)

    fit_rows: list[dict[str, Any]] = []
    nuisance_rows: list[dict[str, Any]] = []
    for h_index, h0 in enumerate(h0_values, start=1):
        data = rescale_for_h0(base_data, h0, args.h0_ref, args.distance_scaling_mode)
        specs = filter_model_specs(model_specs(h0), args.models)
        print(f"=== H0={h0:.3f} ({h_index}/{len(h0_values)}) ===")
        for s_index, spec in enumerate(specs, start=1):
            print(f"[{s_index}/{len(specs)}] {spec.name}")
            theta, chi2 = fit_model(spec, data, seed=args.seed + 100 * h_index + s_index, maxiter=args.maxiter)
            row = summarize_fit(
                h0,
                args.distance_scaling_mode,
                spec,
                theta,
                chi2,
                len(data.vobs),
                len(set(data.galaxy.tolist())),
            )
            fit_rows.append(row)
            print(f"    chi2/dof={row['chi2_dof']:.4g} AIC={row['AIC']:.4g} BIC={row['BIC']:.4g}")

            if not args.skip_galaxy_ml and spec.name != "BARYON_FREE_ML":
                chi2_data, chi2_prior, detail = galaxy_ml_nuisance(
                    spec,
                    theta,
                    data,
                    sigma_ln_eta=args.galaxy_ml_sigma,
                    eta_min=args.eta_min,
                    eta_max=args.eta_max,
                )
                nuisance_total = chi2_data + chi2_prior
                nuisance_row = summarize_fit(
                    h0,
                    args.distance_scaling_mode,
                    spec,
                    theta,
                    nuisance_total,
                    len(data.vobs),
                    len(set(data.galaxy.tolist())),
                    row_type="galaxy_ml_nuisance",
                    chi2_data=chi2_data,
                    chi2_prior=chi2_prior,
                    extra_k=len(set(data.galaxy.tolist())),
                )
                fit_rows.append(nuisance_row)
                for detail_row in detail:
                    detail_row["H0_trial"] = h0
                    detail_row["distance_scaling_mode"] = args.distance_scaling_mode
                nuisance_rows.extend(detail)
                print(
                    f"    +galaxy M/L: chi2_data={chi2_data:.4g} prior={chi2_prior:.4g} "
                    f"AIC_penalized={nuisance_row['AIC']:.4g}"
                )

    fieldnames = [
        "row_type",
        "H0_trial",
        "distance_scaling_mode",
        "model",
        "N_points",
        "N_galaxies",
        "k",
        "chi2",
        "chi2_data",
        "chi2_prior",
        "dof",
        "chi2_dof",
        "AIC",
        "BIC",
        "ydisk",
        "ybul",
            "log10_gdagger",
            "log10_g0_cosmo",
            "log10_g0",
        "log10_kappa",
        "bridge_exponent",
        "inertia_amp",
        "inertia_index",
        "screen_index",
        "note",
    ]
    OUTDIR.mkdir(parents=True, exist_ok=True)
    write_csv(OUTDIR / "sparc_h0_nuisance_fit_summary.csv", sorted(fit_rows, key=lambda row: float(row["AIC"])), fieldnames)
    write_csv(
        OUTDIR / "sparc_h0_nuisance_galaxy_ml.csv",
        nuisance_rows,
        [
            "H0_trial",
            "distance_scaling_mode",
            "model",
            "Galaxy",
            "eta_star_ml",
            "log_eta",
            "chi2_data",
            "chi2_prior",
            "N_points",
            "Q",
            "f_D",
        ],
    )
    make_h0_plot(fit_rows, OUTDIR / "sparc_h0_nuisance_aic_vs_h0.png")
    write_report(fit_rows, nuisance_rows, args, len(base_data.vobs), len(set(base_data.galaxy.tolist())))

    best_global = min([row for row in fit_rows if row["row_type"] == "global"], key=lambda row: float(row["AIC"]))
    print(f"Saved fit summary: {OUTDIR / 'sparc_h0_nuisance_fit_summary.csv'}")
    print(f"Saved report: {OUTDIR / 'sparc_h0_nuisance_report.md'}")
    print(
        f"Best global AIC: H0={best_global['H0_trial']:.3f} {best_global['model']} "
        f"AIC={best_global['AIC']:.4g} chi2/dof={best_global['chi2_dof']:.4g}"
    )


if __name__ == "__main__":
    main()
