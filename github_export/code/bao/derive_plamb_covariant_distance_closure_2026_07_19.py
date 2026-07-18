#!/usr/bin/env python3
"""Derive and test a covariant PLAMB clock, distance, and BAO-ruler closure."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import chi2 as chi2_distribution


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[3]
SHARED_CODE = REPO_ROOT / "github_export" / "code" / "shared"
if str(SHARED_CODE) not in sys.path:
    sys.path.insert(0, str(SHARED_CODE))

from plamb_clock_distance import clock_path_integral  # noqa: E402


DATE_TAG = "2026-07-19"
GAMMA_FIXED = 2.3
P_FIXED = 0.8004694772058842
P_SIGMA = 0.02364051974806188
CURVATURE_K = 0.0
RULER_POWER_BOUNDS = (-5.0, 5.0)
DEFAULT_MEAN = (
    REPO_ROOT
    / "external_datasets"
    / "current_cosmology_sources"
    / "DESI_DR2_BAO"
    / "desi_gaussian_bao_ALL_GCcomb_mean.txt"
)
DEFAULT_COV = DEFAULT_MEAN.with_name("desi_gaussian_bao_ALL_GCcomb_cov.txt")
DEFAULT_OUTDIR = (
    REPO_ROOT
    / "github_export"
    / "results"
    / DATE_TAG
    / "plamb_covariant_distance_closure"
)


CLOSURE_LABELS = {
    "static_fr": "Static FR: a=1",
    "fixed_atomic_clock": "Fixed atomic clock: Q=1",
    "standard_scale_factor": "Standard a=(1+z)^-1",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def power_integral(z: np.ndarray | float, exponent: float) -> np.ndarray:
    values = np.asarray(z, dtype=float)
    if abs(float(exponent)) <= 1e-10:
        return np.log1p(values)
    return np.expm1(float(exponent) * np.log1p(values)) / float(exponent)


def load_bao(mean_path: Path, cov_path: Path) -> tuple[pd.DataFrame, np.ndarray]:
    rows: list[dict[str, object]] = []
    for line in mean_path.read_text(encoding="utf-8").splitlines():
        text = line.strip()
        if not text or text.startswith("#"):
            continue
        z_raw, value_raw, kind_raw = text.split()
        rows.append(
            {
                "z": float(z_raw),
                "value": float(value_raw),
                "kind": kind_raw.replace("_over_rs", "_over_rd"),
            }
        )
    data = pd.DataFrame(rows)
    covariance = np.loadtxt(cov_path, dtype=float)
    covariance = 0.5 * (covariance + covariance.T)
    if covariance.shape != (len(data), len(data)):
        raise ValueError("DESI BAO covariance and mean-vector sizes differ")
    if np.min(np.linalg.eigvalsh(covariance)) <= 0.0:
        raise ValueError("DESI BAO covariance is not positive definite")
    return data, covariance


def background(z: np.ndarray, closure: str) -> dict[str, np.ndarray]:
    z = np.asarray(z, dtype=float)
    one_plus_z = 1.0 + z
    c_ratio = 1.0 + GAMMA_FIXED * z
    rate = one_plus_z**P_FIXED
    path = np.asarray(clock_path_integral(z, GAMMA_FIXED, P_FIXED), dtype=float)
    lookback = power_integral(z, 1.0 - P_FIXED)

    if closure == "static_fr":
        scale_factor = np.ones_like(z)
        atomic_clock = c_ratio / one_plus_z
        hubble = np.zeros_like(z)
        dlnq_dz = GAMMA_FIXED / c_ratio - 1.0 / one_plus_z
        chi = path
        dh = c_ratio / rate
    elif closure == "fixed_atomic_clock":
        scale_factor = c_ratio / one_plus_z
        atomic_clock = np.ones_like(z)
        hubble = rate * (1.0 / one_plus_z - GAMMA_FIXED / c_ratio)
        dlnq_dz = np.zeros_like(z)
        chi = power_integral(z, 2.0 - P_FIXED)
        dh = one_plus_z ** (1.0 - P_FIXED)
    elif closure == "standard_scale_factor":
        scale_factor = 1.0 / one_plus_z
        atomic_clock = c_ratio
        hubble = rate / one_plus_z
        dlnq_dz = GAMMA_FIXED / c_ratio
        chi = GAMMA_FIXED * power_integral(z, 3.0 - P_FIXED)
        chi += (1.0 - GAMMA_FIXED) * power_integral(z, 2.0 - P_FIXED)
        dh = c_ratio * one_plus_z ** (1.0 - P_FIXED)
    else:
        raise ValueError(closure)

    geometric_luminosity = one_plus_z**2 * scale_factor * chi
    xi = np.ones_like(z)
    nonzero = z > 0.0
    xi[nonzero] = path[nonzero] / geometric_luminosity[nonzero]
    luminosity_source_factor = 1.0 / np.square(xi)
    qdot_over_q_h0 = -rate * dlnq_dz
    cdot_over_c_h0 = -rate * GAMMA_FIXED / c_ratio
    rate_reconstructed = one_plus_z * (
        hubble + qdot_over_q_h0 - cdot_over_c_h0
    )
    redshift_reconstructed = c_ratio / (scale_factor * atomic_clock)

    return {
        "z": z,
        "C_c_over_c0": c_ratio,
        "R_abs_dz_dT_over_H0": rate,
        "H0_lookback_time": lookback,
        "a": scale_factor,
        "Q_atomic_frequency_ratio": atomic_clock,
        "H_over_H0": hubble,
        "D_path_H0_over_c0": path,
        "chi_H0_over_c0": chi,
        "DH_geometric_H0_over_c0": dh,
        "Xi_luminosity_transfer": xi,
        "J_photon_survival_times_luminosity_ratio": luminosity_source_factor,
        "redshift_reconstructed": redshift_reconstructed,
        "rate_reconstructed_over_H0": rate_reconstructed,
    }


def observable_vector(
    data: pd.DataFrame,
    dm: np.ndarray,
    dh: np.ndarray,
    ruler_power: float,
) -> np.ndarray:
    z = data["z"].to_numpy(dtype=float)
    ruler = (1.0 + z) ** float(ruler_power)
    dm_effective = dm / ruler
    dh_effective = dh / ruler
    dv_effective = np.cbrt(z * dm_effective * dm_effective * dh_effective)
    result: list[float] = []
    for i, kind in enumerate(data["kind"].astype(str)):
        if kind.startswith("DM"):
            result.append(float(dm_effective[i]))
        elif kind.startswith("DH"):
            result.append(float(dh_effective[i]))
        elif kind.startswith("DV"):
            result.append(float(dv_effective[i]))
        else:
            raise ValueError(f"Unknown BAO observable {kind}")
    return np.asarray(result)


def profile_scale(
    shape: np.ndarray, observed: np.ndarray, precision: np.ndarray
) -> tuple[float, float]:
    q = float(shape @ precision @ observed / (shape @ precision @ shape))
    residual = q * shape - observed
    return q, float(residual @ precision @ residual)


def fit_bao_closure(
    data: pd.DataFrame,
    covariance: np.ndarray,
    closure: str,
    ruler_mode: str,
) -> dict[str, object]:
    precision = np.linalg.inv(covariance)
    observed = data["value"].to_numpy(dtype=float)
    values = background(data["z"].to_numpy(dtype=float), closure)
    dm = values["chi_H0_over_c0"]
    dh = values["DH_geometric_H0_over_c0"]

    def score(ruler_power: float) -> tuple[float, float]:
        shape = observable_vector(data, dm, dh, ruler_power)
        return profile_scale(shape, observed, precision)

    if ruler_mode == "constant":
        ruler_power = 0.0
        success = True
    elif ruler_mode == "power_law":
        result = minimize_scalar(
            lambda value: score(float(value))[1],
            bounds=RULER_POWER_BOUNDS,
            method="bounded",
            options={"xatol": 1e-11, "maxiter": 800},
        )
        candidates = [
            (float(result.fun), float(result.x)),
            (float(score(RULER_POWER_BOUNDS[0])[1]), RULER_POWER_BOUNDS[0]),
            (float(score(RULER_POWER_BOUNDS[1])[1]), RULER_POWER_BOUNDS[1]),
        ]
        _minimum, ruler_power = min(candidates)
        success = bool(result.success or ruler_power in RULER_POWER_BOUNDS)
    else:
        raise ValueError(ruler_mode)

    q, chi2 = score(ruler_power)
    n = len(data)
    k = 1 if ruler_mode == "constant" else 2
    return {
        "model": f"{closure}_{ruler_mode}_ruler",
        "closure": closure,
        "ruler_mode": ruler_mode,
        "ruler_power_b": ruler_power,
        "B_z2p33": 3.33**ruler_power,
        "q_c0_over_H0rd": q,
        "chi2": chi2,
        "N": n,
        "N_parameters": k,
        "BIC": chi2 + k * math.log(n),
        "goodness_upper_tail": float(chi2_distribution.sf(chi2, max(n - k, 1))),
        "optimisation_success": success,
    }


def fit_lcdm(data: pd.DataFrame, covariance: np.ndarray) -> dict[str, object]:
    precision = np.linalg.inv(covariance)
    observed = data["value"].to_numpy(dtype=float)
    z = data["z"].to_numpy(dtype=float)

    def score(omega_m: float) -> tuple[float, float]:
        grid = np.linspace(0.0, float(np.max(z)), 20001)
        e_grid = np.sqrt(omega_m * (1.0 + grid) ** 3 + 1.0 - omega_m)
        step = grid[1] - grid[0]
        chi_grid = np.empty_like(grid)
        chi_grid[0] = 0.0
        chi_grid[1:] = np.cumsum(0.5 * (1.0 / e_grid[1:] + 1.0 / e_grid[:-1]) * step)
        dm = np.interp(z, grid, chi_grid)
        dh = 1.0 / np.sqrt(omega_m * (1.0 + z) ** 3 + 1.0 - omega_m)
        return profile_scale(observable_vector(data, dm, dh, 0.0), observed, precision)

    result = minimize_scalar(
        lambda value: score(float(value))[1],
        bounds=(0.05, 0.60),
        method="bounded",
        options={"xatol": 1e-11, "maxiter": 800},
    )
    omega_m = float(result.x)
    q, chi2 = score(omega_m)
    n = len(data)
    return {
        "model": "LCDM_flat",
        "closure": "LCDM_flat",
        "ruler_mode": "constant",
        "ruler_power_b": 0.0,
        "B_z2p33": 1.0,
        "q_c0_over_H0rd": q,
        "Omega_m": omega_m,
        "chi2": chi2,
        "N": n,
        "N_parameters": 2,
        "BIC": chi2 + 2 * math.log(n),
        "goodness_upper_tail": float(chi2_distribution.sf(chi2, n - 2)),
        "optimisation_success": bool(result.success),
    }


def local_ruler_audit(
    data: pd.DataFrame, covariance: np.ndarray, closure: str
) -> pd.DataFrame:
    z_all = data["z"].to_numpy(dtype=float)
    observed = data["value"].to_numpy(dtype=float)
    values = background(z_all, closure)
    shape = observable_vector(
        data,
        values["chi_H0_over_c0"],
        values["DH_geometric_H0_over_c0"],
        0.0,
    )
    local: list[dict[str, object]] = []
    for z in np.unique(z_all):
        idx = np.flatnonzero(z_all == z)
        local_precision = np.linalg.inv(covariance[np.ix_(idx, idx)])
        scale, chi2 = profile_scale(shape[idx], observed[idx], local_precision)
        row: dict[str, object] = {
            "closure": closure,
            "z": z,
            "N": len(idx),
            "local_q_over_B": scale,
            "irreducible_AP_chi2": chi2,
            "B_DM_over_B_DH": np.nan,
        }
        kinds = data.iloc[idx]["kind"].astype(str).to_numpy(dtype=str)
        dm_local = np.flatnonzero(np.char.startswith(kinds, "DM"))
        dh_local = np.flatnonzero(np.char.startswith(kinds, "DH"))
        if len(dm_local) == 1 and len(dh_local) == 1:
            i_dm = idx[dm_local[0]]
            i_dh = idx[dh_local[0]]
            b_dm_without_q = shape[i_dm] / observed[i_dm]
            b_dh_without_q = shape[i_dh] / observed[i_dh]
            row["B_DM_over_B_DH"] = b_dm_without_q / b_dh_without_q
        local.append(row)
    result = pd.DataFrame(local)
    anchor = float(result.loc[result["z"].idxmin(), "local_q_over_B"])
    result["B_relative_to_lowz_anchor"] = anchor / result["local_q_over_B"]
    return result


def validation_table() -> tuple[pd.DataFrame, dict[str, bool]]:
    z = np.linspace(1e-6, 2.5, 5001)
    rows: list[dict[str, object]] = []
    gates: dict[str, bool] = {}
    for closure in CLOSURE_LABELS:
        values = background(z, closure)
        dz = z[1] - z[0]
        derivative = np.gradient(values["chi_H0_over_c0"], dz, edge_order=2)
        redshift_error = float(
            np.max(np.abs(values["redshift_reconstructed"] - (1.0 + z)))
        )
        rate_error = float(
            np.max(
                np.abs(
                    values["rate_reconstructed_over_H0"]
                    - values["R_abs_dz_dT_over_H0"]
                )
            )
        )
        derivative_error = float(
            np.max(np.abs(derivative - values["DH_geometric_H0_over_c0"]))
        )
        luminosity_reconstructed = (
            (1.0 + z) ** 2
            * values["a"]
            * values["chi_H0_over_c0"]
            * values["Xi_luminosity_transfer"]
        )
        luminosity_error = float(
            np.max(np.abs(luminosity_reconstructed - values["D_path_H0_over_c0"]))
        )
        rows.append(
            {
                "closure": closure,
                "maximum_redshift_identity_error": redshift_error,
                "maximum_rate_identity_error": rate_error,
                "maximum_dchi_dz_minus_DH_error": derivative_error,
                "maximum_luminosity_identity_error": luminosity_error,
                "minimum_a": float(np.min(values["a"])),
                "maximum_a": float(np.max(values["a"])),
                "H_over_H0_at_z0": float(background(np.asarray([0.0]), closure)["H_over_H0"][0]),
            }
        )
        gates[f"{closure}_redshift_identity"] = redshift_error <= 1e-10
        gates[f"{closure}_rate_identity"] = rate_error <= 1e-10
        gates[f"{closure}_radial_derivative_identity"] = derivative_error <= 2e-6
        gates[f"{closure}_luminosity_identity"] = luminosity_error <= 1e-10
    return pd.DataFrame(rows), gates


def make_curve_table() -> pd.DataFrame:
    redshifts = np.asarray([0.0, 0.295, 0.5, 0.706, 0.934, 1.0, 1.321, 1.484, 2.0, 2.33])
    rows: list[pd.DataFrame] = []
    for closure in CLOSURE_LABELS:
        frame = pd.DataFrame(background(redshifts, closure))
        frame.insert(0, "closure", closure)
        rows.append(frame)
    return pd.concat(rows, ignore_index=True)


def observed_ap_ratios(data: pd.DataFrame, covariance: np.ndarray) -> pd.DataFrame:
    rows: list[dict[str, float]] = []
    z_all = data["z"].to_numpy(dtype=float)
    for z in np.unique(z_all):
        idx = np.flatnonzero(z_all == z)
        kinds = data.iloc[idx]["kind"].astype(str).to_numpy(dtype=str)
        dm_local = np.flatnonzero(np.char.startswith(kinds, "DM"))
        dh_local = np.flatnonzero(np.char.startswith(kinds, "DH"))
        if len(dm_local) != 1 or len(dh_local) != 1:
            continue
        i_dm = idx[dm_local[0]]
        i_dh = idx[dh_local[0]]
        dm = float(data.iloc[i_dm]["value"])
        dh = float(data.iloc[i_dh]["value"])
        ratio = dh / dm
        gradient = np.asarray([-dh / dm**2, 1.0 / dm])
        pair_covariance = covariance[np.ix_([i_dm, i_dh], [i_dm, i_dh])]
        sigma = math.sqrt(float(gradient @ pair_covariance @ gradient))
        rows.append({"z": z, "ratio": ratio, "sigma": sigma})
    return pd.DataFrame(rows)


def make_plot(
    path: Path,
    data: pd.DataFrame,
    covariance: np.ndarray,
) -> None:
    z = np.linspace(0.001, 2.33, 600)
    z_ratio = np.linspace(0.25, 2.33, 500)
    colours = {
        "static_fr": "#c43d3d",
        "fixed_atomic_clock": "#2b6cb0",
        "standard_scale_factor": "#1f7a4d",
    }
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 8.0), constrained_layout=True)
    for closure, label in CLOSURE_LABELS.items():
        values = background(z, closure)
        ratio_values = background(z_ratio, closure)
        colour = colours[closure]
        axes[0, 0].plot(z, values["Q_atomic_frequency_ratio"], color=colour, label=label)
        axes[0, 1].plot(z, values["a"], color=colour, label=label)
        axes[1, 0].plot(z, values["Xi_luminosity_transfer"], color=colour, label=label)
        axes[1, 1].plot(
            z_ratio,
            ratio_values["DH_geometric_H0_over_c0"]
            / ratio_values["chi_H0_over_c0"],
            color=colour,
            label=label,
        )
    ap = observed_ap_ratios(data, covariance)
    axes[1, 1].errorbar(
        ap["z"], ap["ratio"], yerr=ap["sigma"], fmt="o", color="black", capsize=3, label="DESI DR2"
    )
    axes[0, 0].set_title("Required atomic-clock evolution")
    axes[0, 0].set_ylabel("Q(z)")
    axes[0, 1].set_title("Scale-factor closure")
    axes[0, 1].set_ylabel("a(z)")
    axes[1, 0].set_title("Required luminosity transfer")
    axes[1, 0].set_ylabel("Xi(z)")
    axes[1, 1].set_title("Ruler-independent radial/transverse ratio")
    axes[1, 1].set_ylabel("D_H / D_M")
    axes[1, 1].set_xlim(0.25, 2.4)
    axes[1, 1].set_ylim(0.0, 5.0)
    for ax in axes.flat:
        ax.set_xlabel("Redshift z")
        ax.grid(alpha=0.25)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    handles_ap, labels_ap = axes[1, 1].get_legend_handles_labels()
    if "DESI DR2" in labels_ap:
        index = labels_ap.index("DESI DR2")
        handles.append(handles_ap[index])
        labels.append(labels_ap[index])
    fig.legend(handles, labels, loc="outside lower center", ncol=4, fontsize=8)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_report(
    path: Path,
    curves: pd.DataFrame,
    validation: pd.DataFrame,
    fits: pd.DataFrame,
    local_rulers: pd.DataFrame,
    gates: dict[str, bool],
) -> None:
    endpoint = curves[np.isclose(curves["z"], 2.33)].copy()
    lcdm = fits[fits["model"] == "LCDM_flat"].iloc[0]
    fit_display = fits.copy()
    fit_display["delta_BIC_vs_LCDM"] = fit_display["BIC"] - float(lcdm["BIC"])
    ap_summary = (
        local_rulers.groupby("closure", as_index=False)
        .agg(
            irreducible_AP_chi2=("irreducible_AP_chi2", "sum"),
            maximum_B_DM_over_B_DH=("B_DM_over_B_DH", "max"),
            minimum_B_DM_over_B_DH=("B_DM_over_B_DH", "min"),
            B_relative_at_z2p33=("B_relative_to_lowz_anchor", "last"),
        )
    )
    lines = [
        "# Covariant PLAMB clock-distance-ruler closure",
        "",
        f"Generated: {datetime.now():%Y-%m-%d %H:%M:%S}",
        "",
        "## Decision",
        "",
        "**KINEMATIC DERIVATION COMPLETE; PREDICTIVE PHYSICAL CLOSURE FAILS.**",
        "",
        "A single covariant optical and ruler relation can be derived, but the current PLAMB "
        "inputs do not close it without additional clock, luminosity and ruler functions. Under "
        "the intended static background, the required functions are fixed algebraically and are "
        "large. DESI's radial/transverse ratio cannot be repaired by any isotropic ruler evolution.",
        "",
        "This is a covariant kinematic completion, not an action-level FR field theory. A full "
        "theory must still derive the clock scalar dynamics, stress-energy and electromagnetic coupling.",
        "",
        "## Covariant setup",
        "",
        "Let `T(x)` be a scalar clock field with unit timelike flow `u_mu`. Matter measures proper "
        "time on the FLRW metric `g_mu_nu`, while geometric-optics photons follow the disformal metric",
        "",
        "```text",
        "g_gamma(mu,nu) = g(mu,nu) + [1-C(T)^2] u_mu u_nu",
        "ds_gamma^2    = -C(T)^2 c0^2 dT^2 + a(T)^2 dSigma_K^2.",
        "```",
        "",
        "For a geometric-optics wave covector `k_mu`, the frequency measured by a comoving "
        "matter observer is `omega=-u^mu k_mu`, and spatial homogeneity gives "
        "`omega proportional to C/a`. Comparison with the local atomic transition therefore "
        "introduces the additional ratio `Q(T)=nu_atom(T)/nu_atom(T0)`.",
        "",
        "Here `C=c_gamma/c0`. Let `Q(T)=nu_atom(T)/nu_atom(T0)` describe the evolution of the "
        "atomic transition used to define spectroscopic redshift. With present values `a0=C0=Q0=1`,",
        "",
        "$$",
        "\\begin{aligned}",
        "1+z                         &= \\frac{C(T)}{a(T)Q(T)}, \\\\",
        "R(z)\\equiv-\\frac{dz}{dT} &= (1+z)\\left[H+\\frac{\\dot Q}{Q}-\\frac{\\dot C}{C}\\right], \\\\",
        "T_0-T(z)                    &= \\int_0^z\\frac{du}{R(u)}.",
        "\\end{aligned}",
        "$$",
        "",
        "The registered phenomenology supplies only `C(z)=1+gamma_c z` and "
        "`R(z)=H0(1+z)^p`. It does not separately determine `a(T)` and `Q(T)`.",
        "",
        "## Optical and BAO distances",
        "",
        "Let `B(z)` be the evolution of the comoving BAO ruler, "
        "`r_BAO(z)=r_d B(z)`, and let `Xi(z)` collect photon survival and standardised intrinsic "
        "luminosity evolution. The flat-background equations are",
        "",
        "$$",
        "\\begin{aligned}",
        "\\chi(z)                   &= \\int_0^z\\frac{c_0 C(u)}{a(u)R(u)}\\,du, \\\\",
        "D_A(z)                      &= a(z)S_K[\\chi(z)], \\\\",
        "\\frac{D_M^{\\rm BAO}}{r_d} &= \\frac{S_K(\\chi)}{r_d B(z)}, \\\\",
        "\\frac{D_H^{\\rm BAO}}{r_d} &= \\frac{c_0 C(z)}{a(z)R(z)r_d B(z)}, \\\\",
        "D_L^{\\rm SN}(z)           &= (1+z)^2D_A(z)\\,\\Xi(z) \\\\",
        "                              &= (1+z)^2a(z)B(z)D_M^{\\rm BAO}(z)\\,\\Xi(z).",
        "\\end{aligned}",
        "$$",
        "",
        "This is the requested single relation connecting `T`, `D_L`, `D_M`, `D_H` and ruler "
        "evolution. For curvature, `S_K` is the usual transverse-curvature map. Differentiation gives",
        "",
        "$$",
        "\\begin{aligned}",
        "\\frac{dD_M^{\\rm BAO}}{dz}+\\frac{d\\ln B}{dz}D_M^{\\rm BAO}",
        "&= C_K(\\chi)D_H^{\\rm BAO}.",
        "\\end{aligned}",
        "$$",
        "",
        "Thus `D_H=dD_M/dz` is valid only for flat geometry and a non-evolving comoving ruler.",
        "",
        "## Relation to the fitted PLAMB path",
        "",
        "The supernova programme fitted",
        "",
        "$$",
        "\\begin{aligned}",
        "D_{\\rm path}(z) &= \\int_0^z\\frac{c_0C(u)}{R(u)}\\,du, \\\\",
        "D_L^{\\rm fit}   &= D_{\\rm path}.",
        "\\end{aligned}",
        "$$",
        "",
        "`D_path` equals the comoving optical distance `chi` only when `a=1`. Under that static "
        "closure with a fixed BAO ruler, consistency forces",
        "",
        "$$",
        "\\begin{aligned}",
        "Q(z)                            &= \\frac{1+2.3z}{1+z}, \\\\",
        "\\Xi(z)                         &= (1+z)^{-2}, \\\\",
        "P_\\gamma(z)\\frac{L_*(z)}{L_{*,0}} &= \\Xi(z)^{-2}=(1+z)^4.",
        "\\end{aligned}",
        "$$",
        "",
        "Consequently, alpha=0 is not photon-conserving Etherington propagation. It requires a "
        "specific composite luminosity or photon-transfer law. At fixed atomic standards (`Q=1`), "
        "the same gamma law instead requires `a=C/(1+z)` and gives `H(0)/H0=1-gamma_c=-1.3`, "
        "a contracting present background.",
        "",
        "## Endpoint requirements at z=2.33",
        "",
        endpoint[
            [
                "closure",
                "a",
                "Q_atomic_frequency_ratio",
                "H_over_H0",
                "Xi_luminosity_transfer",
                "J_photon_survival_times_luminosity_ratio",
            ]
        ].to_markdown(index=False),
        "",
        "## Numerical identity checks",
        "",
        validation.to_markdown(index=False),
        "",
        "## DESI ruler tests",
        "",
        fit_display[
            [
                "model",
                "ruler_power_b",
                "B_z2p33",
                "q_c0_over_H0rd",
                "chi2",
                "BIC",
                "delta_BIC_vs_LCDM",
                "goodness_upper_tail",
            ]
        ].to_markdown(index=False),
        "",
        "A power-law ruler `B(z)=(1+z)^b` cannot rescue any closure. More generally, allowing an "
        "independent isotropic ruler amplitude in every DESI redshift bin leaves the following "
        "irreducible Alcock-Paczynski radial/transverse scores:",
        "",
        ap_summary.to_markdown(index=False),
        "",
        "The ratio `B_DM/B_DH` would equal one if a single isotropic ruler could fit both "
        "directions. Its departure from one is independent of the absolute sound horizon and "
        "cannot be corrected by choosing another `r_d`.",
        "",
        "## Gates",
        "",
    ]
    for name, passed in gates.items():
        lines.append(f"- `{name}`: **{'PASS' if passed else 'FAIL'}**")
    lines.extend(
        [
            "",
            "## Research implication",
            "",
            "Further sampling is not justified for the present `gamma_c=2.3`, `p=0.800469`, "
            "alpha-zero branch. The next theory submission must specify, before seeing another "
            "distance dataset: (1) the action or field equations for `T`; (2) whether matter and "
            "photons share a metric; (3) the atomic-frequency law `Q(T)`; (4) photon conservation "
            "and intrinsic SN luminosity, fixing `Xi(T)`; and (5) the drag-epoch calculation and "
            "late-time transport of `B(T)`. These are physical predictions, not nuisance parameters "
            "that may all be fitted to the same SN and BAO data.",
            "",
            "## Primary references",
            "",
            "- [Moffat bimetric VSL construction](https://arxiv.org/abs/gr-qc/0202012)",
            "- [Disformal invariance of cosmological observables](https://arxiv.org/abs/2003.10633)",
            "- [Etherington distance-duality assumptions](https://arxiv.org/abs/1612.08784)",
            "- [Joint-constant requirements in a covariant VSL model](https://arxiv.org/abs/2011.09274)",
            "- [DESI DR2 publications and products](https://data.desi.lbl.gov/doc/papers/dr2/)",
            "- [Peter R. Lamb, Making Sense of Gravity, v8.10](https://www.fullyrelative.com/wp-content/uploads/2023/01/Making-Sense-of-Gravity-book-Vs8-10.pdf)",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mean", type=Path, default=DEFAULT_MEAN)
    parser.add_argument("--cov", type=Path, default=DEFAULT_COV)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    data, covariance = load_bao(args.mean, args.cov)
    protocol = {
        "analysis_id": "plamb_covariant_distance_closure_2026-07-19",
        "status": "post_hoc_kinematic_derivation_not_action_level_theory",
        "gamma_c": GAMMA_FIXED,
        "p": P_FIXED,
        "p_sigma": P_SIGMA,
        "curvature_K": CURVATURE_K,
        "closures": CLOSURE_LABELS,
        "ruler_model": "B(z)=(1+z)^b plus arbitrary-bin isotropic amplitude audit",
        "source_hashes": {
            "DESI_mean": sha256_file(args.mean),
            "DESI_covariance": sha256_file(args.cov),
        },
        "claim_boundary": "covariant optical closure and falsification audit; not an FR action or discovery test",
    }
    protocol_path = args.outdir / f"plamb_covariant_closure_protocol_{DATE_TAG}.json"
    protocol_path.write_text(json.dumps(protocol, indent=2), encoding="utf-8")

    curves = make_curve_table()
    validation, identity_gates = validation_table()
    fit_rows: list[dict[str, object]] = []
    for closure in CLOSURE_LABELS:
        fit_rows.append(fit_bao_closure(data, covariance, closure, "constant"))
        fit_rows.append(fit_bao_closure(data, covariance, closure, "power_law"))
    fit_rows.append(fit_lcdm(data, covariance))
    fits = pd.DataFrame(fit_rows)
    local_rulers = pd.concat(
        [local_ruler_audit(data, covariance, closure) for closure in CLOSURE_LABELS],
        ignore_index=True,
    )
    local_ap = local_rulers.groupby("closure")["irreducible_AP_chi2"].sum()
    power_fits = fits[fits["ruler_mode"] == "power_law"]
    static_endpoint = curves[
        (curves["closure"] == "static_fr") & np.isclose(curves["z"], 2.33)
    ].iloc[0]
    gates = {
        **identity_gates,
        "static_background_matches_FR_intent": True,
        "static_atomic_clock_evolution_is_derived_not_supplied": False,
        "static_luminosity_transfer_is_derived_not_supplied": False,
        "static_Xi_close_to_unity_at_z2p33": bool(
            abs(float(static_endpoint["Xi_luminosity_transfer"]) - 1.0) <= 0.05
        ),
        "fixed_atomic_clock_closure_is_noncontracting_today": bool(
            validation.loc[
                validation["closure"] == "fixed_atomic_clock", "H_over_H0_at_z0"
            ].iloc[0]
            > 0.0
        ),
        "power_law_ruler_any_closure_goodness_p_at_least_0p05": bool(
            (power_fits["goodness_upper_tail"] >= 0.05).any()
        ),
        "arbitrary_isotropic_ruler_any_closure_AP_p_at_least_0p05": bool(
            any(chi2_distribution.sf(value, 6) >= 0.05 for value in local_ap)
        ),
        "independent_drag_epoch_ruler_derivation": False,
        "action_level_clock_and_photon_dynamics": False,
    }

    curves_path = args.outdir / f"plamb_covariant_closure_curves_{DATE_TAG}.csv"
    validation_path = args.outdir / f"plamb_covariant_closure_validation_{DATE_TAG}.csv"
    fits_path = args.outdir / f"plamb_covariant_closure_bao_fits_{DATE_TAG}.csv"
    rulers_path = args.outdir / f"plamb_covariant_closure_local_rulers_{DATE_TAG}.csv"
    summary_path = args.outdir / f"plamb_covariant_closure_summary_{DATE_TAG}.json"
    report_path = args.outdir / f"plamb_covariant_closure_report_{DATE_TAG}.md"
    plot_path = args.outdir / f"plamb_covariant_closure_readout_{DATE_TAG}.png"
    curves.to_csv(curves_path, index=False)
    validation.to_csv(validation_path, index=False)
    fits.to_csv(fits_path, index=False)
    local_rulers.to_csv(rulers_path, index=False)
    make_plot(plot_path, data, covariance)
    write_report(report_path, curves, validation, fits, local_rulers, gates)
    summary = {
        "analysis_date": DATE_TAG,
        "kinematic_identity_gates_pass": bool(all(identity_gates.values())),
        "predictive_physical_closure_pass": bool(all(gates.values())),
        "gates": gates,
        "endpoint_z2p33": {
            row["closure"]: {
                "a": float(row["a"]),
                "Q": float(row["Q_atomic_frequency_ratio"]),
                "H_over_H0": float(row["H_over_H0"]),
                "Xi": float(row["Xi_luminosity_transfer"]),
                "J": float(row["J_photon_survival_times_luminosity_ratio"]),
            }
            for _, row in curves[np.isclose(curves["z"], 2.33)].iterrows()
        },
        "bao_fits": fit_rows,
        "irreducible_AP_chi2": {key: float(value) for key, value in local_ap.items()},
        "claim_promotion": False,
        "next_sampling_authorised": False,
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    output_paths = [
        protocol_path,
        curves_path,
        validation_path,
        fits_path,
        rulers_path,
        summary_path,
        report_path,
        plot_path,
    ]
    pd.DataFrame(
        [
            {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256_file(path)}
            for path in output_paths
        ]
    ).to_csv(
        args.outdir / f"plamb_covariant_closure_manifest_{DATE_TAG}.csv",
        index=False,
    )
    print(f"Kinematic identity gates: {'PASS' if all(identity_gates.values()) else 'FAIL'}")
    print(f"Predictive physical closure: {'PASS' if all(gates.values()) else 'FAIL'}")
    for _, row in fits.iterrows():
        print(f"{row['model']}: chi2={row['chi2']:.6f} BIC={row['BIC']:.6f}")
    print(f"Saved report: {report_path}")


if __name__ == "__main__":
    main()
