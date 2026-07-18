#!/usr/bin/env python3
"""Derive and test a late-time PLAMB clock-law adapter against DESI DR2 BAO."""

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
from scipy.special import logsumexp
from scipy.stats import chi2 as chi2_distribution


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[3]
SHARED_CODE = REPO_ROOT / "github_export" / "code" / "shared"
if str(SHARED_CODE) not in sys.path:
    sys.path.insert(0, str(SHARED_CODE))

from plamb_clock_distance import clock_path_integral, clock_path_integrand  # noqa: E402


DATE_TAG = "2026-07-19"
GAMMA_FIXED = 2.3
P_SN_MEAN = 0.8004694772058842
P_SN_SIGMA = 0.02364051974806188
P_BOUNDS = (-0.5, 2.5)
OMEGA_BOUNDS = (0.05, 0.60)
DEFAULT_MEAN = (
    REPO_ROOT
    / "external_datasets"
    / "current_cosmology_sources"
    / "DESI_DR2_BAO"
    / "desi_gaussian_bao_ALL_GCcomb_mean.txt"
)
DEFAULT_COV = DEFAULT_MEAN.with_name("desi_gaussian_bao_ALL_GCcomb_cov.txt")
DEFAULT_FETCH_MANIFEST = (
    REPO_ROOT
    / "external_datasets"
    / "current_cosmology_sources"
    / "dataset_fetch_manifest.json"
)
DEFAULT_OUTDIR = (
    REPO_ROOT
    / "github_export"
    / "results"
    / DATE_TAG
    / "plamb_gamma2p3_desi_dr2_bao"
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_data(mean_path: Path, cov_path: Path) -> tuple[pd.DataFrame, np.ndarray]:
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
        raise ValueError(
            f"Covariance shape {covariance.shape} does not match {len(data)} BAO rows"
        )
    return data, covariance


def lcdm_dimensionless(z: np.ndarray, omega_m: float) -> tuple[np.ndarray, np.ndarray]:
    grid = np.linspace(0.0, float(np.max(z)), 20001)
    e_grid = np.sqrt(float(omega_m) * (1.0 + grid) ** 3 + 1.0 - float(omega_m))
    dz = grid[1] - grid[0]
    integral = np.empty_like(grid)
    integral[0] = 0.0
    integral[1:] = np.cumsum(0.5 * (1.0 / e_grid[1:] + 1.0 / e_grid[:-1]) * dz)
    dm = np.interp(z, grid, integral)
    dh = 1.0 / np.sqrt(float(omega_m) * (1.0 + z) ** 3 + 1.0 - float(omega_m))
    return dm, dh


def path_dimensionless(
    z: np.ndarray, p: float, branch: str
) -> tuple[np.ndarray, np.ndarray]:
    integral = np.asarray(clock_path_integral(z, GAMMA_FIXED, float(p)), dtype=float)
    kernel = np.asarray(clock_path_integrand(z, GAMMA_FIXED, float(p)), dtype=float)
    if branch == "stationary_path":
        return integral, kernel
    if branch == "metric_duality_control":
        one_plus_z = 1.0 + z
        dm = integral / one_plus_z
        dh = kernel / one_plus_z - integral / one_plus_z**2
        return dm, dh
    raise ValueError(branch)


def observable_vector(
    data: pd.DataFrame, dm: np.ndarray, dh: np.ndarray
) -> np.ndarray:
    z = data["z"].to_numpy(dtype=float)
    dv = np.cbrt(z * dm * dm * dh)
    values: list[float] = []
    for i, kind in enumerate(data["kind"].astype(str)):
        if kind.startswith("DM"):
            values.append(float(dm[i]))
        elif kind.startswith("DH"):
            values.append(float(dh[i]))
        elif kind.startswith("DV"):
            values.append(float(dv[i]))
        else:
            raise ValueError(f"Unknown BAO kind: {kind}")
    return np.asarray(values)


def profile_scale(
    shape: np.ndarray, observed: np.ndarray, precision: np.ndarray
) -> tuple[float, float, np.ndarray]:
    denominator = float(shape @ precision @ shape)
    if denominator <= 0.0:
        raise ValueError("BAO scale information is not positive")
    q = float(shape @ precision @ observed / denominator)
    residual = q * shape - observed
    chi2 = float(residual @ precision @ residual)
    return q, chi2, residual


def fit_path(
    data: pd.DataFrame,
    observed: np.ndarray,
    precision: np.ndarray,
    branch: str,
    p_mode: str,
) -> dict[str, object]:
    def at_p(p: float) -> tuple[float, float, np.ndarray, np.ndarray]:
        dm, dh = path_dimensionless(data["z"].to_numpy(dtype=float), p, branch)
        if np.any(dm <= 0.0) or np.any(dh <= 0.0):
            return float("nan"), float("inf"), np.full(len(data), np.nan), np.full(len(data), np.nan)
        shape = observable_vector(data, dm, dh)
        q, chi2, residual = profile_scale(shape, observed, precision)
        return q, chi2, residual, shape

    if p_mode == "fixed_sn":
        p = P_SN_MEAN
        success = True
        prior_chi2 = 0.0
    else:
        def objective(value: float) -> float:
            chi2 = at_p(float(value))[1]
            if p_mode == "sn_gaussian_prior":
                chi2 += ((float(value) - P_SN_MEAN) / P_SN_SIGMA) ** 2
            return chi2

        result = minimize_scalar(
            objective,
            bounds=P_BOUNDS,
            method="bounded",
            options={"xatol": 1e-11, "maxiter": 800},
        )
        candidates = [
            (float(result.fun), float(result.x)),
            (float(objective(P_BOUNDS[0])), P_BOUNDS[0]),
            (float(objective(P_BOUNDS[1])), P_BOUNDS[1]),
        ]
        _objective, p = min(candidates)
        success = bool(result.success or p in P_BOUNDS)
        prior_chi2 = (
            ((p - P_SN_MEAN) / P_SN_SIGMA) ** 2
            if p_mode == "sn_gaussian_prior"
            else 0.0
        )
    q, chi2, residual, shape = at_p(p)
    n = len(data)
    k = 1 if p_mode == "fixed_sn" else 2
    return {
        "model": f"PLAMB_{branch}_{p_mode}",
        "family": "PLAMB",
        "branch": branch,
        "p_mode": p_mode,
        "gamma_c": GAMMA_FIXED,
        "p": p,
        "p_shift_from_SN_sigma": (p - P_SN_MEAN) / P_SN_SIGMA,
        "q_c0_over_H0rd": q,
        "chi2_data": chi2,
        "chi2_SN_prior": prior_chi2,
        "posterior_objective": chi2 + prior_chi2,
        "N": n,
        "N_fitted_BAO_parameters": k,
        "BIC_data_only": chi2 + k * math.log(n),
        "goodness_upper_tail": float(chi2_distribution.sf(chi2, max(n - k, 1))),
        "optimisation_success": success,
        "at_p_bound": bool(min(p - P_BOUNDS[0], P_BOUNDS[1] - p) <= 0.01 * (P_BOUNDS[1] - P_BOUNDS[0])),
        "prediction": q * shape,
        "residual": residual,
    }


def fit_lcdm(
    data: pd.DataFrame, observed: np.ndarray, precision: np.ndarray
) -> dict[str, object]:
    z = data["z"].to_numpy(dtype=float)

    def at_omega(omega_m: float) -> tuple[float, float, np.ndarray, np.ndarray]:
        dm, dh = lcdm_dimensionless(z, omega_m)
        shape = observable_vector(data, dm, dh)
        q, chi2, residual = profile_scale(shape, observed, precision)
        return q, chi2, residual, shape

    result = minimize_scalar(
        lambda value: at_omega(float(value))[1],
        bounds=OMEGA_BOUNDS,
        method="bounded",
        options={"xatol": 1e-11, "maxiter": 800},
    )
    candidates = [
        (float(result.fun), float(result.x)),
        (float(at_omega(OMEGA_BOUNDS[0])[1]), OMEGA_BOUNDS[0]),
        (float(at_omega(OMEGA_BOUNDS[1])[1]), OMEGA_BOUNDS[1]),
    ]
    _objective, omega_m = min(candidates)
    q, chi2, residual, shape = at_omega(omega_m)
    n = len(data)
    return {
        "model": "LCDM_flat_Omega_m_free",
        "family": "LCDM",
        "branch": "flat_metric",
        "p_mode": "not_applicable",
        "gamma_c": np.nan,
        "p": np.nan,
        "p_shift_from_SN_sigma": np.nan,
        "Omega_m": omega_m,
        "q_c0_over_H0rd": q,
        "chi2_data": chi2,
        "chi2_SN_prior": 0.0,
        "posterior_objective": chi2,
        "N": n,
        "N_fitted_BAO_parameters": 2,
        "BIC_data_only": chi2 + 2 * math.log(n),
        "goodness_upper_tail": float(chi2_distribution.sf(chi2, max(n - 2, 1))),
        "optimisation_success": bool(result.success or omega_m in OMEGA_BOUNDS),
        "at_p_bound": False,
        "prediction": q * shape,
        "residual": residual,
    }


def posterior_grid(
    data: pd.DataFrame,
    observed: np.ndarray,
    precision: np.ndarray,
    branch: str,
    points: int = 3001,
) -> tuple[pd.DataFrame, dict[str, float]]:
    grid = np.linspace(P_BOUNDS[0], P_BOUNDS[1], points)
    chi2_values = np.empty(points)
    q_values = np.empty(points)
    for i, p in enumerate(grid):
        dm, dh = path_dimensionless(data["z"].to_numpy(dtype=float), p, branch)
        if np.any(dm <= 0.0) or np.any(dh <= 0.0):
            chi2_values[i] = np.inf
            q_values[i] = np.nan
            continue
        q, chi2, _ = profile_scale(observable_vector(data, dm, dh), observed, precision)
        q_values[i] = q
        chi2_values[i] = chi2
    prior_chi2 = ((grid - P_SN_MEAN) / P_SN_SIGMA) ** 2
    log_weights = -0.5 * (chi2_values + prior_chi2)
    log_weights -= logsumexp(log_weights)
    weights = np.exp(log_weights)
    mean = float(np.sum(grid * weights))
    variance = float(np.sum((grid - mean) ** 2 * weights))
    cumulative = np.cumsum(weights)
    summary = {
        "posterior_mean_p": mean,
        "posterior_standard_deviation_p": math.sqrt(max(variance, 0.0)),
        "posterior_q2p5_p": float(np.interp(0.025, cumulative, grid)),
        "posterior_q97p5_p": float(np.interp(0.975, cumulative, grid)),
    }
    return (
        pd.DataFrame(
            {
                "branch": branch,
                "p": grid,
                "q_c0_over_H0rd": q_values,
                "chi2_data": chi2_values,
                "chi2_SN_prior": prior_chi2,
                "posterior_weight": weights,
            }
        ),
        summary,
    )


def numerical_validation(data: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, bool]]:
    z_test = np.unique(np.r_[np.linspace(1e-5, 2.5, 1000), data["z"].to_numpy(dtype=float)])
    rows: list[dict[str, object]] = []
    gates: dict[str, bool] = {}
    step = 1e-6
    for branch in ("stationary_path", "metric_duality_control"):
        dm, dh = path_dimensionless(z_test, P_SN_MEAN, branch)
        dm_plus = path_dimensionless(z_test + step, P_SN_MEAN, branch)[0]
        dm_minus = path_dimensionless(np.maximum(z_test - step, 0.0), P_SN_MEAN, branch)[0]
        denominator = (z_test + step) - np.maximum(z_test - step, 0.0)
        derivative = (dm_plus - dm_minus) / denominator
        maximum_error = float(np.max(np.abs(derivative - dh)))
        positive = bool(np.all(dm > 0.0) and np.all(dh > 0.0))
        rows.append(
            {
                "branch": branch,
                "maximum_abs_dDMdz_minus_DH": maximum_error,
                "minimum_DM_dimensionless": float(np.min(dm)),
                "minimum_DH_dimensionless": float(np.min(dh)),
                "positive_to_z2p5": positive,
            }
        )
        gates[f"{branch}_derivative_identity"] = maximum_error <= 2e-7
        gates[f"{branch}_positive_to_z2p5"] = positive
    return pd.DataFrame(rows), gates


def per_redshift_contributions(
    data: pd.DataFrame,
    covariance: np.ndarray,
    fits: list[dict[str, object]],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    z_values = data["z"].to_numpy(dtype=float)
    observed = data["value"].to_numpy(dtype=float)
    for fit in fits:
        prediction = np.asarray(fit["prediction"])
        for z in np.unique(z_values):
            idx = np.flatnonzero(z_values == z)
            precision = np.linalg.inv(covariance[np.ix_(idx, idx)])
            residual = prediction[idx] - observed[idx]
            rows.append(
                {
                    "model": fit["model"],
                    "z": z,
                    "kinds": "+".join(data.iloc[idx]["kind"].astype(str)),
                    "N": len(idx),
                    "chi2_contribution": float(residual @ precision @ residual),
                }
            )
    return pd.DataFrame(rows)


def make_prediction_table(data: pd.DataFrame, fits: list[dict[str, object]]) -> pd.DataFrame:
    result = data.copy()
    result["sigma"] = np.nan
    for fit in fits:
        key = str(fit["model"])
        result[f"prediction__{key}"] = np.asarray(fit["prediction"])
        result[f"residual__{key}"] = np.asarray(fit["residual"])
    return result


def make_plot(path: Path, data: pd.DataFrame, covariance: np.ndarray, fits: list[dict[str, object]]) -> None:
    sigma = np.sqrt(np.diag(covariance))
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.7), constrained_layout=True)
    colours = {
        "PLAMB_stationary_path_fixed_sn": "#c43d3d",
        "PLAMB_metric_duality_control_fixed_sn": "#2b6cb0",
        "LCDM_flat_Omega_m_free": "#1f7a4d",
    }
    for ax, prefix, title in zip(axes, ("DM", "DH"), ("Transverse BAO", "Radial BAO")):
        mask = data["kind"].astype(str).str.startswith(prefix).to_numpy()
        ax.errorbar(
            data.loc[mask, "z"],
            data.loc[mask, "value"],
            yerr=sigma[mask],
            fmt="o",
            color="black",
            capsize=3,
            label="DESI DR2",
        )
        for fit in fits:
            if fit["model"] not in colours:
                continue
            order = np.argsort(data.loc[mask, "z"].to_numpy())
            ax.plot(
                data.loc[mask, "z"].to_numpy()[order],
                np.asarray(fit["prediction"])[mask][order],
                marker="s",
                linewidth=1.5,
                color=colours[fit["model"]],
                label=fit["model"].replace("PLAMB_", "").replace("_fixed_sn", "").replace("_", " "),
            )
        ax.set_title(title)
        ax.set_xlabel("Redshift z")
        ax.set_ylabel(f"{prefix} / r_d")
        ax.grid(alpha=0.25)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=4, fontsize=8)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def serialisable_fit(fit: dict[str, object]) -> dict[str, object]:
    return {
        key: value
        for key, value in fit.items()
        if key not in {"prediction", "residual"}
    }


def write_report(
    path: Path,
    fit_table: pd.DataFrame,
    validation: pd.DataFrame,
    contributions: pd.DataFrame,
    posterior_summaries: dict[str, dict[str, float]],
    theory_gates: dict[str, bool],
    statistical_gates: dict[str, bool],
    source_url: str,
) -> None:
    lcdm = fit_table[fit_table["model"] == "LCDM_flat_Omega_m_free"].iloc[0]
    stationary = fit_table[fit_table["model"] == "PLAMB_stationary_path_fixed_sn"].iloc[0]
    duality = fit_table[fit_table["model"] == "PLAMB_metric_duality_control_fixed_sn"].iloc[0]
    worst = contributions[
        contributions["model"].isin(
            ["PLAMB_stationary_path_fixed_sn", "PLAMB_metric_duality_control_fixed_sn"]
        )
    ].sort_values("chi2_contribution", ascending=False).head(10)
    lines = [
        "# PLAMB gamma = 2.3 observable adapter and DESI DR2 BAO audit",
        "",
        f"Generated: {datetime.now():%Y-%m-%d %H:%M:%S}",
        "",
        "## Decision",
        "",
        "**THE PHYSICAL BAO-ADAPTER GATE FAILS. DO NOT START AN OVERNIGHT SAMPLER.**",
        "",
        "A mathematically self-consistent late-time shape adapter can be written, but no tested "
        "branch simultaneously preserves the registered alpha=0 supernova distance, standard "
        "distance duality, the direct clock-law radial distance and an independently predicted "
        "sound horizon. The DESI diagnostic is nevertheless evaluated because its scale can be "
        "profiled analytically and it provides a useful falsification readout.",
        "",
        "## Registered clock law",
        "",
        "```text",
        "c(z) / c0                    = 1 + gamma_c z",
        "|dz/dT| / H0                = (1 + z)^p",
        "gamma_c                     = 2.3",
        "p_SN                        = 0.800469 +/- 0.023641",
        "I(z)                        = integral[0,z] (1+gamma_c u)/(1+u)^p du",
        "K(z)                        = dI/dz = (1+gamma_c z)/(1+z)^p",
        "q                           = c0/(H0 r_d)",
        "```",
        "",
        "The stationary-ruler branch uses `D_M/r_d=q I` and `D_H/r_d=q K`. This is the direct "
        "clock-path mapping, but with alpha=0 it implies `D_L=D_M`, which is incompatible with "
        "standard Etherington duality `D_L=(1+z)D_M`.",
        "",
        "The metric-duality control instead uses `D_M/r_d=q I/(1+z)` and differentiates this "
        "quantity for `D_H`. It preserves standard duality and radial/transverse derivative "
        "consistency, but its radial distance is no longer the direct clock kernel `q K`.",
        "",
        "The early-Universe sound horizon is not calculated. A common late-time scale `q` is "
        "profiled, so this is a BAO-shape test and does not extrapolate `c/c0=1+2.3z` to the drag epoch.",
        "",
        "## Numerical validation",
        "",
        validation.to_markdown(index=False),
        "",
        "## DESI DR2 likelihood",
        "",
        fit_table[
            [
                "model",
                "p",
                "p_shift_from_SN_sigma",
                "q_c0_over_H0rd",
                "chi2_data",
                "chi2_SN_prior",
                "BIC_data_only",
                "goodness_upper_tail",
                "at_p_bound",
            ]
        ].to_markdown(index=False),
        "",
        f"For the external fixed-SN prediction, the stationary branch gives chi2={stationary['chi2_data']:.3f} "
        f"and delta BIC={stationary['BIC_data_only'] - lcdm['BIC_data_only']:.3f} relative to flat Lambda-CDM.",
        f"The metric-duality control gives chi2={duality['chi2_data']:.3f} and "
        f"delta BIC={duality['BIC_data_only'] - lcdm['BIC_data_only']:.3f}.",
        "",
        "Allowing `p` to be selected by DESI is not an independent validation. It is reported only "
        "to diagnose whether the failure is a fixed-parameter mismatch; boundary solutions are gating failures.",
        "",
        "## SN-prior posterior updates",
        "",
    ]
    for branch, summary in posterior_summaries.items():
        lines.append(
            f"- `{branch}`: p = {summary['posterior_mean_p']:.6f} +/- "
            f"{summary['posterior_standard_deviation_p']:.6f}, 95% "
            f"[{summary['posterior_q2p5_p']:.6f}, {summary['posterior_q97p5_p']:.6f}]."
        )
    lines.extend(
        [
            "",
            "## Largest redshift-block contributions",
            "",
            worst.to_markdown(index=False),
            "",
            "## Theoretical gates",
            "",
        ]
    )
    for name, passed in theory_gates.items():
        lines.append(f"- `{name}`: **{'PASS' if passed else 'FAIL'}**")
    lines.extend(["", "## Statistical gates", ""])
    for name, passed in statistical_gates.items():
        lines.append(f"- `{name}`: **{'PASS' if passed else 'FAIL'}**")
    lines.extend(
        [
            "",
            "## Claim boundary",
            "",
            "The gamma value is outcome-aware and the p prior comes from the same supernova family "
            "being tested. DESI is an independent dataset, but the adapter lacks a complete physical "
            "ruler and reciprocity derivation. These results therefore reject the present adapter "
            "branches; they do not test every possible PLAMB completion.",
            "",
            "## Data source",
            "",
            f"- [DESI DR2 Gaussian BAO data product]({source_url})",
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
    data, covariance = load_data(args.mean, args.cov)
    observed = data["value"].to_numpy(dtype=float)
    eigenvalues = np.linalg.eigvalsh(covariance)
    if np.min(eigenvalues) <= 0.0:
        raise ValueError("DESI covariance is not positive definite")
    precision = np.linalg.inv(covariance)
    validation, numerical_gates = numerical_validation(data)

    protocol = {
        "analysis_id": "plamb_gamma2p3_desi_dr2_bao_2026-07-19",
        "status": "post_hoc_adapter_gate_and_external_diagnostic",
        "gamma_c_fixed": GAMMA_FIXED,
        "p_SN_prior": {"mean": P_SN_MEAN, "sigma": P_SN_SIGMA},
        "p_bounds": P_BOUNDS,
        "branches": {
            "stationary_path": "DM=qI; DH=q dI/dz",
            "metric_duality_control": "DM=qI/(1+z); DH=dDM/dz",
        },
        "scale": "q=c0/(H0 r_d), analytically profiled; no early-time r_d prediction",
        "overnight_rule": "launch only if one branch passes every physical theory gate",
        "data": {
            "mean_path": str(args.mean),
            "mean_sha256": sha256_file(args.mean),
            "cov_path": str(args.cov),
            "cov_sha256": sha256_file(args.cov),
        },
    }
    protocol_path = args.outdir / f"plamb_gamma2p3_desi_dr2_protocol_{DATE_TAG}.json"
    protocol_path.write_text(json.dumps(protocol, indent=2), encoding="utf-8")

    fits: list[dict[str, object]] = []
    for branch in ("stationary_path", "metric_duality_control"):
        for p_mode in ("fixed_sn", "sn_gaussian_prior", "free_bao"):
            fits.append(fit_path(data, observed, precision, branch, p_mode))
    fits.append(fit_lcdm(data, observed, precision))
    fit_table = pd.DataFrame([serialisable_fit(fit) for fit in fits])
    lcdm = next(fit for fit in fits if fit["model"] == "LCDM_flat_Omega_m_free")
    fit_table["delta_BIC_vs_LCDM"] = fit_table["BIC_data_only"] - float(lcdm["BIC_data_only"])

    grids: list[pd.DataFrame] = []
    posterior_summaries: dict[str, dict[str, float]] = {}
    for branch in ("stationary_path", "metric_duality_control"):
        grid, summary = posterior_grid(data, observed, precision, branch)
        grids.append(grid)
        posterior_summaries[branch] = summary
    posterior_table = pd.concat(grids, ignore_index=True)

    validation["covariance_minimum_eigenvalue"] = float(np.min(eigenvalues))
    validation["covariance_condition_number"] = float(np.max(eigenvalues) / np.min(eigenvalues))
    theory_gates = {
        **numerical_gates,
        "single_common_BAO_scale": True,
        "no_early_time_gamma_extrapolation": True,
        "stationary_branch_standard_distance_duality": False,
        "duality_branch_direct_clock_radial_mapping": False,
        "independent_sound_horizon_prediction": False,
        "one_branch_satisfies_all_physical_requirements": False,
    }
    fixed_stationary = next(
        fit for fit in fits if fit["model"] == "PLAMB_stationary_path_fixed_sn"
    )
    fixed_duality = next(
        fit for fit in fits if fit["model"] == "PLAMB_metric_duality_control_fixed_sn"
    )
    free_stationary = next(
        fit for fit in fits if fit["model"] == "PLAMB_stationary_path_free_bao"
    )
    statistical_gates = {
        "fixed_stationary_goodness_p_at_least_0p05": bool(fixed_stationary["goodness_upper_tail"] >= 0.05),
        "fixed_duality_goodness_p_at_least_0p05": bool(fixed_duality["goodness_upper_tail"] >= 0.05),
        "fixed_stationary_delta_BIC_vs_LCDM_at_most_2": bool(
            fixed_stationary["BIC_data_only"] - lcdm["BIC_data_only"] <= 2.0
        ),
        "fixed_duality_delta_BIC_vs_LCDM_at_most_2": bool(
            fixed_duality["BIC_data_only"] - lcdm["BIC_data_only"] <= 2.0
        ),
        "free_stationary_p_not_at_bound": bool(not free_stationary["at_p_bound"]),
        "covariance_positive_definite": True,
    }
    theory_pass = bool(all(theory_gates.values()))

    contributions = per_redshift_contributions(data, covariance, fits)
    predictions = make_prediction_table(data, fits)
    predictions["sigma"] = np.sqrt(np.diag(covariance))
    plot_fits = [fixed_stationary, fixed_duality, lcdm]
    plot_path = args.outdir / f"plamb_gamma2p3_desi_dr2_readout_{DATE_TAG}.png"
    make_plot(plot_path, data, covariance, plot_fits)

    fits_path = args.outdir / f"plamb_gamma2p3_desi_dr2_fits_{DATE_TAG}.csv"
    posterior_path = args.outdir / f"plamb_gamma2p3_desi_dr2_posterior_grid_{DATE_TAG}.csv.gz"
    contribution_path = args.outdir / f"plamb_gamma2p3_desi_dr2_block_chi2_{DATE_TAG}.csv"
    prediction_path = args.outdir / f"plamb_gamma2p3_desi_dr2_predictions_{DATE_TAG}.csv"
    validation_path = args.outdir / f"plamb_gamma2p3_desi_dr2_validation_{DATE_TAG}.csv"
    report_path = args.outdir / f"plamb_gamma2p3_desi_dr2_report_{DATE_TAG}.md"
    summary_path = args.outdir / f"plamb_gamma2p3_desi_dr2_summary_{DATE_TAG}.json"
    fit_table.to_csv(fits_path, index=False)
    posterior_table.to_csv(posterior_path, index=False, compression="gzip")
    contributions.to_csv(contribution_path, index=False)
    predictions.to_csv(prediction_path, index=False)
    validation.to_csv(validation_path, index=False)
    source_url = (
        "https://github.com/CobayaSampler/bao_data/tree/master/desi_bao_dr2"
    )
    write_report(
        report_path,
        fit_table,
        validation,
        contributions,
        posterior_summaries,
        theory_gates,
        statistical_gates,
        source_url,
    )
    summary = {
        "analysis_date": DATE_TAG,
        "theory_gate_pass": theory_pass,
        "overnight_sampler_launched": False,
        "overnight_sampler_reason": (
            "physical_theory_gate_failed" if not theory_pass else "not_launched_by_this_diagnostic"
        ),
        "theory_gates": theory_gates,
        "statistical_gates": statistical_gates,
        "posterior_summaries": posterior_summaries,
        "fits": {fit["model"]: serialisable_fit(fit) for fit in fits},
        "claim_promotion": False,
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    output_paths = [
        protocol_path,
        fits_path,
        posterior_path,
        contribution_path,
        prediction_path,
        validation_path,
        report_path,
        summary_path,
        plot_path,
    ]
    pd.DataFrame(
        [
            {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256_file(path)}
            for path in output_paths
        ]
    ).to_csv(
        args.outdir / f"plamb_gamma2p3_desi_dr2_manifest_{DATE_TAG}.csv",
        index=False,
    )
    print(f"Theory gate: {'PASS' if theory_pass else 'FAIL'}")
    for fit in fits:
        print(
            f"{fit['model']}: chi2={fit['chi2_data']:.6f} "
            f"q={fit['q_c0_over_H0rd']:.6f}",
            flush=True,
        )
    print(f"Saved report: {report_path}")


if __name__ == "__main__":
    main()
