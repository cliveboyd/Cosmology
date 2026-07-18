#!/usr/bin/env python3
"""Analyse the outcome-aware fixed-gamma=2.3 PLAMB target."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from dataclasses import dataclass
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
RAWMU_CODE = SCRIPT_PATH.parent
SHARED_CODE = REPO_ROOT / "github_export" / "code" / "shared"
for code_path in (RAWMU_CODE, SHARED_CODE):
    if str(code_path) not in sys.path:
        sys.path.insert(0, str(code_path))

import run_plamb_nested_clock_fit_2026_07_18 as nested  # noqa: E402
from plamb_clock_distance import clock_path_integral  # noqa: E402


DATE_TAG = "2026-07-19"
GAMMA_FIXED = 2.3
ALPHA_FIXED = 0.0
TRAINING_LIMIT = 0.5
DEFAULT_OUTDIR = (
    REPO_ROOT / "github_export" / "results" / DATE_TAG / "rawmu_plamb_gamma2p3"
)
DEFAULT_PROTOCOL = DEFAULT_OUTDIR / f"plamb_gamma2p3_protocol_{DATE_TAG}.json"
BOUNDS = {
    "PATH_GAMMA2P3": (-0.5, 2.5),
    "LCDM": (0.05, 0.60),
}


@dataclass
class PredictivePartition:
    block: nested.base.Block
    low_indices: np.ndarray
    high_indices: np.ndarray
    conditional_projection: np.ndarray
    intercept_direction: np.ndarray
    intercept_variance: float
    predictive_cholesky: np.ndarray
    predictive_logdet: float
    low_precision: np.ndarray
    low_intercept_denominator: float
    cholesky_jitter: float

    @property
    def n_high(self) -> int:
        return int(len(self.high_indices))


def stable_cholesky(matrix: np.ndarray) -> tuple[np.ndarray, float]:
    matrix = 0.5 * (matrix + matrix.T)
    scale = max(float(np.median(np.diag(matrix))), 1e-16)
    for multiplier in (0.0, 1e-12, 1e-10, 1e-8):
        jitter = multiplier * scale
        try:
            return np.linalg.cholesky(matrix + jitter * np.eye(len(matrix))), jitter
        except np.linalg.LinAlgError:
            continue
    raise np.linalg.LinAlgError("Predictive covariance is not positive definite")


def make_low_blocks(blocks: list[nested.base.Block]) -> list[nested.base.Block]:
    result = []
    for block in blocks:
        keep = block.z < TRAINING_LIMIT
        if not np.any(keep) or np.all(keep):
            raise ValueError(f"{block.label} does not span the z=0.5 split")
        result.append(nested.base.subset_block(block, keep, "train_z_lt_0p5"))
    return result


def make_predictive_partitions(
    blocks: list[nested.base.Block],
) -> list[PredictivePartition]:
    result = []
    for block in blocks:
        low_indices = np.flatnonzero(block.z < TRAINING_LIMIT)
        high_unsorted = np.flatnonzero(block.z >= TRAINING_LIMIT)
        order = np.lexsort(
            (block.source_indices[high_unsorted], block.z[high_unsorted])
        )
        high_indices = high_unsorted[order]
        c_ll = block.covariance[np.ix_(low_indices, low_indices)]
        c_lh = block.covariance[np.ix_(low_indices, high_indices)]
        c_hh = block.covariance[np.ix_(high_indices, high_indices)]
        low_precision = nested.base.stable_inverse(c_ll)
        ones_low = np.ones(len(low_indices))
        denominator = float(ones_low @ low_precision @ ones_low)
        intercept_variance = 1.0 / denominator
        solved = np.linalg.solve(c_ll, c_lh)
        projection = solved.T
        conditional_covariance = c_hh - c_lh.T @ solved
        intercept_direction = (
            np.ones(len(high_indices)) - projection @ ones_low
        )
        predictive_covariance = conditional_covariance + intercept_variance * np.outer(
            intercept_direction, intercept_direction
        )
        cholesky, jitter = stable_cholesky(predictive_covariance)
        logdet = float(2.0 * np.sum(np.log(np.diag(cholesky))))
        result.append(
            PredictivePartition(
                block=block,
                low_indices=low_indices,
                high_indices=high_indices,
                conditional_projection=projection,
                intercept_direction=intercept_direction,
                intercept_variance=intercept_variance,
                predictive_cholesky=cholesky,
                predictive_logdet=logdet,
                low_precision=low_precision,
                low_intercept_denominator=denominator,
                cholesky_jitter=jitter,
            )
        )
    return result


def model_parameters(model_key: str, theta: float) -> dict[str, float]:
    if model_key == "PATH_GAMMA2P3":
        return {
            "gamma_c": GAMMA_FIXED,
            "p": float(theta),
            "alpha": ALPHA_FIXED,
        }
    if model_key == "LCDM":
        return {"Omega_m": float(theta)}
    raise ValueError(model_key)


def model_spec(model_key: str) -> nested.ModelSpec:
    return (
        nested.MODEL_SPECS["GENERAL_FREE"]
        if model_key == "PATH_GAMMA2P3"
        else nested.MODEL_SPECS["LCDM"]
    )


def score_blocks(
    blocks: list[nested.base.Block],
    model_key: str,
    theta: float,
    h0: float,
) -> nested.Score:
    if model_key == "PATH_GAMMA2P3":
        return nested.score_clock_profiled_alpha(
            blocks, GAMMA_FIXED, float(theta), ALPHA_FIXED, h0
        )
    if model_key == "LCDM":
        return nested.score_lcdm(blocks, float(theta), h0)
    raise ValueError(model_key)


def fit_scalar(
    blocks: list[nested.base.Block],
    model_key: str,
    h0: float,
) -> dict[str, object]:
    bounds = BOUNDS[model_key]
    objective = lambda value: score_blocks(blocks, model_key, value, h0).chi2
    result = minimize_scalar(
        objective,
        bounds=bounds,
        method="bounded",
        options={"xatol": 1e-10, "maxiter": 600},
    )
    candidates = [
        (float(result.fun), float(result.x)),
        (float(objective(bounds[0])), bounds[0]),
        (float(objective(bounds[1])), bounds[1]),
    ]
    chi2, theta = min(candidates)
    score = score_blocks(blocks, model_key, theta, h0)
    n = int(sum(block.n for block in blocks))
    k = len(blocks) + 1
    return {
        "model": model_key,
        "theta": theta,
        "chi2": float(score.chi2),
        "BIC": float(score.chi2 + k * math.log(n)),
        "N": n,
        "N_parameters": k,
        "parameters": score.parameters,
        "offsets": score.offsets,
        "optimisation_success": bool(result.success or theta in bounds),
    }


def posterior_quantile(grid: np.ndarray, weights: np.ndarray, probability: float) -> float:
    cumulative = np.cumsum(weights)
    cumulative[-1] = 1.0
    return float(np.interp(probability, cumulative, grid))


def build_posterior_grid(
    blocks: list[nested.base.Block],
    model_key: str,
    h0: float,
    adaptive_points: int,
    fallback_points: int,
    sigma_width: float,
    edge_delta_required: float,
) -> dict[str, object]:
    fit = fit_scalar(blocks, model_key, h0)
    optimum = float(fit["theta"])
    minimum = float(fit["chi2"])
    lower_bound, upper_bound = BOUNDS[model_key]
    parameter_range = upper_bound - lower_bound
    step = 1e-4 * parameter_range
    if optimum - step <= lower_bound or optimum + step >= upper_bound:
        local_sigma = parameter_range / 8.0
    else:
        f_plus = score_blocks(blocks, model_key, optimum + step, h0).chi2
        f_minus = score_blocks(blocks, model_key, optimum - step, h0).chi2
        curvature = (f_plus - 2.0 * minimum + f_minus) / step**2
        local_sigma = math.sqrt(2.0 / curvature) if curvature > 0.0 else parameter_range / 8.0
    lower = max(lower_bound, optimum - sigma_width * local_sigma)
    upper = min(upper_bound, optimum + sigma_width * local_sigma)
    grid = np.linspace(lower, upper, adaptive_points)
    chi2_values = np.asarray(
        [score_blocks(blocks, model_key, value, h0).chi2 for value in grid]
    )
    edge_deltas = (
        float(chi2_values[0] - chi2_values.min()),
        float(chi2_values[-1] - chi2_values.min()),
    )
    used_fallback = bool(min(edge_deltas) < edge_delta_required)
    if used_fallback:
        grid = np.linspace(lower_bound, upper_bound, fallback_points)
        chi2_values = np.asarray(
            [score_blocks(blocks, model_key, value, h0).chi2 for value in grid]
        )
        edge_deltas = (
            float(chi2_values[0] - chi2_values.min()),
            float(chi2_values[-1] - chi2_values.min()),
        )
    log_weights = -0.5 * (chi2_values - chi2_values.min())
    log_weights -= logsumexp(log_weights)
    weights = np.exp(log_weights)
    mean = float(np.sum(weights * grid))
    variance = float(np.sum(weights * (grid - mean) ** 2))
    summary = {
        "model": model_key,
        "parameter": "p" if model_key == "PATH_GAMMA2P3" else "Omega_m",
        "MLE": optimum,
        "posterior_mean": mean,
        "posterior_median": posterior_quantile(grid, weights, 0.5),
        "posterior_standard_deviation": math.sqrt(max(variance, 0.0)),
        "q16": posterior_quantile(grid, weights, 0.16),
        "q84": posterior_quantile(grid, weights, 0.84),
        "q2p5": posterior_quantile(grid, weights, 0.025),
        "q97p5": posterior_quantile(grid, weights, 0.975),
        "grid_lower": float(grid[0]),
        "grid_upper": float(grid[-1]),
        "grid_points": int(len(grid)),
        "edge_delta_chi2_lower": edge_deltas[0],
        "edge_delta_chi2_upper": edge_deltas[1],
        "used_full_bound_fallback": used_fallback,
        "local_sigma_estimate": float(local_sigma),
    }
    return {
        "fit": fit,
        "grid": grid,
        "chi2": chi2_values,
        "weights": weights,
        "summary": summary,
    }


def predictive_at_theta(
    partitions: list[PredictivePartition],
    model_key: str,
    theta: float,
    h0: float,
) -> tuple[float, float, dict[str, tuple[float, float]]]:
    parameters = model_parameters(model_key, theta)
    spec = model_spec(model_key)
    total_q = 0.0
    total_logpdf = 0.0
    release_values: dict[str, tuple[float, float]] = {}
    for partition in partitions:
        block = partition.block
        residual = block.mu - nested.model_mu_from_parameters(
            block.z, spec, parameters, h0
        )
        low_residual = residual[partition.low_indices]
        ones_low = np.ones(len(partition.low_indices))
        offset = float(
            ones_low @ partition.low_precision @ low_residual
            / partition.low_intercept_denominator
        )
        predictive_residual = (
            residual[partition.high_indices]
            - partition.conditional_projection @ low_residual
            - offset * partition.intercept_direction
        )
        whitened = np.linalg.solve(
            partition.predictive_cholesky, predictive_residual
        )
        q = float(whitened @ whitened)
        logpdf = float(
            -0.5
            * (
                q
                + partition.predictive_logdet
                + partition.n_high * math.log(2.0 * math.pi)
            )
        )
        total_q += q
        total_logpdf += logpdf
        release_values[block.label] = (q, logpdf)
    return total_q, total_logpdf, release_values


def integrate_predictive(
    posterior: dict[str, object],
    partitions: list[PredictivePartition],
    model_key: str,
    h0: float,
) -> tuple[dict[str, object], pd.DataFrame]:
    grid = np.asarray(posterior["grid"], dtype=float)
    weights = np.asarray(posterior["weights"], dtype=float)
    total_n = sum(partition.n_high for partition in partitions)
    rows = []
    release_labels = [partition.block.label for partition in partitions]
    release_n = {partition.block.label: partition.n_high for partition in partitions}
    for theta, weight in zip(grid, weights):
        total_q, total_logpdf, release_values = predictive_at_theta(
            partitions, model_key, float(theta), h0
        )
        row: dict[str, object] = {
            "model": model_key,
            "theta": float(theta),
            "posterior_weight": float(weight),
            "high_z_chi2": total_q,
            "high_z_logpdf": total_logpdf,
        }
        for label in release_labels:
            row[f"{label}_chi2"] = release_values[label][0]
            row[f"{label}_logpdf"] = release_values[label][1]
        rows.append(row)
    table = pd.DataFrame(rows)
    log_weights = np.log(np.clip(weights, 1e-300, None))
    mixture_logpdf = float(logsumexp(log_weights + table["high_z_logpdf"]))
    tail_area = float(
        np.sum(
            weights
            * chi2_distribution.sf(table["high_z_chi2"].to_numpy(), total_n)
        )
    )
    two_sided_tail = float(2.0 * min(tail_area, 1.0 - tail_area))
    summary: dict[str, object] = {
        "model": model_key,
        "N_high_z": total_n,
        "integrated_log_predictive_density": mixture_logpdf,
        "posterior_mean_high_z_chi2": float(
            np.sum(weights * table["high_z_chi2"])
        ),
        "posterior_predictive_upper_tail_area": tail_area,
        "posterior_predictive_two_sided_tail_area": two_sided_tail,
        "release_readouts": {},
    }
    for label in release_labels:
        logpdf_values = table[f"{label}_logpdf"].to_numpy()
        q_values = table[f"{label}_chi2"].to_numpy()
        release_tail = float(
            np.sum(weights * chi2_distribution.sf(q_values, release_n[label]))
        )
        summary["release_readouts"][label] = {
            "N_high_z": release_n[label],
            "integrated_log_predictive_density": float(
                logsumexp(log_weights + logpdf_values)
            ),
            "posterior_mean_high_z_chi2": float(np.sum(weights * q_values)),
            "posterior_predictive_upper_tail_area": release_tail,
            "posterior_predictive_two_sided_tail_area": float(
                2.0 * min(release_tail, 1.0 - release_tail)
            ),
        }
    return summary, table


def release_specific_posteriors(
    low_blocks: list[nested.base.Block],
    h0: float,
    protocol: dict[str, object],
) -> tuple[pd.DataFrame, float]:
    rows = []
    for block in low_blocks:
        posterior = build_posterior_grid(
            [block],
            "PATH_GAMMA2P3",
            h0,
            int(protocol["integration"]["adaptive_grid_points"]),
            int(protocol["integration"]["fallback_full_bound_points"]),
            float(protocol["integration"]["adaptive_half_width_local_sigma"]),
            float(protocol["integration"]["minimum_edge_delta_chi2"]),
        )
        summary = posterior["summary"]
        rows.append(
            {
                "release": block.label,
                "N_low_z": block.n,
                **summary,
            }
        )
    table = pd.DataFrame(rows)
    max_tension = 0.0
    for i in range(len(table)):
        for j in range(i):
            difference = abs(
                float(table.iloc[i]["posterior_mean"])
                - float(table.iloc[j]["posterior_mean"])
            )
            denominator = math.sqrt(
                float(table.iloc[i]["posterior_standard_deviation"]) ** 2
                + float(table.iloc[j]["posterior_standard_deviation"]) ** 2
            )
            if denominator > 0.0:
                max_tension = max(max_tension, difference / denominator)
    return table, float(max_tension)


def physical_readout(
    p_mean: float,
    omega_m_full: float,
) -> pd.DataFrame:
    redshifts = np.asarray([0.5, 1.0, 1.5, 2.26226], dtype=float)
    anchor = 0.1
    path = np.asarray(clock_path_integral(redshifts, GAMMA_FIXED, p_mean))
    path_anchor = float(clock_path_integral(anchor, GAMMA_FIXED, p_mean))
    lcdm = (1.0 + redshifts) * nested.base.lcdm_integral(redshifts, omega_m_full)
    lcdm_anchor = (1.0 + anchor) * float(
        nested.base.lcdm_integral(np.asarray([anchor]), omega_m_full)[0]
    )
    delta_mu = 5.0 * np.log10((path / path_anchor) / (lcdm / lcdm_anchor))
    return pd.DataFrame(
        {
            "z": redshifts,
            "c_over_c0": 1.0 + GAMMA_FIXED * redshifts,
            "abs_dz_dT_over_H0": (1.0 + redshifts) ** p_mean,
            "dimensionless_PATH_integral": path,
            "dimensionless_LCDM_luminosity_distance": lcdm,
            "anchored_delta_mu_PATH_minus_LCDM_mag": delta_mu,
        }
    )


def make_plot(
    path_posterior: dict[str, object],
    path_predictive_grid: pd.DataFrame,
    release_posteriors: pd.DataFrame,
    physical: pd.DataFrame,
    lcdm_predictive_summary: dict[str, object],
    output: Path,
) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(12.5, 8.8))
    grid = np.asarray(path_posterior["grid"])
    weights = np.asarray(path_posterior["weights"])
    density = weights / np.trapezoid(weights, grid)
    path_summary = path_posterior["summary"]
    axes[0, 0].plot(grid, density, color="#1971c2", linewidth=1.8)
    axes[0, 0].axvline(path_summary["posterior_mean"], color="#c92a2a", linewidth=1.2)
    axes[0, 0].axvspan(
        path_summary["q2p5"], path_summary["q97p5"], color="#74c0fc", alpha=0.25
    )
    axes[0, 0].set_xlabel("p at fixed gamma_c=2.3")
    axes[0, 0].set_ylabel("Low-z posterior density")
    axes[0, 0].set_title("Low-redshift PATH posterior")
    axes[0, 0].grid(alpha=0.2)

    axes[0, 1].plot(
        path_predictive_grid["theta"],
        path_predictive_grid["high_z_chi2"],
        color="#c92a2a",
        linewidth=1.6,
        label="PATH high-z chi-squared",
    )
    axes[0, 1].axhline(
        lcdm_predictive_summary["posterior_mean_high_z_chi2"],
        color="#2b8a3e",
        linestyle="--",
        label="LCDM posterior mean",
    )
    axes[0, 1].axvspan(
        path_summary["q16"], path_summary["q84"], color="#ffa8a8", alpha=0.25
    )
    axes[0, 1].set_xlabel("p")
    axes[0, 1].set_ylabel("Integrated-intercept high-z chi-squared")
    axes[0, 1].set_title("Predictive sensitivity to p")
    axes[0, 1].legend(fontsize=8)
    axes[0, 1].grid(alpha=0.2)

    positions = np.arange(len(release_posteriors))
    means = release_posteriors["posterior_mean"].to_numpy()
    lower = means - release_posteriors["q2p5"].to_numpy()
    upper = release_posteriors["q97p5"].to_numpy() - means
    axes[1, 0].errorbar(
        means,
        positions,
        xerr=np.vstack([lower, upper]),
        fmt="o",
        color="#1971c2",
        capsize=4,
    )
    axes[1, 0].axvline(path_summary["posterior_mean"], color="#c92a2a", linestyle="--")
    axes[1, 0].set_yticks(positions)
    axes[1, 0].set_yticklabels(release_posteriors["release"])
    axes[1, 0].set_xlabel("Release-specific p posterior mean and 95% interval")
    axes[1, 0].set_title("Low-redshift release stability")
    axes[1, 0].grid(axis="x", alpha=0.2)

    curve_z = np.linspace(0.02, 2.3, 500)
    p_mean = float(path_summary["posterior_mean"])
    omega_m = float(physical.attrs["omega_m_full"])
    path_curve = np.asarray(clock_path_integral(curve_z, GAMMA_FIXED, p_mean))
    path_anchor = float(clock_path_integral(0.1, GAMMA_FIXED, p_mean))
    lcdm_curve = (1.0 + curve_z) * nested.base.lcdm_integral(curve_z, omega_m)
    lcdm_anchor = (1.0 + 0.1) * float(
        nested.base.lcdm_integral(np.asarray([0.1]), omega_m)[0]
    )
    delta_mu = 5.0 * np.log10(
        (path_curve / path_anchor) / (lcdm_curve / lcdm_anchor)
    )
    axes[1, 1].plot(curve_z, delta_mu, color="#6741d9", linewidth=1.8)
    axes[1, 1].axhline(0.0, color="black", linewidth=0.8)
    axes[1, 1].axvline(0.5, color="#868e96", linestyle=":")
    axes[1, 1].scatter(
        physical["z"],
        physical["anchored_delta_mu_PATH_minus_LCDM_mag"],
        color="#c92a2a",
        zorder=3,
    )
    axes[1, 1].set_xlabel("Redshift")
    axes[1, 1].set_ylabel("Anchored Delta mu: PATH - Lambda-CDM (mag)")
    axes[1, 1].set_title("Fixed-gamma shape target")
    axes[1, 1].grid(alpha=0.2)
    figure.tight_layout()
    figure.savefig(output, dpi=180)
    plt.close(figure)


def safe_json(value):
    if isinstance(value, dict):
        return {str(key): safe_json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [safe_json(item) for item in value]
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return None if not np.isfinite(value) else float(value)
    return value


def write_report(
    output: Path,
    full_fits: dict[str, dict[str, object]],
    low_fits: dict[str, dict[str, object]],
    posterior_summaries: dict[str, dict[str, object]],
    predictive_summaries: dict[str, dict[str, object]],
    release_posteriors: pd.DataFrame,
    physical: pd.DataFrame,
    gates: dict[str, bool],
    metrics: dict[str, float],
) -> None:
    path_post = posterior_summaries["PATH_GAMMA2P3"]
    lines = [
        "# Fixed-gamma=2.3 PLAMB conditional analysis",
        "",
        f"**Analysis date:** {DATE_TAG}  ",
        "**Status:** outcome-aware conditional target; cannot promote a PLAMB claim  ",
        "**Primary split:** fit z<0.5, predict z>=0.5 with released covariance",
        "",
        "## Headline result",
        "",
        f"With gamma_c fixed at 2.3, the full-sample PATH fit gives p="
        f"{full_fits['PATH_GAMMA2P3']['theta']:.6f} and Delta BIC="
        f"{metrics['full_delta_BIC_PATH_minus_LCDM']:+.6f} relative to flat Lambda-CDM. "
        "The two curves are statistically indistinguishable by this BIC threshold when gamma_c "
        "is treated as externally fixed.",
        "",
        f"The z<0.5 PATH posterior gives p={path_post['posterior_mean']:.6f} +/- "
        f"{path_post['posterior_standard_deviation']:.6f}, with 95% interval "
        f"[{path_post['q2p5']:.6f}, {path_post['q97p5']:.6f}]. Its low-redshift "
        f"Delta BIC is {metrics['low_delta_BIC_PATH_minus_LCDM']:+.6f}.",
        "",
        f"After integrating p and each release intercept, the high-redshift predictive "
        f"Delta log density is {metrics['high_z_delta_log_predictive_density']:+.6f}, "
        "where positive values favour PATH. This is a negligible preference for Lambda-CDM. "
        "The PATH upper-tail area is "
        f"{predictive_summaries['PATH_GAMMA2P3']['posterior_predictive_upper_tail_area']:.6f}.",
        "",
        "This is a conditional theory target, not independent support: gamma_c=2.3 was chosen "
        "after inspecting these same high-redshift outcomes.",
        "",
        "## Equal-parameter fits",
        "",
        "| Sample | Model | Shape parameter | chi-squared | BIC | Delta BIC vs LCDM |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for sample, fits in (("full", full_fits), ("z<0.5", low_fits)):
        lcdm_bic = float(fits["LCDM"]["BIC"])
        for model_key in ("PATH_GAMMA2P3", "LCDM"):
            fit = fits[model_key]
            lines.append(
                f"| {sample} | `{model_key}` | {fit['theta']:.6f} | {fit['chi2']:.6f} | "
                f"{fit['BIC']:.6f} | {float(fit['BIC'] - lcdm_bic):+.6f} |"
            )
    lines.extend(
        [
            "",
            "Both cells contain one shape parameter and one release intercept per included "
            "release. The BIC differences therefore equal the chi-squared differences.",
            "",
            "## Low-redshift posteriors",
            "",
            "| Model | Mean | Standard deviation | 68% interval | 95% interval | Edge Delta chi-squared |",
            "|---|---:|---:|---|---|---|",
        ]
    )
    for model_key in ("PATH_GAMMA2P3", "LCDM"):
        row = posterior_summaries[model_key]
        lines.append(
            f"| `{model_key}` | {row['posterior_mean']:.6f} | "
            f"{row['posterior_standard_deviation']:.6f} | "
            f"[{row['q16']:.6f}, {row['q84']:.6f}] | "
            f"[{row['q2p5']:.6f}, {row['q97p5']:.6f}] | "
            f"{row['edge_delta_chi2_lower']:.2f}, {row['edge_delta_chi2_upper']:.2f} |"
        )
    lines.extend(
        [
            "",
            "## Integrated high-redshift prediction",
            "",
            "| Model | N high-z | Log predictive density | Mean predictive chi-squared | Upper-tail area | Two-sided tail |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for model_key in ("PATH_GAMMA2P3", "LCDM"):
        row = predictive_summaries[model_key]
        lines.append(
            f"| `{model_key}` | {row['N_high_z']} | "
            f"{row['integrated_log_predictive_density']:.6f} | "
            f"{row['posterior_mean_high_z_chi2']:.6f} | "
            f"{row['posterior_predictive_upper_tail_area']:.6f} | "
            f"{row['posterior_predictive_two_sided_tail_area']:.6f} |"
        )
    lines.extend(
        [
            "",
            "### Release-level predictive readout",
            "",
            "| Release | PATH log density | LCDM log density | Delta log density | PATH upper tail | PATH two-sided tail |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    path_releases = predictive_summaries["PATH_GAMMA2P3"]["release_readouts"]
    lcdm_releases = predictive_summaries["LCDM"]["release_readouts"]
    for label in path_releases:
        path_row = path_releases[label]
        lcdm_row = lcdm_releases[label]
        lines.append(
            f"| {label} | {path_row['integrated_log_predictive_density']:.6f} | "
            f"{lcdm_row['integrated_log_predictive_density']:.6f} | "
            f"{path_row['integrated_log_predictive_density'] - lcdm_row['integrated_log_predictive_density']:+.6f} | "
            f"{path_row['posterior_predictive_upper_tail_area']:.6f} | "
            f"{path_row['posterior_predictive_two_sided_tail_area']:.6f} |"
        )
    lines.extend(
        [
            "",
            "Joint mixture log density is the registered comparison. Release-level mixture "
            "densities share the same joint p posterior and are descriptive; they are not "
            "algebraically additive after parameter integration.",
            "The registered gate uses only a minimum upper-tail area. The additional two-sided "
            "readout is a post hoc calibration warning: values near zero indicate chi-squared "
            "that is unusually high or unusually low. Both cosmologies have unusually low joint "
            "quadratic scores, with Pantheon+ providing the strongest effect.",
            "",
            "## Release-specific p stability",
            "",
            "| Release | N low-z | Mean p | Standard deviation | 95% interval |",
            "|---|---:|---:|---:|---|",
        ]
    )
    for row in release_posteriors.itertuples():
        lines.append(
            f"| {row.release} | {row.N_low_z} | {row.posterior_mean:.6f} | "
            f"{row.posterior_standard_deviation:.6f} | "
            f"[{row.q2p5:.6f}, {row.q97p5:.6f}] |"
        )
    lines.extend(
        [
            "",
            f"Maximum pairwise release tension: {metrics['maximum_pairwise_release_p_tension_sigma']:.6f} standard deviations.",
            "",
            "## Physical consequences",
            "",
            f"Using the joint low-redshift posterior mean p={path_post['posterior_mean']:.6f}:",
            "",
            "| z | c(z)/c0 | abs(dz/dT)/H0 | PATH integral | Anchored Delta mu vs LCDM (mag) |",
            "|---:|---:|---:|---:|---:|",
        ]
    )
    for row in physical.itertuples():
        lines.append(
            f"| {row.z:.6f} | {row.c_over_c0:.6f} | {row.abs_dz_dT_over_H0:.6f} | "
            f"{row.dimensionless_PATH_integral:.6f} | "
            f"{row.anchored_delta_mu_PATH_minus_LCDM_mag:+.6f} |"
        )
    lines.extend(
        [
            "",
            "The light-speed and clock-rate ratios are implications of the phenomenological "
            "parameterisation. They are not independently measured quantities. In particular, "
            "the large c(z)/c0 values require a separate covariant derivation and consistency "
            "checks against non-supernova observables.",
            "",
            "## Conditional gates",
            "",
            "| Gate | Pass |",
            "|---|---|",
        ]
    )
    for name, passed in gates.items():
        lines.append(f"| `{name}` | **{bool(passed)}** |")
    statistical_gates = [value for key, value in gates.items() if key != "independent_gamma_derivation"]
    lines.extend(
        [
            "",
            f"Statistical gates passed: {sum(statistical_gates)}/{len(statistical_gates)}.",
            "",
            "**Claim decision: DO NOT PROMOTE.** The independent-gamma gate is false by "
            "construction, and the value 2.3 remains an outcome-aware target for theory rather "
            "than an externally predicted constant.",
            "",
            "## Reproduction",
            "",
            "```powershell",
            "python github_export/code/rawmu/analyze_plamb_gamma2p3_2026_07_19.py",
            "```",
        ]
    )
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_manifest(outdir: Path, protocol_path: Path) -> None:
    paths = [
        SCRIPT_PATH,
        protocol_path,
        nested.base.PANTHEON_DATA,
        nested.base.PANTHEON_TOTAL,
        nested.base.DES_DATA,
        nested.base.DES_TOTAL,
        nested.base.UNION_PATH,
    ]
    paths.extend(
        path
        for path in sorted(outdir.iterdir())
        if path.is_file() and "manifest" not in path.name
    )
    rows = []
    for path in paths:
        path = Path(path)
        rows.append(
            {
                "path": str(path),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    pd.DataFrame(rows).drop_duplicates("path").to_csv(
        outdir / f"plamb_gamma2p3_manifest_{DATE_TAG}.csv", index=False
    )


def run_self_tests() -> None:
    z = np.linspace(0.04, 1.05, 30)
    p_true = 0.8
    mu = nested.clock_mu_alpha0(z, GAMMA_FIXED, p_true, 67.5) + 0.14
    correlation = 0.2 ** np.abs(
        np.subtract.outer(np.arange(len(z)), np.arange(len(z)))
    )
    covariance = 0.01**2 * correlation
    block = nested.base.make_block(
        "synthetic",
        z,
        mu,
        covariance,
        np.asarray(["S"] * len(z), dtype=object),
        np.arange(len(z)),
        "synthetic",
    )
    low = make_low_blocks([block])
    fit = fit_scalar(low, "PATH_GAMMA2P3", 67.5)
    assert abs(float(fit["theta"]) - p_true) < 1e-6
    partitions = make_predictive_partitions([block])
    q, logpdf, _release = predictive_at_theta(
        partitions, "PATH_GAMMA2P3", p_true, 67.5
    )
    assert q < 1e-12 and np.isfinite(logpdf)
    assert partitions[0].cholesky_jitter == 0.0
    print("All fixed-gamma=2.3 self-tests passed.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--protocol", type=Path, default=DEFAULT_PROTOCOL)
    parser.add_argument("--h0", type=float, default=67.5)
    parser.add_argument("--min-survey-n", type=int, default=20)
    parser.add_argument("--self-test-only", action="store_true")
    cli = parser.parse_args()

    run_self_tests()
    if cli.self_test_only:
        return
    protocol = json.loads(cli.protocol.read_text(encoding="utf-8"))
    if protocol["claim_status"] != "post_hoc_conditional_target_cannot_promote":
        raise ValueError("Protocol claim boundary changed")
    cli.outdir.mkdir(parents=True, exist_ok=True)
    variants, metadata = nested.base.load_blocks(cli.min_survey_n)
    blocks = variants["released_total_primary"]
    if sum(block.n for block in blocks) != protocol["sample"]["expected_N"]:
        raise ValueError("Locked sample size changed")
    low_blocks = make_low_blocks(blocks)
    partitions = make_predictive_partitions(blocks)

    full_fits = {
        model_key: fit_scalar(blocks, model_key, cli.h0)
        for model_key in ("PATH_GAMMA2P3", "LCDM")
    }
    low_fits = {
        model_key: fit_scalar(low_blocks, model_key, cli.h0)
        for model_key in ("PATH_GAMMA2P3", "LCDM")
    }
    print(
        f"Full PATH p={full_fits['PATH_GAMMA2P3']['theta']:.6f} "
        f"DeltaBIC={full_fits['PATH_GAMMA2P3']['BIC'] - full_fits['LCDM']['BIC']:+.6f}",
        flush=True,
    )

    posteriors = {}
    predictive_summaries = {}
    predictive_tables = []
    posterior_tables = []
    for model_key in ("PATH_GAMMA2P3", "LCDM"):
        print(f"Building low-z posterior for {model_key}", flush=True)
        posterior = build_posterior_grid(
            low_blocks,
            model_key,
            cli.h0,
            int(protocol["integration"]["adaptive_grid_points"]),
            int(protocol["integration"]["fallback_full_bound_points"]),
            float(protocol["integration"]["adaptive_half_width_local_sigma"]),
            float(protocol["integration"]["minimum_edge_delta_chi2"]),
        )
        posteriors[model_key] = posterior
        posterior_tables.append(
            pd.DataFrame(
                {
                    "model": model_key,
                    "theta": posterior["grid"],
                    "chi2": posterior["chi2"],
                    "posterior_weight": posterior["weights"],
                }
            )
        )
        print(f"Integrating high-z prediction for {model_key}", flush=True)
        predictive_summary, predictive_table = integrate_predictive(
            posterior, partitions, model_key, cli.h0
        )
        predictive_summaries[model_key] = predictive_summary
        predictive_tables.append(predictive_table)

    posterior_grid_table = pd.concat(posterior_tables, ignore_index=True)
    posterior_grid_table.to_csv(
        cli.outdir / f"plamb_gamma2p3_posterior_grids_{DATE_TAG}.csv", index=False
    )
    predictive_grid_table = pd.concat(predictive_tables, ignore_index=True)
    predictive_grid_table.to_csv(
        cli.outdir / f"plamb_gamma2p3_predictive_grids_{DATE_TAG}.csv", index=False
    )

    release_posteriors, max_release_tension = release_specific_posteriors(
        low_blocks, cli.h0, protocol
    )
    release_posteriors.to_csv(
        cli.outdir / f"plamb_gamma2p3_release_posteriors_{DATE_TAG}.csv", index=False
    )
    physical = physical_readout(
        float(posteriors["PATH_GAMMA2P3"]["summary"]["posterior_mean"]),
        float(full_fits["LCDM"]["theta"]),
    )
    physical.attrs["omega_m_full"] = float(full_fits["LCDM"]["theta"])
    physical.to_csv(
        cli.outdir / f"plamb_gamma2p3_physical_readout_{DATE_TAG}.csv", index=False
    )

    posterior_summaries = {
        model_key: posteriors[model_key]["summary"]
        for model_key in posteriors
    }
    full_delta_bic = float(
        full_fits["PATH_GAMMA2P3"]["BIC"] - full_fits["LCDM"]["BIC"]
    )
    low_delta_bic = float(
        low_fits["PATH_GAMMA2P3"]["BIC"] - low_fits["LCDM"]["BIC"]
    )
    delta_log_predictive = float(
        predictive_summaries["PATH_GAMMA2P3"]["integrated_log_predictive_density"]
        - predictive_summaries["LCDM"]["integrated_log_predictive_density"]
    )
    path_release_tails = [
        value["posterior_predictive_upper_tail_area"]
        for value in predictive_summaries["PATH_GAMMA2P3"]["release_readouts"].values()
    ]
    path_summary = posterior_summaries["PATH_GAMMA2P3"]
    lower_bound, upper_bound = BOUNDS["PATH_GAMMA2P3"]
    bound_fraction = min(
        (path_summary["posterior_mean"] - lower_bound) / (upper_bound - lower_bound),
        (upper_bound - path_summary["posterior_mean"]) / (upper_bound - lower_bound),
    )
    edge_minimum = min(
        path_summary["edge_delta_chi2_lower"],
        path_summary["edge_delta_chi2_upper"],
    )
    gates = {
        "full_absolute_delta_BIC_at_most_2": bool(abs(full_delta_bic) <= 2.0),
        "low_z_absolute_delta_BIC_at_most_2": bool(abs(low_delta_bic) <= 2.0),
        "high_z_delta_log_predictive_density_nonnegative": bool(
            delta_log_predictive >= 0.0
        ),
        "joint_PATH_tail_area_at_least_0p05": bool(
            predictive_summaries["PATH_GAMMA2P3"][
                "posterior_predictive_upper_tail_area"
            ]
            >= 0.05
        ),
        "every_release_PATH_tail_area_at_least_0p05": bool(
            min(path_release_tails) >= 0.05
        ),
        "maximum_pairwise_release_p_tension_at_most_2_sigma": bool(
            max_release_tension <= 2.0
        ),
        "p_posterior_interior_and_grid_edges_delta_chi2_at_least_30": bool(
            bound_fraction >= 0.01 and edge_minimum >= 30.0
        ),
        "independent_gamma_derivation": False,
    }
    metrics = {
        "full_delta_BIC_PATH_minus_LCDM": full_delta_bic,
        "low_delta_BIC_PATH_minus_LCDM": low_delta_bic,
        "high_z_delta_log_predictive_density": delta_log_predictive,
        "maximum_pairwise_release_p_tension_sigma": max_release_tension,
        "PATH_posterior_bound_fraction": float(bound_fraction),
        "PATH_minimum_grid_edge_delta_chi2": float(edge_minimum),
        "maximum_predictive_cholesky_jitter": float(
            max(partition.cholesky_jitter for partition in partitions)
        ),
    }
    fit_table = []
    for sample, fits in (("full", full_fits), ("low_z", low_fits)):
        for model_key, fit in fits.items():
            fit_table.append(
                {
                    "sample": sample,
                    "model": model_key,
                    "theta": fit["theta"],
                    "chi2": fit["chi2"],
                    "BIC": fit["BIC"],
                    "N": fit["N"],
                    "parameters": json.dumps(fit["parameters"], sort_keys=True),
                    "offsets": json.dumps(fit["offsets"], sort_keys=True),
                }
            )
    pd.DataFrame(fit_table).to_csv(
        cli.outdir / f"plamb_gamma2p3_equal_parameter_fits_{DATE_TAG}.csv",
        index=False,
    )
    predictive_rows = []
    for model_key, summary in predictive_summaries.items():
        predictive_rows.append(
            {
                "model": model_key,
                "release": "JOINT",
                "N_high_z": summary["N_high_z"],
                "integrated_log_predictive_density": summary[
                    "integrated_log_predictive_density"
                ],
                "posterior_mean_high_z_chi2": summary[
                    "posterior_mean_high_z_chi2"
                ],
                "posterior_predictive_upper_tail_area": summary[
                    "posterior_predictive_upper_tail_area"
                ],
            }
        )
        for release, values in summary["release_readouts"].items():
            predictive_rows.append(
                {
                    "model": model_key,
                    "release": release,
                    **values,
                }
            )
    pd.DataFrame(predictive_rows).to_csv(
        cli.outdir / f"plamb_gamma2p3_predictive_summary_{DATE_TAG}.csv",
        index=False,
    )
    summary = {
        "analysis_date": DATE_TAG,
        "claim_status": protocol["claim_status"],
        "selection_disclosure": protocol["selection_disclosure"],
        "sample_metadata": metadata,
        "full_fits": full_fits,
        "low_z_fits": low_fits,
        "posterior_summaries": posterior_summaries,
        "predictive_summaries": predictive_summaries,
        "metrics": metrics,
        "gates": gates,
        "statistical_gates_passed": int(
            sum(value for key, value in gates.items() if key != "independent_gamma_derivation")
        ),
        "statistical_gates_total": 7,
        "claim_promotion": False,
    }
    (cli.outdir / f"plamb_gamma2p3_summary_{DATE_TAG}.json").write_text(
        json.dumps(safe_json(summary), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_report(
        cli.outdir / f"plamb_gamma2p3_report_{DATE_TAG}.md",
        full_fits,
        low_fits,
        posterior_summaries,
        predictive_summaries,
        release_posteriors,
        physical,
        gates,
        metrics,
    )
    path_predictive_grid = predictive_grid_table[
        predictive_grid_table["model"] == "PATH_GAMMA2P3"
    ]
    make_plot(
        posteriors["PATH_GAMMA2P3"],
        path_predictive_grid,
        release_posteriors,
        physical,
        predictive_summaries["LCDM"],
        cli.outdir / f"plamb_gamma2p3_readout_{DATE_TAG}.png",
    )
    write_manifest(cli.outdir, cli.protocol)
    print(
        f"Saved report: {cli.outdir / f'plamb_gamma2p3_report_{DATE_TAG}.md'}",
        flush=True,
    )


if __name__ == "__main__":
    main()
