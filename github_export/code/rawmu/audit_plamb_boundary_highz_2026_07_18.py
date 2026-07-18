#!/usr/bin/env python3
"""Audit the PLAMB gamma boundary and conditional high-redshift failure."""

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
from scipy.optimize import minimize, minimize_scalar


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[3]
RAWMU_CODE = SCRIPT_PATH.parent
if str(RAWMU_CODE) not in sys.path:
    sys.path.insert(0, str(RAWMU_CODE))

import run_plamb_nested_clock_fit_2026_07_18 as nested  # noqa: E402


DATE_TAG = "2026-07-18"
DEFAULT_OUTDIR = (
    REPO_ROOT
    / "github_export"
    / "results"
    / DATE_TAG
    / "rawmu_plamb_clock_fit"
    / "boundary_highz_followup"
)
DEFAULT_PROTOCOL = DEFAULT_OUTDIR / f"plamb_boundary_highz_protocol_{DATE_TAG}.json"
P_BOUNDS = (-0.5, 2.5)
ALPHA_BOUNDS = (-0.5, 1.5)
GAMMA_CAPS = (2.0, 4.0, 8.0)
TRAINING_LIMIT = 0.5


@dataclass
class WhiteningPartition:
    block: nested.base.Block
    low_indices: np.ndarray
    high_indices: np.ndarray
    conditional_projection: np.ndarray
    conditional_cholesky: np.ndarray
    cholesky_jitter: float

    @property
    def n_high(self) -> int:
        return int(len(self.high_indices))


def registered_gamma_grid() -> np.ndarray:
    first = np.round(np.arange(0.0, 2.0 + 0.025, 0.05), 10)
    second = np.round(np.arange(2.1, 4.0 + 0.05, 0.1), 10)
    third = np.round(np.arange(4.25, 8.0 + 0.125, 0.25), 10)
    grid = np.unique(np.concatenate([first, second, third]))
    if len(grid) != 77 or grid[0] != 0.0 or grid[-1] != 8.0:
        raise AssertionError("Registered gamma grid is malformed")
    return grid


def redshift_bin(z: float) -> str:
    if 0.5 <= z < 0.65:
        return "z0p50_0p65"
    if 0.65 <= z < 0.8:
        return "z0p65_0p80"
    if 0.8 <= z < 1.0:
        return "z0p80_1p00"
    if z >= 1.0:
        return "z1p00_plus"
    raise ValueError(f"Redshift {z} is not in the held-out range")


def make_low_blocks(blocks: list[nested.base.Block]) -> list[nested.base.Block]:
    low_blocks = []
    for block in blocks:
        keep = block.z < TRAINING_LIMIT
        if not np.any(keep) or np.all(keep):
            raise ValueError(f"{block.label} does not span the z=0.5 split")
        low_blocks.append(nested.base.subset_block(block, keep, "train_z_lt_0p5"))
    return low_blocks


def stable_cholesky(matrix: np.ndarray) -> tuple[np.ndarray, float]:
    matrix = 0.5 * (matrix + matrix.T)
    scale = max(float(np.median(np.diag(matrix))), 1e-16)
    for multiplier in (0.0, 1e-12, 1e-10, 1e-8):
        jitter = multiplier * scale
        try:
            return np.linalg.cholesky(matrix + jitter * np.eye(len(matrix))), jitter
        except np.linalg.LinAlgError:
            continue
    raise np.linalg.LinAlgError("Conditional covariance is not positive definite")


def make_whitening_partitions(
    blocks: list[nested.base.Block],
) -> list[WhiteningPartition]:
    partitions = []
    for block in blocks:
        low_indices = np.flatnonzero(block.z < TRAINING_LIMIT)
        high_unsorted = np.flatnonzero(block.z >= TRAINING_LIMIT)
        order = np.lexsort(
            (block.source_indices[high_unsorted], block.z[high_unsorted])
        )
        high_indices = high_unsorted[order]
        c_ll = block.covariance[np.ix_(low_indices, low_indices)]
        c_lh = block.covariance[np.ix_(low_indices, high_indices)]
        c_hl = c_lh.T
        c_hh = block.covariance[np.ix_(high_indices, high_indices)]
        solved = np.linalg.solve(c_ll, c_lh)
        conditional_covariance = c_hh - c_hl @ solved
        conditional_covariance = 0.5 * (
            conditional_covariance + conditional_covariance.T
        )
        cholesky, jitter = stable_cholesky(conditional_covariance)
        partitions.append(
            WhiteningPartition(
                block=block,
                low_indices=low_indices,
                high_indices=high_indices,
                conditional_projection=solved.T,
                conditional_cholesky=cholesky,
                cholesky_jitter=jitter,
            )
        )
    return partitions


def profile_p_at_gamma(
    blocks: list[nested.base.Block],
    model_key: str,
    gamma_c: float,
    h0: float,
) -> dict[str, object]:
    if model_key not in {"PATH_FREE", "GENERAL_FREE"}:
        raise ValueError(model_key)
    alpha_fixed = 0.0 if model_key == "PATH_FREE" else None

    def objective(p: float) -> float:
        return nested.score_clock_profiled_alpha(
            blocks, gamma_c, p, alpha_fixed, h0
        ).chi2

    result = minimize_scalar(
        objective,
        bounds=P_BOUNDS,
        method="bounded",
        options={"xatol": 1e-9, "maxiter": 500},
    )
    candidates = [
        (float(result.fun), float(result.x)),
        (float(objective(P_BOUNDS[0])), P_BOUNDS[0]),
        (float(objective(P_BOUNDS[1])), P_BOUNDS[1]),
    ]
    chi2, p = min(candidates)
    score = nested.score_clock_profiled_alpha(
        blocks, gamma_c, p, alpha_fixed, h0
    )
    return {
        "model": model_key,
        "chi2": float(score.chi2),
        "parameters": score.parameters,
        "offsets": score.offsets,
        "optimisation_success": bool(result.success or p in P_BOUNDS),
        "method": "fixed_gamma_bounded_scalar_p",
        "N": int(sum(block.n for block in blocks)),
    }


def fit_cap(
    blocks: list[nested.base.Block],
    model_key: str,
    gamma_cap: float,
    h0: float,
    grid_starts: pd.DataFrame | None = None,
    extra_start: dict[str, float] | None = None,
) -> dict[str, object]:
    alpha_fixed = 0.0 if model_key == "PATH_FREE" else None

    def objective(values: np.ndarray) -> float:
        gamma_c, p = (float(values[0]), float(values[1]))
        return nested.score_clock_profiled_alpha(
            blocks, gamma_c, p, alpha_fixed, h0
        ).chi2

    starts: list[np.ndarray] = []
    if grid_starts is not None and len(grid_starts):
        eligible = grid_starts[grid_starts["gamma_c"] <= gamma_cap + 1e-12]
        best = eligible.loc[eligible["chi2"].idxmin()]
        starts.append(np.asarray([best["gamma_c"], best["p"]], dtype=float))
    if extra_start is not None:
        starts.append(
            np.asarray([extra_start["gamma_c"], extra_start["p"]], dtype=float)
        )
    starts.extend(
        [
            np.asarray([min(1.0, gamma_cap), 0.0]),
            np.asarray([0.5 * gamma_cap, 0.5]),
            np.asarray([gamma_cap, 0.75]),
            np.asarray([min(1.8, gamma_cap), 1.5]),
        ]
    )
    bounds = [(0.0, gamma_cap), P_BOUNDS]
    solutions: list[tuple[float, np.ndarray, bool, str]] = []
    seen: set[tuple[float, float]] = set()
    for start in starts:
        start = np.asarray(
            [
                np.clip(start[0], 0.0, gamma_cap),
                np.clip(start[1], *P_BOUNDS),
            ],
            dtype=float,
        )
        key = (round(float(start[0]), 10), round(float(start[1]), 10))
        if key in seen:
            continue
        seen.add(key)
        result = minimize(
            objective,
            start,
            method="L-BFGS-B",
            bounds=bounds,
            options={"ftol": 1e-12, "gtol": 1e-8, "maxiter": 600},
        )
        solutions.append(
            (
                float(result.fun),
                np.asarray(result.x, dtype=float),
                bool(result.success),
                str(result.message),
            )
        )
    solutions = [item for item in solutions if np.isfinite(item[0])]
    if len(solutions) < 2:
        raise RuntimeError("Fewer than two finite cap-fit solutions")
    solutions.sort(key=lambda item: item[0])
    best_chi2, best_values, best_success, best_message = solutions[0]
    gamma_c, p = (float(best_values[0]), float(best_values[1]))
    score = nested.score_clock_profiled_alpha(
        blocks, gamma_c, p, alpha_fixed, h0
    )
    n = int(sum(block.n for block in blocks))
    shape_k = 2 if model_key == "PATH_FREE" else 3
    return {
        "model": model_key,
        "gamma_cap": float(gamma_cap),
        "chi2": float(score.chi2),
        "BIC": float(score.chi2 + (len(blocks) + shape_k) * math.log(n)),
        "parameters": score.parameters,
        "offsets": score.offsets,
        "N": n,
        "shape_k": shape_k,
        "optimisation_success": bool(best_success),
        "optimisation_message": best_message,
        "best_two_chi2_disagreement": float(abs(solutions[1][0] - best_chi2)),
        "all_objectives": [float(item[0]) for item in solutions],
        "boundary_fraction_within_cap": float(
            min(gamma_c / gamma_cap, (gamma_cap - gamma_c) / gamma_cap)
        ),
    }


def fit_to_nested_dict(fit: dict[str, object]) -> dict[str, object]:
    return {
        "parameters": fit["parameters"],
        "offsets": fit["offsets"],
    }


def innovations_for_fit(
    partitions: list[WhiteningPartition],
    model_key: str,
    fit: dict[str, object],
    h0: float,
    fit_label: str,
) -> pd.DataFrame:
    spec = nested.MODEL_SPECS[model_key]
    parameters = fit["parameters"]
    offsets = fit["offsets"]
    assert isinstance(parameters, dict) and isinstance(offsets, dict)
    rows = []
    for partition in partitions:
        block = partition.block
        model_mu = nested.model_mu_from_parameters(block.z, spec, parameters, h0)
        residual = block.mu - model_mu - float(offsets[block.label])
        conditional_residual = (
            residual[partition.high_indices]
            - partition.conditional_projection @ residual[partition.low_indices]
        )
        whitened = np.linalg.solve(
            partition.conditional_cholesky, conditional_residual
        )
        for order, (index, innovation) in enumerate(
            zip(partition.high_indices, whitened), start=1
        ):
            rows.append(
                {
                    "fit_label": fit_label,
                    "model": model_key,
                    "release": block.label,
                    "ordered_innovation_index": order,
                    "source_index": int(block.source_indices[index]),
                    "z": float(block.z[index]),
                    "redshift_bin": redshift_bin(float(block.z[index])),
                    "survey": str(block.survey[index]),
                    "whitened_residual": float(innovation),
                    "chi2_contribution": float(innovation**2),
                    "conditional_cholesky_jitter": float(partition.cholesky_jitter),
                }
            )
    return pd.DataFrame(rows)


def conditional_score_from_fit(
    partitions: list[WhiteningPartition],
    model_key: str,
    fit: dict[str, object],
    h0: float,
) -> float:
    values = innovations_for_fit(partitions, model_key, fit, h0, "temporary")
    return float(values["chi2_contribution"].sum())


def marginal_high_score(
    blocks: list[nested.base.Block],
    model_key: str,
    fit: dict[str, object],
    h0: float,
) -> float:
    spec = nested.MODEL_SPECS[model_key]
    parameters = fit["parameters"]
    offsets = fit["offsets"]
    assert isinstance(parameters, dict) and isinstance(offsets, dict)
    total = 0.0
    for block in blocks:
        high = np.flatnonzero(block.z >= TRAINING_LIMIT)
        covariance = block.covariance[np.ix_(high, high)]
        precision = nested.base.stable_inverse(covariance)
        residual = (
            block.mu[high]
            - nested.model_mu_from_parameters(block.z[high], spec, parameters, h0)
            - float(offsets[block.label])
        )
        total += float(residual @ precision @ residual)
    return total


def diagonal_blocks(
    blocks: list[nested.base.Block],
) -> tuple[list[nested.base.Block], list[nested.base.Block]]:
    originals = []
    low = []
    for block in blocks:
        covariance = np.diag(np.diag(block.covariance))
        diagonal = nested.base.make_block(
            block.label,
            block.z,
            block.mu,
            covariance,
            block.survey,
            block.source_indices,
            "released_total_diagonal",
        )
        originals.append(diagonal)
        low.append(
            nested.base.subset_block(
                diagonal, diagonal.z < TRAINING_LIMIT, "train_z_lt_0p5"
            )
        )
    return originals, low


def diagonal_high_score(
    blocks: list[nested.base.Block],
    model_key: str,
    fit: dict[str, object],
    h0: float,
) -> float:
    spec = nested.MODEL_SPECS[model_key]
    parameters = fit["parameters"]
    offsets = fit["offsets"]
    assert isinstance(parameters, dict) and isinstance(offsets, dict)
    total = 0.0
    for block in blocks:
        high = block.z >= TRAINING_LIMIT
        residual = (
            block.mu[high]
            - nested.model_mu_from_parameters(block.z[high], spec, parameters, h0)
            - float(offsets[block.label])
        )
        variances = np.diag(block.covariance)[high]
        total += float(np.sum(residual**2 / variances))
    return total


def serialise(value: object) -> str:
    if isinstance(value, dict):
        return json.dumps(
            {str(key): float(item) for key, item in value.items()}, sort_keys=True
        )
    return json.dumps(value)


def build_profile_scan(
    blocks: list[nested.base.Block],
    low_blocks: list[nested.base.Block],
    partitions: list[WhiteningPartition],
    lcdm_high_chi2: float,
    h0: float,
) -> tuple[pd.DataFrame, dict[tuple[str, str, float], dict[str, object]]]:
    rows = []
    fits: dict[tuple[str, str, float], dict[str, object]] = {}
    grid = registered_gamma_grid()
    for model_key in ("PATH_FREE", "GENERAL_FREE"):
        for position, gamma_c in enumerate(grid, start=1):
            full_fit = profile_p_at_gamma(blocks, model_key, gamma_c, h0)
            low_fit = profile_p_at_gamma(low_blocks, model_key, gamma_c, h0)
            fits[(model_key, "full", float(gamma_c))] = full_fit
            fits[(model_key, "low_z", float(gamma_c))] = low_fit
            high_chi2 = conditional_score_from_fit(
                partitions, model_key, low_fit, h0
            )
            rows.append(
                {
                    "model": model_key,
                    "gamma_c": float(gamma_c),
                    "full_chi2": float(full_fit["chi2"]),
                    "full_p": float(full_fit["parameters"]["p"]),
                    "full_alpha": float(full_fit["parameters"]["alpha"]),
                    "low_z_chi2": float(low_fit["chi2"]),
                    "low_z_p": float(low_fit["parameters"]["p"]),
                    "low_z_alpha": float(low_fit["parameters"]["alpha"]),
                    "high_z_conditional_chi2": high_chi2,
                    "high_z_delta_chi2_vs_LCDM": high_chi2 - lcdm_high_chi2,
                    "optimisation_success": bool(
                        full_fit["optimisation_success"]
                        and low_fit["optimisation_success"]
                    ),
                }
            )
            if position % 10 == 0 or position == len(grid):
                print(
                    f"{model_key} profile {position}/{len(grid)} gamma={gamma_c:.2f}",
                    flush=True,
                )
    return pd.DataFrame(rows), fits


def build_cap_fits(
    blocks: list[nested.base.Block],
    low_blocks: list[nested.base.Block],
    partitions: list[WhiteningPartition],
    profile: pd.DataFrame,
    lcdm_full: dict[str, object],
    lcdm_low: dict[str, object],
    h0: float,
) -> tuple[pd.DataFrame, dict[tuple[str, str, float], dict[str, object]]]:
    cap_fits: dict[tuple[str, str, float], dict[str, object]] = {}
    rows = []
    lcdm_high = conditional_score_from_fit(
        partitions, "LCDM", lcdm_low, h0
    )
    for model_key in ("PATH_FREE", "GENERAL_FREE"):
        model_profile = profile[profile["model"] == model_key]
        for sample_key, sample_blocks, chi2_column, p_column in (
            ("full", blocks, "full_chi2", "full_p"),
            ("low_z", low_blocks, "low_z_chi2", "low_z_p"),
        ):
            starts = model_profile[["gamma_c", chi2_column, p_column]].rename(
                columns={chi2_column: "chi2", p_column: "p"}
            )
            previous: dict[str, float] | None = None
            for cap in GAMMA_CAPS:
                fit = fit_cap(
                    sample_blocks,
                    model_key,
                    cap,
                    h0,
                    grid_starts=starts,
                    extra_start=previous,
                )
                cap_fits[(model_key, sample_key, cap)] = fit
                previous = {
                    "gamma_c": float(fit["parameters"]["gamma_c"]),
                    "p": float(fit["parameters"]["p"]),
                }
                high_chi2 = (
                    conditional_score_from_fit(partitions, model_key, fit, h0)
                    if sample_key == "low_z"
                    else float("nan")
                )
                rows.append(
                    {
                        "model": model_key,
                        "sample": sample_key,
                        "gamma_cap": cap,
                        "N": fit["N"],
                        "chi2": fit["chi2"],
                        "BIC": fit["BIC"],
                        "gamma_c": fit["parameters"]["gamma_c"],
                        "p": fit["parameters"]["p"],
                        "alpha": fit["parameters"]["alpha"],
                        "boundary_fraction_within_cap": fit[
                            "boundary_fraction_within_cap"
                        ],
                        "best_two_chi2_disagreement": fit[
                            "best_two_chi2_disagreement"
                        ],
                        "high_z_conditional_chi2": high_chi2,
                        "high_z_delta_chi2_vs_LCDM": high_chi2 - lcdm_high
                        if np.isfinite(high_chi2)
                        else float("nan"),
                        "optimisation_success": fit["optimisation_success"],
                    }
                )
                print(
                    f"{model_key} {sample_key} cap={cap:g} "
                    f"gamma={fit['parameters']['gamma_c']:.6f} "
                    f"p={fit['parameters']['p']:.6f} chi2={fit['chi2']:.6f}",
                    flush=True,
                )
    lcdm_full_row = {
        "model": "LCDM",
        "sample": "full",
        "gamma_cap": float("nan"),
        "N": lcdm_full["N"],
        "chi2": lcdm_full["chi2"],
        "BIC": lcdm_full["BIC"],
        "gamma_c": float("nan"),
        "p": float("nan"),
        "alpha": float("nan"),
        "boundary_fraction_within_cap": float("nan"),
        "best_two_chi2_disagreement": float("nan"),
        "high_z_conditional_chi2": float("nan"),
        "high_z_delta_chi2_vs_LCDM": float("nan"),
        "optimisation_success": lcdm_full["optimisation_success"],
    }
    lcdm_low_row = dict(lcdm_full_row)
    lcdm_low_row.update(
        {
            "sample": "low_z",
            "N": lcdm_low["N"],
            "chi2": lcdm_low["chi2"],
            "BIC": lcdm_low["BIC"],
            "high_z_conditional_chi2": lcdm_high,
            "high_z_delta_chi2_vs_LCDM": 0.0,
            "optimisation_success": lcdm_low["optimisation_success"],
        }
    )
    rows.extend([lcdm_full_row, lcdm_low_row])
    return pd.DataFrame(rows), cap_fits


def delta_table(
    grouped: pd.DataFrame,
    candidate_label: str,
    keys: list[str],
) -> pd.DataFrame:
    candidate = grouped[grouped["fit_label"] == candidate_label].copy()
    lcdm = grouped[grouped["fit_label"] == "LCDM_LOWZ"][
        keys + ["chi2_contribution"]
    ].rename(columns={"chi2_contribution": "LCDM_chi2_contribution"})
    merged = candidate.merge(lcdm, on=keys, how="outer")
    merged["chi2_contribution"] = merged["chi2_contribution"].fillna(0.0)
    merged["LCDM_chi2_contribution"] = merged[
        "LCDM_chi2_contribution"
    ].fillna(0.0)
    merged["delta_chi2_candidate_minus_LCDM"] = (
        merged["chi2_contribution"] - merged["LCDM_chi2_contribution"]
    )
    return merged


def build_covariance_sensitivities(
    blocks: list[nested.base.Block],
    low_blocks: list[nested.base.Block],
    partitions: list[WhiteningPartition],
    path_cap8_low: dict[str, object],
    general_cap8_low: dict[str, object],
    lcdm_low: dict[str, object],
    h0: float,
) -> pd.DataFrame:
    model_fits = {
        "PATH_FREE": path_cap8_low,
        "GENERAL_FREE": general_cap8_low,
        "LCDM": lcdm_low,
    }
    rows = []
    for model_key, fit in model_fits.items():
        rows.append(
            {
                "covariance_readout": "released_conditional",
                "model": model_key,
                "high_z_chi2": conditional_score_from_fit(
                    partitions, model_key, fit, h0
                ),
                "training_parameters": serialise(fit["parameters"]),
            }
        )
        rows.append(
            {
                "covariance_readout": "released_marginal",
                "model": model_key,
                "high_z_chi2": marginal_high_score(
                    blocks, model_key, fit, h0
                ),
                "training_parameters": serialise(fit["parameters"]),
            }
        )

    diagonal_originals, diagonal_low = diagonal_blocks(blocks)
    path_diagonal = fit_cap(
        diagonal_low,
        "PATH_FREE",
        8.0,
        h0,
        extra_start=path_cap8_low["parameters"],
    )
    general_diagonal = fit_cap(
        diagonal_low,
        "GENERAL_FREE",
        8.0,
        h0,
        extra_start=general_cap8_low["parameters"],
    )
    lcdm_diagonal = nested.fit_model(
        diagonal_low,
        nested.MODEL_SPECS["LCDM"],
        h0,
        rigorous=True,
    )
    for model_key, fit in (
        ("PATH_FREE", path_diagonal),
        ("GENERAL_FREE", general_diagonal),
        ("LCDM", lcdm_diagonal),
    ):
        rows.append(
            {
                "covariance_readout": "diagonal_refit",
                "model": model_key,
                "high_z_chi2": diagonal_high_score(
                    diagonal_originals, model_key, fit, h0
                ),
                "training_parameters": serialise(fit["parameters"]),
            }
        )
    result = pd.DataFrame(rows)
    lcdm = result[result["model"] == "LCDM"][[
        "covariance_readout",
        "high_z_chi2",
    ]].rename(columns={"high_z_chi2": "LCDM_high_z_chi2"})
    result = result.merge(lcdm, on="covariance_readout", how="left")
    result["delta_chi2_vs_LCDM"] = (
        result["high_z_chi2"] - result["LCDM_high_z_chi2"]
    )
    return result


def build_release_exclusion_sensitivities(
    blocks: list[nested.base.Block],
    path_start: dict[str, object],
    h0: float,
) -> pd.DataFrame:
    rows = []
    exclusions: list[str | None] = [None] + [block.label for block in blocks]
    for excluded in exclusions:
        selected = [block for block in blocks if block.label != excluded]
        low = make_low_blocks(selected)
        partitions = make_whitening_partitions(selected)
        path_fit = fit_cap(
            low,
            "PATH_FREE",
            8.0,
            h0,
            extra_start=path_start["parameters"],
        )
        lcdm_fit = nested.fit_model(
            low,
            nested.MODEL_SPECS["LCDM"],
            h0,
            rigorous=True,
        )
        path_high = conditional_score_from_fit(
            partitions, "PATH_FREE", path_fit, h0
        )
        lcdm_high = conditional_score_from_fit(
            partitions, "LCDM", lcdm_fit, h0
        )
        rows.append(
            {
                "excluded_release": excluded if excluded is not None else "NONE",
                "included_releases": "+".join(block.label for block in selected),
                "N_train": int(sum(block.n for block in low)),
                "N_test": int(sum(partition.n_high for partition in partitions)),
                "PATH_gamma_c": float(path_fit["parameters"]["gamma_c"]),
                "PATH_p": float(path_fit["parameters"]["p"]),
                "PATH_high_z_chi2": path_high,
                "LCDM_Omega_m": float(lcdm_fit["parameters"]["Omega_m"]),
                "LCDM_high_z_chi2": lcdm_high,
                "delta_chi2_PATH_minus_LCDM": path_high - lcdm_high,
                "status": "post_hoc_influence_sensitivity",
            }
        )
    return pd.DataFrame(rows)


def safe_json(value):
    if isinstance(value, dict):
        return {str(key): safe_json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [safe_json(item) for item in value]
    if isinstance(value, (np.floating, float)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    return value


def make_plot(
    profile: pd.DataFrame,
    caps: pd.DataFrame,
    bin_delta: pd.DataFrame,
    output: Path,
) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(13.0, 9.0))
    colours = {"PATH_FREE": "#c92a2a", "GENERAL_FREE": "#1971c2"}
    for model_key in ("PATH_FREE", "GENERAL_FREE"):
        selected = profile[profile["model"] == model_key]
        axes[0, 0].plot(
            selected["gamma_c"],
            np.log10(1.0 + selected["full_chi2"] - selected["full_chi2"].min()),
            label=f"{model_key} full",
            color=colours[model_key],
            linewidth=1.8,
        )
        axes[0, 0].plot(
            selected["gamma_c"],
            np.log10(1.0 + selected["low_z_chi2"] - selected["low_z_chi2"].min()),
            label=f"{model_key} z<0.5",
            color=colours[model_key],
            linestyle="--",
            linewidth=1.4,
        )
        axes[0, 1].plot(
            selected["gamma_c"],
            selected["high_z_delta_chi2_vs_LCDM"],
            label=model_key,
            color=colours[model_key],
            linewidth=1.8,
        )
    for cap in GAMMA_CAPS:
        axes[0, 0].axvline(cap, color="#868e96", linewidth=0.7, alpha=0.6)
    axes[0, 0].axhline(
        np.log10(1.0 + 3.84), color="black", linewidth=0.8, linestyle=":"
    )
    axes[0, 0].set_xlabel("gamma_c")
    axes[0, 0].set_ylabel("log10(1 + profile Delta chi-squared)")
    axes[0, 0].set_title("Boundary profiles")
    axes[0, 0].legend(fontsize=8)
    axes[0, 0].grid(alpha=0.2)

    axes[0, 1].axhline(0.0, color="black", linewidth=0.8)
    axes[0, 1].set_xlabel("gamma_c fitted on z<0.5")
    axes[0, 1].set_ylabel("High-z Delta chi-squared vs Lambda-CDM")
    axes[0, 1].set_title("Conditional transfer")
    axes[0, 1].set_yscale("symlog", linthresh=10.0)
    axes[0, 1].legend(fontsize=8)
    axes[0, 1].grid(alpha=0.2)

    path_caps = caps[caps["model"] == "PATH_FREE"]
    for sample, marker, colour in (
        ("full", "o", "#c92a2a"),
        ("low_z", "s", "#2b8a3e"),
    ):
        chosen = path_caps[path_caps["sample"] == sample]
        axes[1, 0].plot(
            chosen["gamma_cap"],
            chosen["gamma_c"],
            marker=marker,
            color=colour,
            label=sample,
        )
    axes[1, 0].plot([2, 8], [2, 8], color="#868e96", linestyle=":", label="at cap")
    axes[1, 0].set_xlabel("Registered diagnostic gamma cap")
    axes[1, 0].set_ylabel("Fitted gamma_c")
    axes[1, 0].set_title("PATH cap response")
    axes[1, 0].legend(fontsize=8)
    axes[1, 0].grid(alpha=0.2)

    labels = [f"{row.release}\n{row.redshift_bin}" for row in bin_delta.itertuples()]
    values = bin_delta["delta_chi2_candidate_minus_LCDM"].to_numpy()
    bar_colours = ["#2b8a3e" if value < 0.0 else "#c92a2a" for value in values]
    positions = np.arange(len(values))
    axes[1, 1].bar(positions, values, color=bar_colours)
    axes[1, 1].axhline(0.0, color="black", linewidth=0.8)
    axes[1, 1].set_xticks(positions)
    axes[1, 1].set_xticklabels(labels, rotation=65, ha="right", fontsize=7)
    axes[1, 1].set_ylabel("Delta chi-squared vs Lambda-CDM")
    axes[1, 1].set_title("Cap-8 PATH conditional decomposition")
    axes[1, 1].grid(axis="y", alpha=0.2)
    figure.tight_layout()
    figure.savefig(output, dpi=180)
    plt.close(figure)


def write_report(
    output: Path,
    caps: pd.DataFrame,
    profile: pd.DataFrame,
    bin_delta: pd.DataFrame,
    survey_delta: pd.DataFrame,
    object_delta: pd.DataFrame,
    sensitivities: pd.DataFrame,
    release_exclusions: pd.DataFrame,
    gates: dict[str, bool],
    metrics: dict[str, object],
) -> None:
    path_caps = caps[caps["model"] == "PATH_FREE"]
    cap8_full = path_caps[
        (path_caps["sample"] == "full") & (path_caps["gamma_cap"] == 8.0)
    ].iloc[0]
    cap8_low = path_caps[
        (path_caps["sample"] == "low_z") & (path_caps["gamma_cap"] == 8.0)
    ].iloc[0]
    lines = [
        "# PLAMB boundary and high-redshift audit",
        "",
        f"**Analysis date:** {DATE_TAG}  ",
        "**Status:** post hoc diagnostic; cannot promote a PLAMB claim  ",
        "**Locked data:** 3,422 supernovae with released total covariance",
        "",
        "## Headline result",
        "",
        f"Under the widened cap gamma_c<=8, the full-sample PATH fit gives gamma_c="
        f"{cap8_full.gamma_c:.6f}, p={cap8_full.p:.6f} and chi-squared={cap8_full.chi2:.6f}. "
        f"Its boundary fraction is {cap8_full.boundary_fraction_within_cap:.6f}, and the "
        f"fixed-grid profile rise at gamma_c=8 is {metrics['profile_rise_at_gamma8']:.6f}.",
        "",
        f"The widened PATH curve is almost identical to Lambda-CDM in raw chi-squared, but "
        f"still has Delta BIC={metrics['cap8_path_full_delta_BIC_vs_LCDM']:+.6f} because it "
        "uses one additional shape parameter.",
        "",
        f"The corresponding z<0.5 fit gives gamma_c={cap8_low.gamma_c:.6f}, "
        f"p={cap8_low.p:.6f}. Its conditional z>=0.5 prediction has "
        f"Delta chi-squared={metrics['cap8_path_high_z_delta_chi2_vs_LCDM']:+.6f} "
        "relative to Lambda-CDM.",
        "",
        f"A narrow outcome-visible grid interval, gamma_c="
        f"{metrics['transfer_win_gamma_min']:.2f}-{metrics['transfer_win_gamma_max']:.2f}, "
        "would beat Lambda-CDM on the held-out quadratic score. Its best grid point is "
        f"gamma_c={metrics['best_transfer_gamma']:.2f}, but that point is worse than the "
        f"independently selected low-redshift optimum by Delta chi-squared="
        f"{metrics['best_transfer_low_z_penalty']:.6f}. It is therefore an exploratory "
        "tension diagnostic, not a valid retuning of the predictor.",
        "",
        f"**Diagnostic decision: {'PASS' if all(gates.values()) else 'FAIL'}.** "
        "This decision concerns boundedness and transfer only; a pass would still not be "
        "confirmatory evidence because the wider bound was outcome-motivated.",
        "",
        "## Cap fits",
        "",
        "| Model | Sample | gamma cap | gamma_c | p | alpha | chi-squared | high-z Delta chi-squared | Bound fraction |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in caps[caps["model"] != "LCDM"].itertuples():
        high_delta = (
            f"{row.high_z_delta_chi2_vs_LCDM:+.6f}"
            if np.isfinite(row.high_z_delta_chi2_vs_LCDM)
            else "-"
        )
        lines.append(
            f"| `{row.model}` | {row.sample} | {row.gamma_cap:.0f} | {row.gamma_c:.6f} | "
            f"{row.p:.6f} | {row.alpha:.6f} | {row.chi2:.6f} | "
            f"{high_delta} | {row.boundary_fraction_within_cap:.6f} |"
        )
    lines.extend(
        [
            "",
            "## Additive conditional decomposition",
            "",
            "The rows below are exact Cholesky-innovation contributions for the cap-8 PATH "
            "low-redshift fit relative to the independently fitted low-redshift Lambda-CDM control.",
            "",
            "| Release | Redshift band | PATH chi-squared | LCDM chi-squared | Delta chi-squared |",
            "|---|---|---:|---:|---:|",
        ]
    )
    for row in bin_delta.itertuples():
        lines.append(
            f"| {row.release} | `{row.redshift_bin}` | {row.chi2_contribution:.6f} | "
            f"{row.LCDM_chi2_contribution:.6f} | {row.delta_chi2_candidate_minus_LCDM:+.6f} |"
        )
    top_surveys = survey_delta.sort_values(
        "delta_chi2_candidate_minus_LCDM", ascending=False
    ).head(10)
    top_objects = object_delta.sort_values(
        "delta_chi2_candidate_minus_LCDM", ascending=False
    ).head(12)
    lines.extend(
        [
            "",
            "## Largest survey contributions",
            "",
            "| Release | Survey | Delta chi-squared |",
            "|---|---|---:|",
        ]
    )
    for row in top_surveys.itertuples():
        lines.append(
            f"| {row.release} | `{row.survey}` | {row.delta_chi2_candidate_minus_LCDM:+.6f} |"
        )
    lines.extend(
        [
            "",
            "## Largest ordered-innovation contributions",
            "",
            "| Release | Source index | z | Survey | Redshift band | Delta chi-squared |",
            "|---|---:|---:|---|---|---:|",
        ]
    )
    for row in top_objects.itertuples():
        lines.append(
            f"| {row.release} | {row.source_index} | {row.z:.6f} | `{row.survey}` | "
            f"`{row.redshift_bin}` | {row.delta_chi2_candidate_minus_LCDM:+.6f} |"
        )
    lines.extend(
        [
            "",
            "Individual values are ordered Gaussian innovations after sorting by redshift and "
            "source index. Their allocation is order-dependent, while their sum is the exact, "
            "order-invariant conditional quadratic form.",
            "",
            "Union3.1 is a 22-node compressed likelihood. Its source indices are compressed "
            "redshift nodes, not individual supernovae. The dominant z=2.262260 entry therefore "
            "identifies endpoint leverage in the compressed release rather than one aberrant event.",
            "",
            "## Covariance sensitivities",
            "",
            "| Readout | Model | High-z chi-squared | Delta chi-squared vs LCDM |",
            "|---|---|---:|---:|",
        ]
    )
    for row in sensitivities.itertuples():
        lines.append(
            f"| `{row.covariance_readout}` | `{row.model}` | {row.high_z_chi2:.6f} | "
            f"{row.delta_chi2_vs_LCDM:+.6f} |"
        )
    lines.extend(
        [
            "",
            "## Exploratory release-exclusion sensitivity",
            "",
            "This refit was triggered by the dominant Union3.1 compressed endpoint and was not "
            "a registered gate.",
            "",
            "| Excluded release | PATH gamma_c | PATH p | N test | Delta high-z chi-squared |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for row in release_exclusions.itertuples():
        lines.append(
            f"| `{row.excluded_release}` | {row.PATH_gamma_c:.6f} | {row.PATH_p:.6f} | "
            f"{row.N_test} | {row.delta_chi2_PATH_minus_LCDM:+.6f} |"
        )
    lines.extend(
        [
            "",
            "The conditional score is primary. Marginal and diagonal rows diagnose whether "
            "released cross-covariance or individual variances control the sign.",
            "The general low-redshift control reaches alpha=-0.5, its registered lower bound, "
            "so it does not provide an identified alternative explanation of the transfer failure.",
            "",
            "## Registered diagnostic gates",
            "",
            "| Gate | Pass |",
            "|---|---|",
        ]
    )
    for name, passed in gates.items():
        lines.append(f"| `{name}` | **{bool(passed)}** |")
    lines.extend(
        [
            "",
            "## Interpretation boundary",
            "",
            "This audit can establish whether the earlier gamma boundary was artificial and "
            "whether the associated curve transfers beyond z=0.5. It cannot convert a post hoc "
            "bound extension into evidence for clock-rate or light-speed evolution. A physical "
            "claim still requires a pre-data derivation of gamma_c, p and the flux exponent, "
            "followed by a fresh external-data test.",
            "",
            "## Reproduction",
            "",
            "```powershell",
            "python github_export/code/rawmu/audit_plamb_boundary_highz_2026_07_18.py",
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
        outdir / f"plamb_boundary_highz_manifest_{DATE_TAG}.csv", index=False
    )


def run_self_tests() -> None:
    assert len(registered_gamma_grid()) == 77
    assert redshift_bin(0.5) == "z0p50_0p65"
    assert redshift_bin(0.65) == "z0p65_0p80"
    assert redshift_bin(0.8) == "z0p80_1p00"
    assert redshift_bin(1.0) == "z1p00_plus"

    z = np.linspace(0.05, 1.1, 28)
    gamma_c, p = 1.2, 0.35
    mu = nested.clock_mu_alpha0(z, gamma_c, p, 67.5) + 0.17
    correlation = 0.15 ** np.abs(np.subtract.outer(np.arange(len(z)), np.arange(len(z))))
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
    fit = profile_p_at_gamma(low, "PATH_FREE", gamma_c, 67.5)
    assert abs(float(fit["parameters"]["p"]) - p) < 1e-6
    partitions = make_whitening_partitions([block])
    innovations = innovations_for_fit(
        partitions, "PATH_FREE", fit, 67.5, "synthetic"
    )
    direct = conditional_score_from_fit(
        partitions, "PATH_FREE", fit, 67.5
    )
    assert np.isclose(innovations["chi2_contribution"].sum(), direct, rtol=1e-12)
    assert direct < 1e-12
    print("All boundary/high-redshift self-tests passed.")


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
    if protocol["status_at_registration"] != "no_gamma_c_above_2_outcome_run_or_inspected":
        raise ValueError("Protocol status does not preserve the diagnostic boundary")
    cli.outdir.mkdir(parents=True, exist_ok=True)

    variants, metadata = nested.base.load_blocks(cli.min_survey_n)
    blocks = variants["released_total_primary"]
    if sum(block.n for block in blocks) != protocol["sample"]["expected_N"]:
        raise ValueError("Locked sample size changed")
    low_blocks = make_low_blocks(blocks)
    partitions = make_whitening_partitions(blocks)

    lcdm_full = nested.fit_model(
        blocks, nested.MODEL_SPECS["LCDM"], cli.h0, rigorous=True
    )
    lcdm_low = nested.fit_model(
        low_blocks, nested.MODEL_SPECS["LCDM"], cli.h0, rigorous=True
    )
    lcdm_high = conditional_score_from_fit(
        partitions, "LCDM", lcdm_low, cli.h0
    )
    print(
        f"LCDM full chi2={lcdm_full['chi2']:.6f}; high-z conditional={lcdm_high:.6f}",
        flush=True,
    )

    profile, _profile_fits = build_profile_scan(
        blocks, low_blocks, partitions, lcdm_high, cli.h0
    )
    profile.to_csv(
        cli.outdir / f"plamb_boundary_highz_profiles_{DATE_TAG}.csv", index=False
    )
    caps, cap_fits = build_cap_fits(
        blocks,
        low_blocks,
        partitions,
        profile,
        lcdm_full,
        lcdm_low,
        cli.h0,
    )
    caps.to_csv(
        cli.outdir / f"plamb_boundary_highz_cap_fits_{DATE_TAG}.csv", index=False
    )

    labelled_fits = {
        "PATH_CAP2_LOWZ": (
            "PATH_FREE",
            cap_fits[("PATH_FREE", "low_z", 2.0)],
        ),
        "PATH_CAP4_LOWZ": (
            "PATH_FREE",
            cap_fits[("PATH_FREE", "low_z", 4.0)],
        ),
        "PATH_CAP8_LOWZ": (
            "PATH_FREE",
            cap_fits[("PATH_FREE", "low_z", 8.0)],
        ),
        "GENERAL_CAP8_LOWZ": (
            "GENERAL_FREE",
            cap_fits[("GENERAL_FREE", "low_z", 8.0)],
        ),
        "LCDM_LOWZ": ("LCDM", lcdm_low),
    }
    innovation_frames = []
    for label, (model_key, fit) in labelled_fits.items():
        innovation_frames.append(
            innovations_for_fit(partitions, model_key, fit, cli.h0, label)
        )
    innovations = pd.concat(innovation_frames, ignore_index=True)
    innovations.to_csv(
        cli.outdir / f"plamb_boundary_highz_ordered_innovations_{DATE_TAG}.csv",
        index=False,
    )

    grouped_bins = (
        innovations.groupby(
            ["fit_label", "model", "release", "redshift_bin"], as_index=False
        )["chi2_contribution"]
        .sum()
        .sort_values(["fit_label", "release", "redshift_bin"])
    )
    grouped_bins.to_csv(
        cli.outdir / f"plamb_boundary_highz_bin_contributions_{DATE_TAG}.csv",
        index=False,
    )
    grouped_surveys = (
        innovations.groupby(
            ["fit_label", "model", "release", "survey"], as_index=False
        )["chi2_contribution"]
        .sum()
        .sort_values(["fit_label", "release", "survey"])
    )
    grouped_surveys.to_csv(
        cli.outdir / f"plamb_boundary_highz_survey_contributions_{DATE_TAG}.csv",
        index=False,
    )

    candidate_label = "PATH_CAP8_LOWZ"
    bin_delta = delta_table(
        grouped_bins, candidate_label, ["release", "redshift_bin"]
    )
    survey_delta = delta_table(
        grouped_surveys, candidate_label, ["release", "survey"]
    )
    candidate_objects = innovations[innovations["fit_label"] == candidate_label].copy()
    lcdm_objects = innovations[innovations["fit_label"] == "LCDM_LOWZ"][
        ["release", "source_index", "chi2_contribution"]
    ].rename(columns={"chi2_contribution": "LCDM_chi2_contribution"})
    object_delta = candidate_objects.merge(
        lcdm_objects, on=["release", "source_index"], how="inner"
    )
    object_delta["delta_chi2_candidate_minus_LCDM"] = (
        object_delta["chi2_contribution"]
        - object_delta["LCDM_chi2_contribution"]
    )
    bin_delta.to_csv(
        cli.outdir / f"plamb_boundary_highz_bin_deltas_{DATE_TAG}.csv", index=False
    )
    survey_delta.to_csv(
        cli.outdir / f"plamb_boundary_highz_survey_deltas_{DATE_TAG}.csv", index=False
    )
    object_delta.to_csv(
        cli.outdir / f"plamb_boundary_highz_object_deltas_{DATE_TAG}.csv", index=False
    )

    for label in labelled_fits:
        direct = float(
            innovations.loc[
                innovations["fit_label"] == label, "chi2_contribution"
            ].sum()
        )
        grouped = float(
            grouped_bins.loc[
                grouped_bins["fit_label"] == label, "chi2_contribution"
            ].sum()
        )
        if not np.isclose(direct, grouped, rtol=1e-8, atol=1e-8):
            raise AssertionError(f"Conditional decomposition failed for {label}")

    sensitivities = build_covariance_sensitivities(
        blocks,
        low_blocks,
        partitions,
        cap_fits[("PATH_FREE", "low_z", 8.0)],
        cap_fits[("GENERAL_FREE", "low_z", 8.0)],
        lcdm_low,
        cli.h0,
    )
    sensitivities.to_csv(
        cli.outdir / f"plamb_boundary_highz_covariance_sensitivities_{DATE_TAG}.csv",
        index=False,
    )
    release_exclusions = build_release_exclusion_sensitivities(
        blocks,
        cap_fits[("PATH_FREE", "low_z", 8.0)],
        cli.h0,
    )
    release_exclusions.to_csv(
        cli.outdir / f"plamb_boundary_highz_release_exclusions_{DATE_TAG}.csv",
        index=False,
    )

    path_cap8_full = cap_fits[("PATH_FREE", "full", 8.0)]
    path_cap8_low = cap_fits[("PATH_FREE", "low_z", 8.0)]
    path_profile = profile[profile["model"] == "PATH_FREE"]
    gamma8_chi2 = float(
        path_profile.loc[np.isclose(path_profile["gamma_c"], 8.0), "full_chi2"].iloc[0]
    )
    profile_rise = gamma8_chi2 - float(path_cap8_full["chi2"])
    candidate_high = float(
        innovations.loc[
            innovations["fit_label"] == candidate_label, "chi2_contribution"
        ].sum()
    )
    lcdm_high_from_innovations = float(
        innovations.loc[
            innovations["fit_label"] == "LCDM_LOWZ", "chi2_contribution"
        ].sum()
    )
    high_delta = candidate_high - lcdm_high_from_innovations
    net_positive = max(high_delta, 0.0)
    max_object_positive = max(
        float(object_delta["delta_chi2_candidate_minus_LCDM"].max()), 0.0
    )
    genuine_surveys = survey_delta[survey_delta["survey"] != "COMPRESSED"]
    max_survey_positive = max(
        float(genuine_surveys["delta_chi2_candidate_minus_LCDM"].max()), 0.0
    )
    object_fraction = max_object_positive / net_positive if net_positive > 0.0 else 0.0
    survey_fraction = max_survey_positive / net_positive if net_positive > 0.0 else 0.0
    path_sensitivities = sensitivities[sensitivities["model"] == "PATH_FREE"]
    sensitivity_signs = np.sign(path_sensitivities["delta_chi2_vs_LCDM"].to_numpy())
    signs_match = bool(np.all(sensitivity_signs == sensitivity_signs[0]))
    cap8_agreement = max(
        float(path_cap8_full["best_two_chi2_disagreement"]),
        float(path_cap8_low["best_two_chi2_disagreement"]),
    )
    transfer_wins = path_profile[
        path_profile["high_z_delta_chi2_vs_LCDM"] <= 0.0
    ]
    best_transfer = path_profile.loc[
        path_profile["high_z_delta_chi2_vs_LCDM"].idxmin()
    ]
    low_z_profile_minimum = float(path_profile["low_z_chi2"].min())
    full_delta_bic = float(path_cap8_full["BIC"] - lcdm_full["BIC"])
    gates = {
        "cap8_full_optimum_at_least_5_percent_from_boundary": bool(
            float(path_cap8_full["boundary_fraction_within_cap"]) >= 0.05
        ),
        "profile_rise_at_gamma8_at_least_3p84": bool(profile_rise >= 3.84),
        "cap8_best_two_optimisers_agree_within_1e_4": bool(
            cap8_agreement <= 1e-4
        ),
        "cap8_PATH_high_z_delta_chi2_not_positive": bool(high_delta <= 0.0),
        "no_release_redshift_bin_loss_above_10": bool(
            float(bin_delta["delta_chi2_candidate_minus_LCDM"].max()) <= 10.0
        ),
        "no_single_innovation_above_25_percent_of_positive_excess": bool(
            object_fraction <= 0.25
        ),
        "no_survey_above_50_percent_of_positive_excess": bool(
            survey_fraction <= 0.50
        ),
        "covariance_sensitivity_signs_match": signs_match,
    }
    metrics = {
        "profile_rise_at_gamma8": profile_rise,
        "cap8_path_full_delta_BIC_vs_LCDM": full_delta_bic,
        "cap8_path_high_z_delta_chi2_vs_LCDM": high_delta,
        "transfer_win_gamma_min": float(transfer_wins["gamma_c"].min())
        if len(transfer_wins)
        else float("nan"),
        "transfer_win_gamma_max": float(transfer_wins["gamma_c"].max())
        if len(transfer_wins)
        else float("nan"),
        "best_transfer_gamma": float(best_transfer["gamma_c"]),
        "best_transfer_delta_chi2_vs_LCDM": float(
            best_transfer["high_z_delta_chi2_vs_LCDM"]
        ),
        "best_transfer_low_z_penalty": float(
            best_transfer["low_z_chi2"] - low_z_profile_minimum
        ),
        "cap8_best_two_chi2_disagreement_max": cap8_agreement,
        "maximum_release_redshift_bin_loss": float(
            bin_delta["delta_chi2_candidate_minus_LCDM"].max()
        ),
        "maximum_single_innovation_positive_excess_fraction": object_fraction,
        "maximum_survey_positive_excess_fraction": survey_fraction,
        "conditional_decomposition_relative_error": float(
            abs(candidate_high - float(grouped_bins.loc[
                grouped_bins["fit_label"] == candidate_label,
                "chi2_contribution",
            ].sum()))
            / max(candidate_high, 1.0)
        ),
    }
    summary = {
        "analysis_date": DATE_TAG,
        "claim_status": protocol["claim_status"],
        "sample_metadata": metadata,
        "gamma_grid_points": len(registered_gamma_grid()),
        "path_cap8_full": path_cap8_full,
        "path_cap8_low_z": path_cap8_low,
        "lcdm_full": lcdm_full,
        "lcdm_low_z": lcdm_low,
        "metrics": metrics,
        "gates": gates,
        "overall_diagnostic_pass": bool(all(gates.values())),
        "release_exclusion_sensitivities": release_exclusions.to_dict(
            orient="records"
        ),
    }
    (cli.outdir / f"plamb_boundary_highz_summary_{DATE_TAG}.json").write_text(
        json.dumps(safe_json(summary), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_report(
        cli.outdir / f"plamb_boundary_highz_report_{DATE_TAG}.md",
        caps,
        profile,
        bin_delta,
        survey_delta,
        object_delta,
        sensitivities,
        release_exclusions,
        gates,
        metrics,
    )
    make_plot(
        profile,
        caps,
        bin_delta,
        cli.outdir / f"plamb_boundary_highz_readout_{DATE_TAG}.png",
    )
    write_manifest(cli.outdir, cli.protocol)
    print(
        f"Saved report: {cli.outdir / f'plamb_boundary_highz_report_{DATE_TAG}.md'}",
        flush=True,
    )


if __name__ == "__main__":
    main()
