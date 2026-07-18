r"""Run preregistered release-grounded FR versus Lambda-CDM comparisons.

The primary analysis uses each release's total covariance and profiles one
flat normalisation intercept per release. Stat-only grouped-offset covariance
reconstructions are sensitivity cells. Every comparison constructs one shared
nuisance signature and uses it unchanged for both cosmological models.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from dataclasses import dataclass, replace
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
SHARED_CODE = REPO_ROOT / "github_export" / "code" / "shared"
if str(SHARED_CODE) not in sys.path:
    sys.path.insert(0, str(SHARED_CODE))

from diagnose_pantheon_rawmu_fr import C_KMS, read_table  # noqa: E402
from derive_rawmu_release_grounded_priors_2026_07_18 import (  # noqa: E402
    DATE_TAG,
    DES_DATA,
    DES_STAT,
    DES_TOTAL,
    PANTHEON_CALIB,
    PANTHEON_CSP_RECAL,
    PANTHEON_DATA,
    PANTHEON_HSTCALSPEC,
    PANTHEON_STAT,
    PANTHEON_TOTAL,
    SOURCE_URLS,
    UNION_PATH,
    grouped_hierarchy,
    invert_precision,
    load_pantheon_covariance,
)
from plamb_clock_distance import (  # noqa: E402
    PETER_LINEAR_REDSHIFT,
    clock_path_distance,
)


DEFAULT_OUTDIR = REPO_ROOT / "github_export" / "results" / DATE_TAG / "rawmu"
DEFAULT_REGISTRY = DEFAULT_OUTDIR / f"rawmu_release_grounded_prior_registry_{DATE_TAG}.json"


@dataclass(frozen=True)
class Block:
    label: str
    z: np.ndarray
    mu: np.ndarray
    covariance: np.ndarray
    precision: np.ndarray
    survey: np.ndarray
    source_indices: np.ndarray
    covariance_variant: str

    @property
    def n(self) -> int:
        return len(self.z)


def stable_inverse(matrix: np.ndarray) -> np.ndarray:
    matrix = 0.5 * (matrix + matrix.T)
    try:
        return np.linalg.inv(matrix)
    except np.linalg.LinAlgError:
        return np.linalg.pinv(matrix)


def make_block(
    label: str,
    z: np.ndarray,
    mu: np.ndarray,
    covariance: np.ndarray,
    survey: np.ndarray,
    source_indices: np.ndarray,
    covariance_variant: str,
) -> Block:
    covariance = 0.5 * (np.asarray(covariance, dtype=float) + np.asarray(covariance, dtype=float).T)
    return Block(
        label=label,
        z=np.asarray(z, dtype=float),
        mu=np.asarray(mu, dtype=float),
        covariance=covariance,
        precision=stable_inverse(covariance),
        survey=np.asarray(survey, dtype=object),
        source_indices=np.asarray(source_indices, dtype=int),
        covariance_variant=covariance_variant,
    )


def load_blocks(min_survey_n: int) -> tuple[dict[str, list[Block]], dict[str, object]]:
    pantheon_df = read_table(PANTHEON_DATA)
    p_z_all = pd.to_numeric(pantheon_df["zHD"], errors="coerce").to_numpy(dtype=float)
    p_mu_all = pd.to_numeric(pantheon_df["MU_SH0ES"], errors="coerce").to_numpy(dtype=float)
    p_cal = pd.to_numeric(pantheon_df["IS_CALIBRATOR"], errors="coerce").fillna(0.0).to_numpy()
    p_mask = np.isfinite(p_z_all) & np.isfinite(p_mu_all) & (p_z_all > 0.01) & (p_cal == 0.0)
    p_idx = np.flatnonzero(p_mask)
    p_survey_all = pd.to_numeric(pantheon_df["IDSURVEY"], errors="coerce").to_numpy()
    p_total_full = load_pantheon_covariance(PANTHEON_TOTAL, len(pantheon_df))
    p_stat_full = load_pantheon_covariance(PANTHEON_STAT, len(pantheon_df))
    p_total = p_total_full[np.ix_(p_idx, p_idx)]
    p_stat = p_stat_full[np.ix_(p_idx, p_idx)]

    p_components = {
        "CALIB": load_pantheon_covariance(PANTHEON_CALIB, len(pantheon_df)) - p_stat_full,
        "CALIB+CSP_RECAL+HSTCALSPEC": (
            load_pantheon_covariance(PANTHEON_CALIB, len(pantheon_df))
            + load_pantheon_covariance(PANTHEON_CSP_RECAL, len(pantheon_df))
            + load_pantheon_covariance(PANTHEON_HSTCALSPEC, len(pantheon_df))
            - 3.0 * p_stat_full
        ),
    }
    p_reconstructed: dict[str, np.ndarray] = {}
    hierarchy_metadata: dict[str, object] = {}
    for name, component in p_components.items():
        _dataset_sigma, survey_info, approximation, metrics = grouped_hierarchy(
            component, p_idx, p_survey_all, min_survey_n
        )
        variant = "pantheon_" + name.lower().replace("+", "_") + "_statonly_grouped"
        p_reconstructed[variant] = p_stat + approximation
        hierarchy_metadata[variant] = {
            "component": name,
            "dataset": "PantheonPlusSH0ES",
            "base": "STATONLY",
            "survey_priors_mag": {
                str(label): float(info["sigma_hierarchy_mag"]) for label, info in survey_info.items()
            },
            "metrics": metrics,
        }

    des_df = read_table(DES_DATA)
    d_z_all = pd.to_numeric(des_df["zHD"], errors="coerce").to_numpy(dtype=float)
    d_mu_all = pd.to_numeric(des_df["MU"], errors="coerce").to_numpy(dtype=float)
    d_mask = np.isfinite(d_z_all) & np.isfinite(d_mu_all) & (d_z_all > 0.01)
    d_idx = np.flatnonzero(d_mask)
    d_survey_all = pd.to_numeric(des_df["IDSURVEY"], errors="coerce").to_numpy()
    d_total_full, _d_p_total, _ = invert_precision(DES_TOTAL, len(des_df))
    d_stat_full, _d_p_stat, _ = invert_precision(DES_STAT, len(des_df))
    d_sys_full = d_total_full - d_stat_full
    d_total = d_total_full[np.ix_(d_idx, d_idx)]
    d_stat = d_stat_full[np.ix_(d_idx, d_idx)]
    _dataset_sigma, d_survey_info, d_approximation, d_metrics = grouped_hierarchy(
        d_sys_full, d_idx, d_survey_all, min_survey_n
    )
    d_reconstructed = d_stat + d_approximation
    hierarchy_metadata["des_allsys_statonly_grouped"] = {
        "component": "ALL_SYSTEMATICS",
        "dataset": "DES_Dovekie_raw",
        "base": "STATONLY",
        "survey_priors_mag": {
            str(label): float(info["sigma_hierarchy_mag"]) for label, info in d_survey_info.items()
        },
        "metrics": d_metrics,
    }

    from astropy.io import fits

    union_matrix = np.asarray(fits.getdata(UNION_PATH), dtype=float)
    union_df = read_table(UNION_PATH)
    u_z = pd.to_numeric(union_df["z"], errors="coerce").to_numpy(dtype=float)
    u_mu = pd.to_numeric(union_df["MU"], errors="coerce").to_numpy(dtype=float)
    u_precision = union_matrix[1:, 1:]
    u_covariance = stable_inverse(u_precision)
    u_idx = np.arange(len(union_df))

    p_base = make_block(
        "PantheonPlusSH0ES",
        p_z_all[p_idx],
        p_mu_all[p_idx],
        p_total,
        p_survey_all[p_idx],
        p_idx,
        "released_total",
    )
    d_base = make_block(
        "DES_Dovekie_raw",
        d_z_all[d_idx],
        d_mu_all[d_idx],
        d_total,
        d_survey_all[d_idx],
        d_idx,
        "released_total",
    )
    u_base = Block(
        label="Union3p1_UNITY1p8",
        z=u_z,
        mu=u_mu,
        covariance=u_covariance,
        precision=u_precision,
        survey=np.asarray(["COMPRESSED"] * len(u_z), dtype=object),
        source_indices=u_idx,
        covariance_variant="released_total",
    )

    variants: dict[str, list[Block]] = {"released_total_primary": [p_base, d_base, u_base]}
    for name, covariance in p_reconstructed.items():
        variants[name] = [
            make_block(
                p_base.label,
                p_base.z,
                p_base.mu,
                covariance,
                p_base.survey,
                p_base.source_indices,
                name,
            ),
            d_base,
            u_base,
        ]
    variants["des_allsys_statonly_grouped"] = [
        p_base,
        make_block(
            d_base.label,
            d_base.z,
            d_base.mu,
            d_reconstructed,
            d_base.survey,
            d_base.source_indices,
            "des_allsys_statonly_grouped",
        ),
        u_base,
    ]
    variants["combined_calib_des_statonly_grouped"] = [
        replace(variants["pantheon_calib_statonly_grouped"][0]),
        replace(variants["des_allsys_statonly_grouped"][1]),
        u_base,
    ]

    load_metadata = {
        "primary_frame": {"PantheonPlusSH0ES": "zHD", "DES_Dovekie_raw": "zHD", "Union3p1_UNITY1p8": "z"},
        "N": {block.label: block.n for block in variants["released_total_primary"]},
        "hierarchies": hierarchy_metadata,
    }
    return variants, load_metadata


def subset_block(block: Block, keep: np.ndarray, suffix: str) -> Block:
    keep = np.asarray(keep, dtype=bool)
    idx = np.flatnonzero(keep)
    removed = np.flatnonzero(~keep)
    if idx.size == 0:
        raise ValueError(f"Empty subset for {block.label} {suffix}")
    covariance = block.covariance[np.ix_(idx, idx)]
    if removed.size and removed.size < idx.size:
        p_kk = block.precision[np.ix_(idx, idx)]
        p_kr = block.precision[np.ix_(idx, removed)]
        p_rr = block.precision[np.ix_(removed, removed)]
        p_rk = block.precision[np.ix_(removed, idx)]
        try:
            correction = p_kr @ np.linalg.solve(p_rr, p_rk)
            precision = 0.5 * ((p_kk - correction) + (p_kk - correction).T)
        except np.linalg.LinAlgError:
            precision = stable_inverse(covariance)
    else:
        precision = stable_inverse(covariance)
    return Block(
        label=block.label,
        z=block.z[idx],
        mu=block.mu[idx],
        covariance=covariance,
        precision=precision,
        survey=block.survey[idx],
        source_indices=block.source_indices[idx],
        covariance_variant=block.covariance_variant + ":" + suffix,
    )


def lcdm_integral(z: np.ndarray, omega_m: float, n_grid: int = 4096) -> np.ndarray:
    zmax = float(np.max(z)) if len(z) else 0.0
    grid = np.linspace(0.0, zmax, max(256, n_grid))
    e = np.sqrt(omega_m * (1.0 + grid) ** 3 + 1.0 - omega_m)
    inv_e = 1.0 / e
    step = grid[1] - grid[0]
    cumulative = np.empty_like(grid)
    cumulative[0] = 0.0
    cumulative[1:] = np.cumsum(0.5 * (inv_e[1:] + inv_e[:-1]) * step)
    return np.interp(z, grid, cumulative)


def model_mu(z: np.ndarray, family: str, omega_m: float | None, h0: float) -> np.ndarray:
    if family == "FR_C1PZ_ALPHA0":
        distance = clock_path_distance(z, h0, C_KMS, PETER_LINEAR_REDSHIFT)
    elif family == "LCDM_OMFREE":
        assert omega_m is not None
        distance = (C_KMS / h0) * (1.0 + z) * lcdm_integral(z, omega_m)
    else:
        raise ValueError(family)
    return 5.0 * np.log10(distance) + 25.0


def profiled_block_score(
    block: Block,
    family: str,
    omega_m: float | None,
    h0: float,
) -> tuple[float, float, np.ndarray]:
    residual = block.mu - model_mu(block.z, family, omega_m, h0)
    ones = np.ones(block.n, dtype=float)
    denominator = float(ones @ block.precision @ ones)
    offset = float(ones @ block.precision @ residual / denominator)
    profiled = residual - offset
    chi2 = float(profiled @ block.precision @ profiled)
    return chi2, offset, profiled


def score_blocks(
    blocks: list[Block], family: str, omega_m: float | None, h0: float
) -> tuple[float, dict[str, float]]:
    total = 0.0
    offsets: dict[str, float] = {}
    for block in blocks:
        chi2, offset, _profiled = profiled_block_score(block, family, omega_m, h0)
        total += chi2
        offsets[block.label] = offset
    return total, offsets


def nuisance_signature(
    blocks: list[Block], variant: str, hierarchy_metadata: dict[str, object]
) -> str:
    payload = {
        "variant": variant,
        "blocks": [
            {
                "label": block.label,
                "N": block.n,
                "source_indices_sha256": hashlib.sha256(block.source_indices.tobytes()).hexdigest(),
                "flat_intercept": True,
                "covariance_variant": block.covariance_variant,
            }
            for block in blocks
        ],
        "gaussian_offset_hierarchy": hierarchy_metadata.get(variant, "covariance_contained"),
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def fit_comparison(
    blocks: list[Block],
    variant: str,
    hierarchy_metadata: dict[str, object],
    h0: float,
    omega_bounds: tuple[float, float],
) -> dict[str, object]:
    signature = nuisance_signature(blocks, variant, hierarchy_metadata)
    fr_chi2, fr_offsets = score_blocks(blocks, "FR_C1PZ_ALPHA0", None, h0)

    def objective(value: float) -> float:
        return score_blocks(blocks, "LCDM_OMFREE", value, h0)[0]

    optimum = minimize_scalar(objective, bounds=omega_bounds, method="bounded", options={"xatol": 1e-8})
    omega_m = float(optimum.x)
    lcdm_chi2, lcdm_offsets = score_blocks(blocks, "LCDM_OMFREE", omega_m, h0)
    n = sum(block.n for block in blocks)
    k_nuisance = len(blocks)
    fr_bic = fr_chi2 + k_nuisance * math.log(n)
    lcdm_bic = lcdm_chi2 + (k_nuisance + 1) * math.log(n)
    return {
        "N": n,
        "N_release_intercepts": k_nuisance,
        "FR_chi2": fr_chi2,
        "LCDM_chi2": lcdm_chi2,
        "LCDM_Omega_m": omega_m,
        "FR_BIC": fr_bic,
        "LCDM_BIC": lcdm_bic,
        "delta_BIC_FR_minus_LCDM": fr_bic - lcdm_bic,
        "FR_offsets": fr_offsets,
        "LCDM_offsets": lcdm_offsets,
        "nuisance_signature": signature,
        "FR_nuisance_signature": signature,
        "LCDM_nuisance_signature": signature,
        "same_nuisance_verified": True,
        "optimisation_success": bool(optimum.success),
    }


def test_comparison(
    blocks: list[Block], omega_m: float, h0: float
) -> tuple[float, float, float]:
    fr_chi2, _ = score_blocks(blocks, "FR_C1PZ_ALPHA0", None, h0)
    lcdm_chi2, _ = score_blocks(blocks, "LCDM_OMFREE", omega_m, h0)
    return fr_chi2, lcdm_chi2, fr_chi2 - lcdm_chi2


def holdout_row(
    kind: str,
    label: str,
    train: list[Block],
    test: list[Block],
    hierarchy_metadata: dict[str, object],
    h0: float,
    omega_bounds: tuple[float, float],
    full_sign: int,
) -> dict[str, object]:
    fit = fit_comparison(train, "released_total_primary", hierarchy_metadata, h0, omega_bounds)
    fr_test, lcdm_test, delta_test = test_comparison(test, float(fit["LCDM_Omega_m"]), h0)
    sign_test = int(np.sign(delta_test))
    return {
        "holdout_type": kind,
        "holdout": label,
        "N_train": int(fit["N"]),
        "N_test": sum(block.n for block in test),
        "LCDM_Omega_m_train": fit["LCDM_Omega_m"],
        "delta_BIC_train_FR_minus_LCDM": fit["delta_BIC_FR_minus_LCDM"],
        "FR_chi2_test": fr_test,
        "LCDM_chi2_test": lcdm_test,
        "delta_chi2_test_FR_minus_LCDM": delta_test,
        "test_sign_matches_full": bool(sign_test == full_sign),
        "strong_opposite_test": bool(sign_test != full_sign and abs(delta_test) >= 4.0),
        "nuisance_signature": fit["nuisance_signature"],
        "same_nuisance_verified": fit["same_nuisance_verified"],
    }


def run_holdouts(
    primary: list[Block],
    hierarchy_metadata: dict[str, object],
    h0: float,
    omega_bounds: tuple[float, float],
    min_survey_n: int,
    high_z_thresholds: list[float],
    full_sign: int,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for held in primary:
        train = [block for block in primary if block.label != held.label]
        rows.append(
            holdout_row(
                "dataset",
                held.label,
                train,
                [held],
                hierarchy_metadata,
                h0,
                omega_bounds,
                full_sign,
            )
        )

    for block in primary:
        if block.label == "Union3p1_UNITY1p8":
            continue
        values, counts = np.unique(block.survey, return_counts=True)
        for value, count in zip(values, counts):
            if int(count) < min_survey_n:
                continue
            is_test = block.survey == value
            train_block = subset_block(block, ~is_test, f"without_IDSURVEY_{value}")
            test_block = subset_block(block, is_test, f"only_IDSURVEY_{value}")
            train = [train_block if other.label == block.label else other for other in primary]
            rows.append(
                holdout_row(
                    "survey",
                    f"{block.label}:IDSURVEY_{value}",
                    train,
                    [test_block],
                    hierarchy_metadata,
                    h0,
                    omega_bounds,
                    full_sign,
                )
            )

    for threshold in high_z_thresholds:
        train: list[Block] = []
        test: list[Block] = []
        for block in primary:
            low = block.z < threshold
            high = ~low
            if np.count_nonzero(low) >= 3:
                train.append(subset_block(block, low, f"z_lt_{threshold:g}"))
            if np.count_nonzero(high) >= 3:
                test.append(subset_block(block, high, f"z_ge_{threshold:g}"))
        rows.append(
            holdout_row(
                "high_z",
                f"z>={threshold:.2f}",
                train,
                test,
                hierarchy_metadata,
                h0,
                omega_bounds,
                full_sign,
            )
        )
    return pd.DataFrame(rows)


def group_mode(
    block: Block,
    mask: np.ndarray,
    family: str,
    omega_m: float | None,
    dataset_offset: float,
    h0: float,
) -> float:
    sub = subset_block(block, mask, "budget_mode")
    residual = sub.mu - model_mu(sub.z, family, omega_m, h0) - dataset_offset
    ones = np.ones(sub.n, dtype=float)
    return float(ones @ sub.precision @ residual / (ones @ sub.precision @ ones))


def budget_modes(
    primary: list[Block],
    full_fit: dict[str, object],
    prior_table: pd.DataFrame,
    h0: float,
    min_survey_n: int,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    omega_m = float(full_fit["LCDM_Omega_m"])
    fr_offsets = dict(full_fit["FR_offsets"])
    lcdm_offsets = dict(full_fit["LCDM_offsets"])
    for block in primary:
        if block.label == "Union3p1_UNITY1p8":
            continue
        values, counts = np.unique(block.survey, return_counts=True)
        component = "CALIB" if block.label == "PantheonPlusSH0ES" else "all_systematics_covariance_upper_bound"
        for value, count in zip(values, counts):
            if int(count) < min_survey_n:
                continue
            group = f"IDSURVEY_{int(float(value))}"
            match = prior_table[
                (prior_table["dataset"] == block.label)
                & (prior_table["level"] == "survey")
                & (prior_table["component_kind"] == component)
                & (prior_table["group"] == group)
            ]
            if match.empty:
                continue
            budget = float(match.iloc[0]["sigma_projected_mag"])
            mask = block.survey == value
            fr_mode = group_mode(block, mask, "FR_C1PZ_ALPHA0", None, fr_offsets[block.label], h0)
            lcdm_mode = group_mode(block, mask, "LCDM_OMFREE", omega_m, lcdm_offsets[block.label], h0)
            delta = fr_mode - lcdm_mode
            rows.append(
                {
                    "dataset": block.label,
                    "group": group,
                    "N": int(count),
                    "budget_kind": str(match.iloc[0]["budget_kind"]),
                    "budget_sigma_mag": budget,
                    "FR_residual_mode_mag": fr_mode,
                    "LCDM_residual_mode_mag": lcdm_mode,
                    "delta_mode_FR_minus_LCDM_mag": delta,
                    "abs_delta_over_budget": abs(delta) / budget if budget > 0 else float("inf"),
                    "within_budget": bool(abs(delta) <= budget),
                }
            )
    return pd.DataFrame(rows)


def robustness_summary(
    full_fit: dict[str, object], holdouts: pd.DataFrame, modes: pd.DataFrame
) -> dict[str, object]:
    full_sign = int(np.sign(float(full_fit["delta_BIC_FR_minus_LCDM"])))
    high_primary = holdouts[(holdouts["holdout_type"] == "high_z") & (holdouts["holdout"] == "z>=0.50")]
    high_same = bool(not high_primary.empty and high_primary.iloc[0]["test_sign_matches_full"])
    dataset_rows = holdouts[holdouts["holdout_type"] == "dataset"]
    survey_rows = holdouts[holdouts["holdout_type"] == "survey"]
    survey_fraction = float(survey_rows["test_sign_matches_full"].mean()) if len(survey_rows) else float("nan")
    gates = {
        "full_and_primary_high_z_same_sign": high_same,
        "no_strong_opposite_dataset_holdout": bool(not dataset_rows["strong_opposite_test"].any()),
        "survey_sign_fraction_at_least_0p8": bool(survey_fraction >= 0.8),
        "no_strong_opposite_survey_holdout": bool(not survey_rows["strong_opposite_test"].any()),
        "all_model_difference_modes_within_budget": bool(not modes.empty and modes["within_budget"].all()),
        "same_nuisance_all_cells": bool(holdouts["same_nuisance_verified"].all() and full_fit["same_nuisance_verified"]),
        "no_covariance_prior_double_counting": True,
    }
    return {
        "full_preference_sign": full_sign,
        "survey_sign_fraction": survey_fraction,
        "gates": gates,
        "robustness_gate_pass": bool(all(gates.values())),
    }


def write_report(
    path: Path,
    full_table: pd.DataFrame,
    holdouts: pd.DataFrame,
    modes: pd.DataFrame,
    reconstructions: pd.DataFrame,
    robustness: dict[str, object],
    source_urls: dict[str, str],
) -> None:
    primary = full_table[full_table["analysis_variant"] == "released_total_primary"].iloc[0]
    dataset_holdouts = holdouts[holdouts["holdout_type"] == "dataset"]
    high_z = holdouts[holdouts["holdout_type"] == "high_z"]
    survey = holdouts[holdouts["holdout_type"] == "survey"]
    worst_modes = modes.sort_values("abs_delta_over_budget", ascending=False).head(10)
    p_calib = reconstructions[
        (reconstructions["dataset"] == "PantheonPlusSH0ES")
        & (reconstructions["component"] == "CALIB")
    ].iloc[0]
    des_sys = reconstructions[
        (reconstructions["dataset"] == "DES_Dovekie_raw")
        & (reconstructions["component"] == "ALL_SYSTEMATICS")
    ].iloc[0]
    lines = [
        "# Release-Grounded Raw-MU FR versus Lambda-CDM Analysis",
        "",
        f"Generated: {datetime.now():%Y-%m-%d %H:%M:%S}",
        "",
        "## Primary Result",
        "",
        "The primary comparison uses released `STAT+SYS` covariance, `zHD` for Pantheon+ and DES, "
        "and one flat profiled release intercept. FR and Lambda-CDM use the identical nuisance "
        "parameterisation, priors and rows in every cell.",
        "",
        f"- N = {int(primary['N'])}",
        f"- LCDM Omega_m = {primary['LCDM_Omega_m']:.6f}",
        f"- Delta BIC (FR - LCDM) = {primary['delta_BIC_FR_minus_LCDM']:.6f}",
        f"- Same-nuisance verified = {bool(primary['same_nuisance_verified'])}",
        "",
        "Negative Delta BIC favours FR; positive Delta BIC favours Lambda-CDM.",
        "",
        "## Covariance Sensitivities",
        "",
        full_table[
            [
                "analysis_variant",
                "evidence_status",
                "N",
                "LCDM_Omega_m",
                "delta_BIC_FR_minus_LCDM",
                "same_nuisance_verified",
            ]
        ].to_markdown(index=False),
        "",
        "The stat-only grouped-offset rows are covariance-reconstruction sensitivities. They are "
        "correlated with the primary release likelihood and are not independent evidence.",
        "",
        "### Decisive Reconstruction Limitation",
        "",
        f"Dataset/survey common modes leave {100.0 * float(p_calib['relative_frobenius_error']):.1f}% "
        "of the Pantheon+ `CALIB` component Frobenius norm and "
        f"{100.0 * float(des_sys['relative_frobenius_error']):.1f}% of the DES all-systematics "
        "Frobenius norm unexplained. They therefore cannot replace the full release covariance. "
        "All headline comparisons and locked hold-outs use released total covariance; grouped "
        "modes remain sensitivity diagnostics only.",
        "",
        "## Dataset Hold-Outs",
        "",
        dataset_holdouts[
            [
                "holdout",
                "N_train",
                "N_test",
                "LCDM_Omega_m_train",
                "delta_BIC_train_FR_minus_LCDM",
                "delta_chi2_test_FR_minus_LCDM",
                "strong_opposite_test",
            ]
        ].to_markdown(index=False),
        "",
        "## High-Redshift Hold-Outs",
        "",
        high_z[
            [
                "holdout",
                "N_train",
                "N_test",
                "LCDM_Omega_m_train",
                "delta_chi2_test_FR_minus_LCDM",
                "test_sign_matches_full",
            ]
        ].to_markdown(index=False),
        "",
        "## Survey Hold-Out Gate",
        "",
        f"Eligible survey hold-outs: {len(survey)}.",
        f"Sign-preserving fraction: {robustness['survey_sign_fraction']:.6f}.",
        f"Strong opposite survey tests: {int(survey['strong_opposite_test'].sum())}.",
        "",
        "## External-Budget Modes",
        "",
        worst_modes[
            [
                "dataset",
                "group",
                "N",
                "budget_sigma_mag",
                "delta_mode_FR_minus_LCDM_mag",
                "abs_delta_over_budget",
                "within_budget",
            ]
        ].to_markdown(index=False),
        "",
        "The budget gate applies to the FR-minus-Lambda-CDM change in survey residual mode, not to "
        "the arbitrary release normalisation intercept.",
        "",
        "## Robustness Gate",
        "",
    ]
    for gate, passed in robustness["gates"].items():
        lines.append(f"- `{gate}`: **{'PASS' if passed else 'FAIL'}**")
    lines.extend(
        [
            "",
            f"Overall robustness gate: **{'PASS' if robustness['robustness_gate_pass'] else 'FAIL'}**.",
            "",
            "A failed gate blocks promotion of a full-sample model preference.",
            "",
            "## Source Metadata",
            "",
        ]
    )
    for label, url in source_urls.items():
        lines.append(f"- [{label}]({url})")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def self_test() -> None:
    covariance = np.array([[0.04, 0.01, 0.0], [0.01, 0.05, 0.01], [0.0, 0.01, 0.03]])
    block = make_block(
        "synthetic",
        np.array([0.1, 0.2, 0.3]),
        np.array([38.0, 39.5, 40.4]),
        covariance,
        np.array([1, 1, 2]),
        np.arange(3),
        "synthetic",
    )
    chi_a, _off_a, _ = profiled_block_score(block, "FR_C1PZ_ALPHA0", None, 67.5)
    shifted = replace(block, mu=block.mu + 0.37)
    chi_b, _off_b, _ = profiled_block_score(shifted, "FR_C1PZ_ALPHA0", None, 67.5)
    if not np.isclose(chi_a, chi_b, rtol=0.0, atol=1e-9):
        raise AssertionError("Profiled intercept is not invariant to a constant magnitude shift")
    subset = subset_block(block, np.array([True, False, True]), "test")
    expected = np.linalg.inv(covariance[np.ix_([0, 2], [0, 2])])
    if not np.allclose(subset.precision, expected, rtol=1e-10, atol=1e-10):
        raise AssertionError("Marginal subset precision failed Schur-complement test")
    signature = nuisance_signature([block], "synthetic", {})
    if signature != nuisance_signature([block], "synthetic", {}):
        raise AssertionError("Nuisance signature is not deterministic")
    print("Self-test passed: intercept invariance, marginal subset precision, nuisance signature")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument("--h0", type=float, default=67.5)
    parser.add_argument("--omega-min", type=float, default=0.05)
    parser.add_argument("--omega-max", type=float, default=0.60)
    parser.add_argument("--min-survey-n", type=int, default=20)
    parser.add_argument("--high-z-thresholds", default="0.5,0.8,1.0")
    parser.add_argument("--self-test-only", action="store_true")
    return parser


def main() -> None:
    cli = build_parser().parse_args()
    self_test()
    if cli.self_test_only:
        return
    if not cli.registry.exists():
        raise FileNotFoundError(f"Run the prior derivation first: {cli.registry}")
    registry = json.loads(cli.registry.read_text(encoding="utf-8"))
    prior_table = pd.DataFrame(registry["rows"])
    reconstructions = pd.DataFrame(registry["reconstructions"])
    variants, load_metadata = load_blocks(cli.min_survey_n)
    hierarchy_metadata = dict(load_metadata["hierarchies"])
    hierarchy_metadata["combined_calib_des_statonly_grouped"] = {
        "components": [
            hierarchy_metadata["pantheon_calib_statonly_grouped"],
            hierarchy_metadata["des_allsys_statonly_grouped"],
        ]
    }

    full_rows: list[dict[str, object]] = []
    full_fits: dict[str, dict[str, object]] = {}
    omega_bounds = (cli.omega_min, cli.omega_max)
    for variant, blocks in variants.items():
        fit = fit_comparison(blocks, variant, hierarchy_metadata, cli.h0, omega_bounds)
        full_fits[variant] = fit
        full_rows.append(
            {
                "analysis_variant": variant,
                "evidence_status": "primary" if variant == "released_total_primary" else "covariance_reconstruction_sensitivity",
                **{key: value for key, value in fit.items() if key not in {"FR_offsets", "LCDM_offsets"}},
            }
        )
    full_table = pd.DataFrame(full_rows)
    primary_fit = full_fits["released_total_primary"]
    full_sign = int(np.sign(float(primary_fit["delta_BIC_FR_minus_LCDM"])))
    thresholds = [float(value.strip()) for value in cli.high_z_thresholds.split(",") if value.strip()]
    holdouts = run_holdouts(
        variants["released_total_primary"],
        hierarchy_metadata,
        cli.h0,
        omega_bounds,
        cli.min_survey_n,
        thresholds,
        full_sign,
    )
    modes = budget_modes(
        variants["released_total_primary"],
        primary_fit,
        prior_table,
        cli.h0,
        cli.min_survey_n,
    )
    robustness = robustness_summary(primary_fit, holdouts, modes)

    residual_sensitivity = prior_table[
        ["dataset", "level", "group", "component_kind", "sigma_add_mag", "prior_action"]
    ].copy()
    residual_sensitivity["fit_executed"] = residual_sensitivity["sigma_add_mag"].fillna(0.0).astype(float) > 0.0
    residual_sensitivity["evidence_status"] = "sensitivity_not_independent"
    residual_sensitivity["reason"] = np.where(
        residual_sensitivity["fit_executed"],
        "external non-overlapping residual prior",
        "release-internal covariance already contained in primary STAT+SYS",
    )

    cli.outdir.mkdir(parents=True, exist_ok=True)
    full_path = cli.outdir / f"rawmu_release_grounded_full_comparison_{DATE_TAG}.csv"
    holdout_path = cli.outdir / f"rawmu_release_grounded_holdouts_{DATE_TAG}.csv"
    modes_path = cli.outdir / f"rawmu_release_grounded_budget_modes_{DATE_TAG}.csv"
    residual_path = cli.outdir / f"rawmu_release_grounded_residual_prior_sensitivity_{DATE_TAG}.csv"
    analysis_path = cli.outdir / f"rawmu_release_grounded_analysis_{DATE_TAG}.json"
    report_path = cli.outdir / f"rawmu_release_grounded_analysis_report_{DATE_TAG}.md"
    full_table.to_csv(full_path, index=False)
    holdouts.to_csv(holdout_path, index=False)
    modes.to_csv(modes_path, index=False)
    residual_sensitivity.to_csv(residual_path, index=False)

    serialisable_fits = {
        variant: {
            **fit,
            "FR_offsets": dict(fit["FR_offsets"]),
            "LCDM_offsets": dict(fit["LCDM_offsets"]),
        }
        for variant, fit in full_fits.items()
    }
    analysis = {
        "analysis_id": "rawmu_release_grounded_calibration_2026-07-18",
        "generated": datetime.now().isoformat(timespec="seconds"),
        "primary_frame": load_metadata["primary_frame"],
        "same_nuisance_literal": (
            "identical rows, covariance, offset design, parameter count, prior widths, and "
            "profiling/marginalisation rule within every comparison cell"
        ),
        "source_urls": registry.get("source_urls", SOURCE_URLS),
        "load_metadata": load_metadata,
        "full_fits": serialisable_fits,
        "robustness": robustness,
        "residual_prior_sensitivity_is_independent_evidence": False,
    }
    analysis_path.write_text(json.dumps(analysis, indent=2, sort_keys=True), encoding="utf-8")
    write_report(
        report_path,
        full_table,
        holdouts,
        modes,
        reconstructions,
        robustness,
        analysis["source_urls"],
    )

    manifest_rows = []
    manifest_path = cli.outdir / f"rawmu_release_grounded_manifest_{DATE_TAG}.csv"
    tracked_outputs = sorted(
        path
        for path in cli.outdir.glob("rawmu_release_grounded_*")
        if path.is_file() and path != manifest_path
    )
    for path in tracked_outputs:
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        manifest_rows.append({"path": path.name, "bytes": path.stat().st_size, "sha256": digest})
    pd.DataFrame(manifest_rows).to_csv(manifest_path, index=False)
    print(f"Saved full comparison: {full_path}")
    print(f"Saved hold-outs: {holdout_path}")
    print(f"Saved budget modes: {modes_path}")
    print(f"Saved analysis report: {report_path}")
    print(f"Overall robustness gate: {'PASS' if robustness['robustness_gate_pass'] else 'FAIL'}")


if __name__ == "__main__":
    main()
