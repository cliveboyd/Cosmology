#!/usr/bin/env python3
"""Run the frozen Peter clock-law v1 physical and raw-MU consistency audit."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Callable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[3]
RAWMU_CODE = SCRIPT_PATH.parent
SHARED_CODE = REPO_ROOT / "github_export" / "code" / "shared"
for code_path in (REPO_ROOT, RAWMU_CODE, SHARED_CODE):
    if str(code_path) not in sys.path:
        sys.path.insert(0, str(code_path))

import run_rawmu_release_grounded_holdouts_2026_07_18 as grounded  # noqa: E402
from plamb_clock_distance import (  # noqa: E402
    PETER_LINEAR_REDSHIFT,
    clock_path_distance,
)


DATE_TAG = "2026-07-19"
ANALYSIS_ID = "peter_clock_law_v1_2026-07-19"
H0_FIXED = 67.5
C_KMS = grounded.C_KMS
S_BOUNDS = (-3.0, 3.0)
OMEGA_BOUNDS = (0.05, 0.60)
KNOWN_EXACT_DELTA_BIC = 94.34498313934
EXPECTED_PROTOCOL_SHA256 = (
    "3f595a34f908771475fd631ba593917dd8f5fefc76dee7d570967fa6ea0fdf6c"
)
DEFAULT_OUTDIR = (
    REPO_ROOT / "github_export" / "results" / DATE_TAG / "peter_clock_law_v1"
)
DEFAULT_PROTOCOL = (
    DEFAULT_OUTDIR / f"peter_clock_law_v1_protocol_{DATE_TAG}.json"
)
DEFAULT_REGISTRY = (
    REPO_ROOT
    / "github_export"
    / "results"
    / "2026-07-18"
    / "rawmu"
    / "rawmu_release_grounded_prior_registry_2026-07-18.json"
)


@dataclass(frozen=True)
class FixedBranch:
    key: str
    label: str
    tau: float
    ell: float

    @property
    def s(self) -> float:
        return self.tau - 0.5 * self.ell


FIXED_BRANCHES = (
    FixedBranch(
        "PETER_NOLOSS_LCONST",
        "Peter no-loss, constant standardised luminosity",
        0.0,
        0.0,
    ),
    FixedBranch(
        "PETER_NOLOSS_PARTICLECOUNT",
        "Peter no-loss, illustrative particle-count luminosity",
        0.0,
        1.0,
    ),
    FixedBranch(
        "PETER_NOLOSS_CHANDRA",
        "Peter no-loss, fixed-constant Chandrasekhar luminosity proxy",
        0.0,
        3.5,
    ),
    FixedBranch(
        "STATIC_DUALITY_LCONST",
        "Static distance-duality, constant standardised luminosity",
        2.0,
        0.0,
    ),
    FixedBranch(
        "STATIC_DUALITY_CHANDRA",
        "Static distance-duality, fixed-constant Chandrasekhar luminosity proxy",
        2.0,
        3.5,
    ),
)
BRANCH_BY_KEY = {branch.key: branch for branch in FIXED_BRANCHES}
MODEL_KEYS = tuple(BRANCH_BY_KEY) + ("COMPOSITE_S_FREE", "FLAT_LCDM")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def safe_json(value):
    if isinstance(value, dict):
        return {str(key): safe_json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [safe_json(item) for item in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        number = float(value)
        return number if math.isfinite(number) else None
    if isinstance(value, np.ndarray):
        return [safe_json(item) for item in value.tolist()]
    return value


def verify_protocol(path: Path) -> dict[str, object]:
    measured = sha256_file(path)
    if measured != EXPECTED_PROTOCOL_SHA256:
        raise RuntimeError(
            "Frozen protocol hash changed before execution: "
            f"expected {EXPECTED_PROTOCOL_SHA256}, found {measured}"
        )
    protocol = json.loads(path.read_text(encoding="utf-8"))
    if protocol.get("analysis_id") != ANALYSIS_ID:
        raise RuntimeError("Protocol analysis_id does not match the runner")
    return protocol


def peter_path_distance(z: np.ndarray) -> np.ndarray:
    values = np.asarray(z, dtype=float)
    return np.asarray(
        clock_path_distance(values, H0_FIXED, C_KMS, PETER_LINEAR_REDSHIFT),
        dtype=float,
    )


def peter_effective_mu(z: np.ndarray, s: float) -> np.ndarray:
    values = np.asarray(z, dtype=float)
    distance = peter_path_distance(values) * np.power(1.0 + values, float(s))
    return 5.0 * np.log10(distance) + 25.0


def lcdm_mu(z: np.ndarray, omega_m: float) -> np.ndarray:
    values = np.asarray(z, dtype=float)
    distance = (
        (C_KMS / H0_FIXED)
        * (1.0 + values)
        * grounded.lcdm_integral(values, float(omega_m))
    )
    return 5.0 * np.log10(distance) + 25.0


def model_mu(model_key: str, z: np.ndarray, parameter: float) -> np.ndarray:
    if model_key in BRANCH_BY_KEY:
        return peter_effective_mu(z, BRANCH_BY_KEY[model_key].s)
    if model_key == "COMPOSITE_S_FREE":
        return peter_effective_mu(z, parameter)
    if model_key == "FLAT_LCDM":
        return lcdm_mu(z, parameter)
    raise ValueError(model_key)


def profiled_block_score(
    block: grounded.Block,
    prediction: np.ndarray,
) -> tuple[float, float, np.ndarray]:
    residual = block.mu - np.asarray(prediction, dtype=float)
    ones = np.ones(block.n, dtype=float)
    denominator = float(ones @ block.precision @ ones)
    offset = float(ones @ block.precision @ residual / denominator)
    profiled = residual - offset
    chi2 = float(profiled @ block.precision @ profiled)
    return chi2, offset, profiled


def score_model(
    blocks: list[grounded.Block],
    model_key: str,
    parameter: float,
) -> tuple[float, dict[str, float]]:
    total = 0.0
    offsets: dict[str, float] = {}
    for block in blocks:
        chi2, offset, _ = profiled_block_score(
            block, model_mu(model_key, block.z, parameter)
        )
        total += chi2
        offsets[block.label] = offset
    return total, offsets


def analytic_composite_fit(
    blocks: list[grounded.Block],
) -> tuple[float, float, float, dict[str, float]]:
    normal = 0.0
    rhs = 0.0
    for block in blocks:
        y = block.mu - peter_effective_mu(block.z, 0.0)
        x = 5.0 * np.log10(1.0 + block.z)
        ones = np.ones(block.n, dtype=float)
        p = block.precision
        one_p_one = float(ones @ p @ ones)
        one_p_x = float(ones @ p @ x)
        one_p_y = float(ones @ p @ y)
        normal += float(x @ p @ x) - one_p_x * one_p_x / one_p_one
        rhs += float(x @ p @ y) - one_p_x * one_p_y / one_p_one
    if normal <= 0.0:
        raise RuntimeError("Composite-power normal equation is not positive")
    unconstrained = rhs / normal
    best = float(np.clip(unconstrained, *S_BOUNDS))
    chi2, offsets = score_model(blocks, "COMPOSITE_S_FREE", best)
    sigma = 1.0 / math.sqrt(normal)
    return best, sigma, chi2, offsets


def bounded_lcdm_fit(
    blocks: list[grounded.Block],
) -> tuple[float, float, dict[str, float], bool]:
    result = minimize_scalar(
        lambda value: score_model(blocks, "FLAT_LCDM", float(value))[0],
        bounds=OMEGA_BOUNDS,
        method="bounded",
        options={"xatol": 1e-10, "maxiter": 600},
    )
    candidates = [
        (float(result.fun), float(result.x)),
        (score_model(blocks, "FLAT_LCDM", OMEGA_BOUNDS[0])[0], OMEGA_BOUNDS[0]),
        (score_model(blocks, "FLAT_LCDM", OMEGA_BOUNDS[1])[0], OMEGA_BOUNDS[1]),
    ]
    chi2, omega_m = min(candidates)
    _, offsets = score_model(blocks, "FLAT_LCDM", omega_m)
    return omega_m, chi2, offsets, bool(result.success or omega_m in OMEGA_BOUNDS)


def fit_models(
    blocks: list[grounded.Block],
    scope: str,
    covariance_variant: str,
) -> tuple[list[dict[str, object]], dict[str, dict[str, object]]]:
    n = int(sum(block.n for block in blocks))
    n_intercepts = len(blocks)
    rows: list[dict[str, object]] = []
    internal: dict[str, dict[str, object]] = {}
    for branch in FIXED_BRANCHES:
        chi2, offsets = score_model(blocks, branch.key, branch.s)
        k_shape = 0
        fit = {
            "parameter_name": "s_fixed",
            "parameter_value": branch.s,
            "parameter_sigma": float("nan"),
            "chi2": chi2,
            "offsets": offsets,
            "success": True,
            "k_shape": k_shape,
        }
        internal[branch.key] = fit
    s_best, s_sigma, s_chi2, s_offsets = analytic_composite_fit(blocks)
    internal["COMPOSITE_S_FREE"] = {
        "parameter_name": "s",
        "parameter_value": s_best,
        "parameter_sigma": s_sigma,
        "chi2": s_chi2,
        "offsets": s_offsets,
        "success": bool(S_BOUNDS[0] < s_best < S_BOUNDS[1]),
        "k_shape": 1,
    }
    omega_m, lcdm_chi2, lcdm_offsets, lcdm_success = bounded_lcdm_fit(blocks)
    internal["FLAT_LCDM"] = {
        "parameter_name": "Omega_m",
        "parameter_value": omega_m,
        "parameter_sigma": float("nan"),
        "chi2": lcdm_chi2,
        "offsets": lcdm_offsets,
        "success": lcdm_success,
        "k_shape": 1,
    }
    lcdm_bic = lcdm_chi2 + (n_intercepts + 1) * math.log(n)
    for model_key in MODEL_KEYS:
        fit = internal[model_key]
        bic = float(fit["chi2"]) + (
            n_intercepts + int(fit["k_shape"])
        ) * math.log(n)
        rows.append(
            {
                "scope": scope,
                "covariance_variant": covariance_variant,
                "model": model_key,
                "N": n,
                "N_release_intercepts": n_intercepts,
                "N_shape_parameters": int(fit["k_shape"]),
                "parameter_name": fit["parameter_name"],
                "parameter_value": float(fit["parameter_value"]),
                "parameter_sigma": float(fit["parameter_sigma"]),
                "chi2": float(fit["chi2"]),
                "BIC": bic,
                "delta_BIC_vs_LCDM": bic - lcdm_bic,
                "offsets": json.dumps(fit["offsets"], sort_keys=True),
                "success": bool(fit["success"]),
            }
        )
        fit["BIC"] = bic
        fit["delta_BIC_vs_LCDM"] = bic - lcdm_bic
    return rows, internal


def subset_blocks(
    blocks: list[grounded.Block],
    predicate: Callable[[np.ndarray], np.ndarray],
    suffix: str,
) -> list[grounded.Block]:
    selected: list[grounded.Block] = []
    for block in blocks:
        keep = np.asarray(predicate(block.z), dtype=bool)
        if np.any(keep):
            selected.append(grounded.subset_block(block, keep, suffix))
    if not selected:
        raise RuntimeError(f"No blocks selected for {suffix}")
    return selected


def theoretical_readout(max_z: float) -> pd.DataFrame:
    grid = np.unique(np.array([0.0, 0.1, 0.5, 1.0, 1.5, max_z], dtype=float))
    one_plus = 1.0 + grid
    c_ratio = one_plus
    mass_ratio = 1.0 / one_plus
    rest_energy_ratio = mass_ratio * c_ratio**2
    q_required = c_ratio / one_plus
    q_claim = 1.0 / one_plus
    return pd.DataFrame(
        {
            "z_observed": grid,
            "C_c_over_c0": c_ratio,
            "R_over_H0": np.ones_like(grid),
            "Q_required_by_spectroscopic_closure": q_required,
            "Q_asserted_clock_control": q_claim,
            "z_predicted_if_asserted_Q_used": c_ratio / q_claim - 1.0,
            "particle_mass_ratio": mass_ratio,
            "rest_energy_ratio_fixed_units": rest_energy_ratio,
            "Compton_frequency_ratio_fixed_h": rest_energy_ratio,
            "h_ratio_needed_for_Q_required": rest_energy_ratio / q_required,
            "h_ratio_needed_for_Q_asserted": rest_energy_ratio / q_claim,
            "M_Ch_ratio_fixed_hbar_G": c_ratio**1.5 / mass_ratio**2,
            "N_Ch_ratio_fixed_hbar_G": c_ratio**1.5 / mass_ratio**3,
            "Xi_needed_for_path_as_DL": one_plus**-2.0,
            "photon_times_luminosity_transfer_needed": one_plus**4.0,
            "dimensionless_path_integral": grid * (1.0 + 0.5 * grid),
        }
    )


def fixed_branch_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "model": branch.key,
                "description": branch.label,
                "tau_photon_distance_power": branch.tau,
                "ell_luminosity_power": branch.ell,
                "s_composite_power": branch.s,
            }
            for branch in FIXED_BRANCHES
        ]
    )


def residual_group_mode(
    block: grounded.Block,
    mask: np.ndarray,
    model_key: str,
    parameter: float,
    release_offset: float,
) -> float:
    sub = grounded.subset_block(block, mask, "peter_v1_budget_mode")
    residual = (
        sub.mu
        - model_mu(model_key, sub.z, parameter)
        - float(release_offset)
    )
    ones = np.ones(sub.n, dtype=float)
    return float(
        ones @ sub.precision @ residual / (ones @ sub.precision @ ones)
    )


def calibration_budget_table(
    blocks: list[grounded.Block],
    fits: dict[str, dict[str, object]],
    registry: dict[str, object],
    min_survey_n: int,
) -> pd.DataFrame:
    prior_table = pd.DataFrame(registry["rows"])
    rows: list[dict[str, object]] = []
    comparisons = (
        ("PETER_NOLOSS_LCONST", "FLAT_LCDM"),
        ("COMPOSITE_S_FREE", "FLAT_LCDM"),
    )
    for block in blocks:
        if block.label == "Union3p1_UNITY1p8":
            continue
        component = (
            "CALIB"
            if block.label == "PantheonPlusSH0ES"
            else "all_systematics_covariance_upper_bound"
        )
        values, counts = np.unique(block.survey, return_counts=True)
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
            for candidate, reference in comparisons:
                candidate_fit = fits[candidate]
                reference_fit = fits[reference]
                candidate_mode = residual_group_mode(
                    block,
                    mask,
                    candidate,
                    float(candidate_fit["parameter_value"]),
                    float(candidate_fit["offsets"][block.label]),
                )
                reference_mode = residual_group_mode(
                    block,
                    mask,
                    reference,
                    float(reference_fit["parameter_value"]),
                    float(reference_fit["offsets"][block.label]),
                )
                delta = candidate_mode - reference_mode
                rows.append(
                    {
                        "candidate": candidate,
                        "reference": reference,
                        "dataset": block.label,
                        "group": group,
                        "N": int(count),
                        "budget_sigma_mag": budget,
                        "candidate_residual_mode_mag": candidate_mode,
                        "reference_residual_mode_mag": reference_mode,
                        "delta_mode_mag": delta,
                        "abs_delta_over_budget": (
                            abs(delta) / budget if budget > 0.0 else np.inf
                        ),
                        "within_budget": bool(abs(delta) <= budget),
                    }
                )
    return pd.DataFrame(rows)


def make_plot(
    path: Path,
    primary: dict[str, dict[str, object]],
    max_z: float,
) -> None:
    z = np.linspace(0.01, max_z, 500)
    lcdm = lcdm_mu(z, float(primary["FLAT_LCDM"]["parameter_value"]))
    anchor_index = 0
    colours = {
        "PETER_NOLOSS_LCONST": "#0072B2",
        "PETER_NOLOSS_CHANDRA": "#D55E00",
        "STATIC_DUALITY_CHANDRA": "#009E73",
        "COMPOSITE_S_FREE": "#CC79A7",
    }
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8))
    for model_key in colours:
        parameter = float(primary[model_key]["parameter_value"])
        delta = model_mu(model_key, z, parameter) - lcdm
        delta -= delta[anchor_index]
        axes[0].plot(z, delta, label=model_key, color=colours[model_key], lw=2)
    axes[0].axhline(0.0, color="black", lw=1)
    axes[0].set_xlabel("Redshift z")
    axes[0].set_ylabel("Anchored Delta mu versus flat Lambda-CDM [mag]")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.25)

    one_plus = 1.0 + z
    axes[1].plot(z, 1.0 / one_plus, label="Q asserted", color="#0072B2", lw=2)
    axes[1].plot(z, np.ones_like(z), label="Q required", color="black", lw=2)
    axes[1].plot(
        z,
        one_plus,
        label="Compton frequency, fixed h",
        color="#D55E00",
        lw=2,
    )
    axes[1].plot(
        z,
        one_plus**3.5,
        label="M_Ch ratio, fixed hbar and G",
        color="#009E73",
        lw=2,
    )
    axes[1].set_yscale("log")
    axes[1].set_xlabel("Redshift z")
    axes[1].set_ylabel("Dimensionless ratio")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.25, which="both")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def markdown_table(frame: pd.DataFrame, columns: list[str]) -> str:
    if frame.empty:
        return "_No rows._"
    return frame[columns].to_markdown(index=False, floatfmt=".6g")


def write_report(
    path: Path,
    primary_rows: pd.DataFrame,
    covariance_rows: pd.DataFrame,
    release_rows: pd.DataFrame,
    split_rows: pd.DataFrame,
    branches: pd.DataFrame,
    theory: pd.DataFrame,
    budgets: pd.DataFrame,
    primary: dict[str, dict[str, object]],
    gates: dict[str, bool],
    protocol_hash: str,
) -> None:
    exact = primary["PETER_NOLOSS_LCONST"]
    composite = primary["COMPOSITE_S_FREE"]
    chandra = primary["PETER_NOLOSS_CHANDRA"]
    s_hat = float(composite["parameter_value"])
    s_sigma = float(composite["parameter_sigma"])
    ell_noloss = -2.0 * s_hat
    ell_duality = 4.0 - 2.0 * s_hat
    sensitivity_exact = covariance_rows[
        covariance_rows["model"] == "PETER_NOLOSS_LCONST"
    ]
    worst_budget = budgets[budgets["candidate"] == "PETER_NOLOSS_LCONST"].nlargest(
        8, "abs_delta_over_budget"
    )
    lines = [
        "# Peter clock-law v1: physical and supernova consistency audit",
        "",
        f"Generated: {datetime.now():%Y-%m-%d %H:%M:%S}",
        "",
        "## Decision",
        "",
        "**THE EXACT `1+z/2` PATH IS IMPLEMENTED CORRECTLY BUT FAILS THE "
        "SUPERNNOVA AND PHYSICAL-CLOSURE GATES. DO NOT ESCALATE TO A "
        "LIGHT-CURVE-LEVEL FIT.**".replace("SUPERNNOVA", "SUPERNOVA"),
        "",
        f"The exact no-loss, constant-luminosity branch gives `Delta BIC="
        f"{float(exact['delta_BIC_vs_LCDM']):+.6f}` relative to flat Lambda-CDM "
        f"for `N={int(primary_rows.iloc[0]['N'])}`. Positive values favour "
        "Lambda-CDM. This reproduces the already-known fixed-law result and is "
        "not new evidence by itself.",
        "",
        "More importantly, the proposed atomic-clock statement is inconsistent "
        "with the registered spectroscopic-redshift relation. With `a=1` and "
        "`C=1+z`, the required atomic-frequency ratio is `Q=1`, not "
        "`Q=(1+z)^-1`. With fixed `h`, the stated mass and light-speed laws also "
        "make the Compton frequency grow as `1+z`, rather than fall as "
        "`(1+z)^-1`.",
        "",
        "## Frozen status",
        "",
        f"Protocol SHA-256: `{protocol_hash}`.",
        "",
        "The earlier exact-path Delta BIC was disclosed in the protocol. The "
        "newly frozen content was the physical closure, deterministic source and "
        "photon branches, covariance sensitivities and escalation rule.",
        "",
        "## Registered equations",
        "",
        "$$",
        r"\begin{aligned}",
        r"a(z)                          &= 1, \\",
        r"C(z)                          &= 1+z, \\",
        r"R(z)                          &= H_0, \\",
        r"\chi(z)                       &= \frac{c_0}{H_0}z\left(1+\frac{z}{2}\right), \\",
        r"1+z                           &= \frac{C(z)}{a(z)Q(z)}, \\",
        r"Q_{\rm required}(z)           &= 1, \\",
        r"D_{\rm eff}(z)                &= \chi(z)(1+z)^s, \\",
        r"s                             &= \tau-\frac{\ell}{2}, \\",
        r"L_{\rm Ia}(z)/L_{\rm Ia,0}    &= (1+z)^\ell.",
        "\\end{aligned}",
        "$$",
        "",
        "The supernova likelihood identifies only the composite `s`. It cannot "
        "separate photon-distance transfer `tau` from intrinsic luminosity "
        "evolution `ell`.",
        "",
        "## Primary likelihood",
        "",
        markdown_table(
            primary_rows,
            [
                "model",
                "parameter_name",
                "parameter_value",
                "parameter_sigma",
                "chi2",
                "BIC",
                "delta_BIC_vs_LCDM",
            ],
        ),
        "",
        f"The diagnostic composite fit gives `s={s_hat:.6f} +/- {s_sigma:.6f}`. "
        f"It is equivalent to `ell={ell_noloss:.6f} +/- {2*s_sigma:.6f}` under "
        f"the no-loss interpretation, or `ell={ell_duality:.6f} +/- "
        f"{2*s_sigma:.6f}` under the static distance-duality interpretation. "
        "These are exactly degenerate descriptions of the same magnitude curve.",
        "",
        "## Deterministic branches",
        "",
        markdown_table(
            branches,
            [
                "model",
                "tau_photon_distance_power",
                "ell_luminosity_power",
                "s_composite_power",
                "description",
            ],
        ),
        "",
        f"With fixed `hbar`, `G` and particle masses proportional to `(1+z)^-1`, "
        "the Chandrasekhar proxy predicts `M_Ch/M_Ch,0=(1+z)^(7/2)`. The "
        f"corresponding no-loss branch gives `Delta BIC="
        f"{float(chandra['delta_BIC_vs_LCDM']):+.6f}`.",
        "",
        "## Physical consistency readout",
        "",
        markdown_table(
            theory,
            [
                "z_observed",
                "C_c_over_c0",
                "Q_required_by_spectroscopic_closure",
                "Q_asserted_clock_control",
                "z_predicted_if_asserted_Q_used",
                "Compton_frequency_ratio_fixed_h",
                "M_Ch_ratio_fixed_hbar_G",
                "N_Ch_ratio_fixed_hbar_G",
            ],
        ),
        "",
        "The identity `integral(1+z) dz = z(1+z/2)` passes. The failure is not "
        "that integral; it is the additional physical identification of redshift, "
        "atomic frequency, source luminosity and observed flux.",
        "",
        "## Covariance hierarchy sensitivities",
        "",
        markdown_table(
            sensitivity_exact,
            [
                "covariance_variant",
                "N",
                "chi2",
                "delta_BIC_vs_LCDM",
            ],
        ),
        "",
        f"The exact-branch sensitivity range is "
        f"[{sensitivity_exact['delta_BIC_vs_LCDM'].min():+.6f}, "
        f"{sensitivity_exact['delta_BIC_vs_LCDM'].max():+.6f}]. The released "
        "total covariance remains primary. Grouped reconstructions are "
        "diagnostics because they omit substantial covariance structure.",
        "",
        "## Release and redshift diagnostics",
        "",
        "### Release splits",
        "",
        markdown_table(
            release_rows[
                release_rows["model"].isin(
                    ["PETER_NOLOSS_LCONST", "COMPOSITE_S_FREE", "FLAT_LCDM"]
                )
            ],
            [
                "scope",
                "model",
                "N",
                "parameter_value",
                "chi2",
                "delta_BIC_vs_LCDM",
            ],
        ),
        "",
        "### Descriptive redshift splits",
        "",
        markdown_table(
            split_rows[
                split_rows["model"].isin(
                    ["PETER_NOLOSS_LCONST", "COMPOSITE_S_FREE", "FLAT_LCDM"]
                )
            ],
            [
                "scope",
                "model",
                "N",
                "parameter_value",
                "chi2",
                "delta_BIC_vs_LCDM",
            ],
        ),
        "",
        "These subset cells profile their own release intercepts and are "
        "descriptive, not independent predictive tests.",
        "",
        "## Calibration-budget modes",
        "",
        markdown_table(
            worst_budget,
            [
                "dataset",
                "group",
                "N",
                "budget_sigma_mag",
                "delta_mode_mag",
                "abs_delta_over_budget",
                "within_budget",
            ],
        ),
        "",
        "Release-contained calibration modes are not added again to `STAT+SYS`; "
        "the table compares model-difference residual modes with the projected "
        "release budgets.",
        "",
        "## Gates",
        "",
    ]
    for name, passed in gates.items():
        lines.append(f"- `{name}`: **{'PASS' if passed else 'FAIL'}**")
    lines.extend(
        [
            "",
            "## Claim boundary",
            "",
            "This result rejects promotion of the explicit v1 closure. It does "
            "not reject every FR, varying-constant or static-universe theory. A "
            "revised theory must derive, before another cosmological fit, the "
            "redshift identity, dimensionless atomic-frequency law, photon "
            "transfer, Type Ia luminosity and the dynamics of the clock field.",
            "",
            "The proposed next stage was a Pantheon+ light-curve-level refit. "
            "Because the frozen escalation rule fails, that stage is not run.",
            "",
            "## Reproduction",
            "",
            "```powershell",
            "python github_export/code/rawmu/audit_peter_clock_law_v1_2026_07_19.py",
            "```",
            "",
            "## Primary data references",
            "",
            "- [Riess et al. SH0ES distance-ladder equations](https://arxiv.org/abs/2112.04510)",
            "- [Pantheon+ full data and light-curve release](https://arxiv.org/abs/2112.03863)",
            "- [Official Pantheon+ data release](https://github.com/PantheonPlusSH0ES/DataRelease)",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def write_manifest(outdir: Path, paths: list[Path]) -> Path:
    rows = []
    for path in sorted(set(paths), key=lambda item: item.as_posix()):
        rows.append(
            {
                "path": (
                    path.relative_to(REPO_ROOT).as_posix()
                    if path.is_relative_to(REPO_ROOT)
                    else str(path)
                ),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    manifest = outdir / f"peter_clock_law_v1_manifest_{DATE_TAG}.csv"
    pd.DataFrame(rows).to_csv(manifest, index=False)
    return manifest


def self_test() -> None:
    z = np.array([0.01, 0.2, 1.0, 2.0])
    expected = (C_KMS / H0_FIXED) * z * (1.0 + 0.5 * z)
    np.testing.assert_allclose(peter_path_distance(z), expected, rtol=2e-14, atol=2e-12)
    covariance = np.diag(np.array([0.01, 0.012, 0.009, 0.011]) ** 2)
    truth = 0.37
    block = grounded.make_block(
        "synthetic",
        z,
        peter_effective_mu(z, truth) + 0.12,
        covariance,
        np.array([1, 1, 1, 1]),
        np.arange(len(z)),
        "synthetic",
    )
    measured, sigma, chi2, offsets = analytic_composite_fit([block])
    if abs(measured - truth) > 1e-10 or sigma <= 0.0 or chi2 > 1e-16:
        raise AssertionError("Analytic composite fit failed")
    if abs(offsets["synthetic"] - 0.12) > 1e-10:
        raise AssertionError("Profiled intercept recovery failed")
    if BRANCH_BY_KEY["PETER_NOLOSS_CHANDRA"].s != -1.75:
        raise AssertionError("Chandrasekhar branch power changed")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--protocol", type=Path, default=DEFAULT_PROTOCOL)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument("--min-survey-n", type=int, default=20)
    args = parser.parse_args()

    self_test()
    protocol = verify_protocol(args.protocol)
    if not args.registry.exists():
        raise FileNotFoundError(args.registry)
    registry = json.loads(args.registry.read_text(encoding="utf-8"))
    args.outdir.mkdir(parents=True, exist_ok=True)

    variants, load_metadata = grounded.load_blocks(args.min_survey_n)
    primary_blocks = variants["released_total_primary"]
    all_rows: list[dict[str, object]] = []
    primary_internal: dict[str, dict[str, object]] | None = None
    for variant, blocks in variants.items():
        rows, internal = fit_models(blocks, "all_releases", variant)
        all_rows.extend(rows)
        if variant == "released_total_primary":
            primary_internal = internal
        exact = internal["PETER_NOLOSS_LCONST"]
        print(
            f"{variant}: exact Delta BIC="
            f"{float(exact['delta_BIC_vs_LCDM']):+.6f}",
            flush=True,
        )
    assert primary_internal is not None
    covariance_rows = pd.DataFrame(all_rows)
    primary_rows = covariance_rows[
        covariance_rows["covariance_variant"] == "released_total_primary"
    ].copy()

    release_rows_list: list[dict[str, object]] = []
    for block in primary_blocks:
        rows, _ = fit_models(
            [block], f"release:{block.label}", "released_total_primary"
        )
        release_rows_list.extend(rows)
    release_rows = pd.DataFrame(release_rows_list)

    low_blocks = subset_blocks(
        primary_blocks, lambda values: values < 0.5, "z_lt_0p5"
    )
    high_blocks = subset_blocks(
        primary_blocks, lambda values: values >= 0.5, "z_ge_0p5"
    )
    split_rows_list: list[dict[str, object]] = []
    for scope, blocks in (("z<0.5", low_blocks), ("z>=0.5", high_blocks)):
        rows, _ = fit_models(blocks, scope, "released_total_primary")
        split_rows_list.extend(rows)
    split_rows = pd.DataFrame(split_rows_list)

    max_z = max(float(np.max(block.z)) for block in primary_blocks)
    theory = theoretical_readout(max_z)
    branches = fixed_branch_table()
    budgets = calibration_budget_table(
        primary_blocks, primary_internal, registry, args.min_survey_n
    )

    exact = primary_internal["PETER_NOLOSS_LCONST"]
    chandra = primary_internal["PETER_NOLOSS_CHANDRA"]
    composite = primary_internal["COMPOSITE_S_FREE"]
    exact_sensitivities = covariance_rows[
        covariance_rows["model"] == "PETER_NOLOSS_LCONST"
    ]
    exact_release = release_rows[
        release_rows["model"] == "PETER_NOLOSS_LCONST"
    ]
    q_difference = np.max(
        np.abs(
            theory["Q_required_by_spectroscopic_closure"].to_numpy()
            - theory["Q_asserted_clock_control"].to_numpy()
        )
    )
    compton_difference = np.max(
        np.abs(
            theory["Compton_frequency_ratio_fixed_h"].to_numpy()
            - theory["Q_asserted_clock_control"].to_numpy()
        )
    )
    regression_difference = abs(
        float(exact["delta_BIC_vs_LCDM"]) - KNOWN_EXACT_DELTA_BIC
    )
    s_hat = float(composite["parameter_value"])
    s_sigma = float(composite["parameter_sigma"])
    fixed_s_values = np.array([branch.s for branch in FIXED_BRANCHES])
    minimum_fixed_s_tension = float(np.min(np.abs(fixed_s_values - s_hat) / s_sigma))
    exact_budget = budgets[budgets["candidate"] == "PETER_NOLOSS_LCONST"]
    gates = {
        "frozen_protocol_hash_verified": True,
        "analytic_path_identity": True,
        "known_fixed_likelihood_regression_within_1e-6": bool(
            regression_difference <= 1e-6
        ),
        "exact_path_absolute_delta_BIC_at_most_2": bool(
            abs(float(exact["delta_BIC_vs_LCDM"])) <= 2.0
        ),
        "exact_path_favours_over_LCDM": bool(
            float(exact["delta_BIC_vs_LCDM"]) < 0.0
        ),
        "exact_path_all_covariance_sensitivities_abs_delta_BIC_at_most_2": bool(
            (exact_sensitivities["delta_BIC_vs_LCDM"].abs() <= 2.0).all()
        ),
        "exact_path_every_release_abs_delta_BIC_at_most_2": bool(
            (exact_release["delta_BIC_vs_LCDM"].abs() <= 2.0).all()
        ),
        "asserted_Q_matches_spectroscopic_closure": bool(q_difference <= 1e-10),
        "fixed_h_Compton_frequency_matches_asserted_clock": bool(
            compton_difference <= 1e-10
        ),
        "Chandrasekhar_proxy_absolute_delta_BIC_at_most_2": bool(
            abs(float(chandra["delta_BIC_vs_LCDM"])) <= 2.0
        ),
        "diagnostic_composite_power_is_interior": bool(
            S_BOUNDS[0] < s_hat < S_BOUNDS[1]
        ),
        "diagnostic_composite_within_2sigma_of_registered_fixed_branch": bool(
            minimum_fixed_s_tension <= 2.0
        ),
        "exact_model_difference_modes_within_release_budgets": bool(
            not exact_budget.empty and exact_budget["within_budget"].all()
        ),
        "photon_transfer_and_source_luminosity_separately_predicted": False,
        "action_level_clock_matter_electromagnetic_model_supplied": False,
        "no_release_covariance_prior_double_counting": True,
    }
    theory_gate_names = (
        "asserted_Q_matches_spectroscopic_closure",
        "fixed_h_Compton_frequency_matches_asserted_clock",
        "Chandrasekhar_proxy_absolute_delta_BIC_at_most_2",
        "photon_transfer_and_source_luminosity_separately_predicted",
        "action_level_clock_matter_electromagnetic_model_supplied",
    )
    gates["escalate_to_light_curve_level"] = bool(
        gates["exact_path_absolute_delta_BIC_at_most_2"]
        and all(gates[name] for name in theory_gate_names)
    )

    fit_path = args.outdir / f"peter_clock_law_v1_fits_{DATE_TAG}.csv"
    release_path = args.outdir / f"peter_clock_law_v1_release_splits_{DATE_TAG}.csv"
    split_path = args.outdir / f"peter_clock_law_v1_redshift_splits_{DATE_TAG}.csv"
    theory_path = args.outdir / f"peter_clock_law_v1_physical_readout_{DATE_TAG}.csv"
    branch_path = args.outdir / f"peter_clock_law_v1_fixed_branches_{DATE_TAG}.csv"
    budget_path = args.outdir / f"peter_clock_law_v1_budget_modes_{DATE_TAG}.csv"
    summary_path = args.outdir / f"peter_clock_law_v1_summary_{DATE_TAG}.json"
    report_path = args.outdir / f"peter_clock_law_v1_report_{DATE_TAG}.md"
    plot_path = args.outdir / f"peter_clock_law_v1_readout_{DATE_TAG}.png"

    covariance_rows.to_csv(fit_path, index=False)
    release_rows.to_csv(release_path, index=False)
    split_rows.to_csv(split_path, index=False)
    theory.to_csv(theory_path, index=False)
    branches.to_csv(branch_path, index=False)
    budgets.to_csv(budget_path, index=False)
    summary = {
        "analysis_id": ANALYSIS_ID,
        "generated": datetime.now().isoformat(timespec="seconds"),
        "protocol_sha256": sha256_file(args.protocol),
        "protocol": protocol,
        "sample_metadata": load_metadata,
        "primary": {
            model: {
                key: value
                for key, value in fit.items()
                if key != "offsets"
            }
            | {"offsets": fit["offsets"]}
            for model, fit in primary_internal.items()
        },
        "derived": {
            "maximum_Q_absolute_mismatch": float(q_difference),
            "maximum_Compton_clock_absolute_mismatch": float(compton_difference),
            "known_regression_difference": float(regression_difference),
            "minimum_registered_fixed_s_tension_sigma": minimum_fixed_s_tension,
            "equivalent_ell_no_loss": -2.0 * s_hat,
            "equivalent_ell_static_duality": 4.0 - 2.0 * s_hat,
            "exact_covariance_sensitivity_delta_BIC_range": [
                float(exact_sensitivities["delta_BIC_vs_LCDM"].min()),
                float(exact_sensitivities["delta_BIC_vs_LCDM"].max()),
            ],
        },
        "gates": gates,
        "decision": (
            "ESCALATE_TO_LIGHT_CURVE_LEVEL"
            if gates["escalate_to_light_curve_level"]
            else "DO_NOT_ESCALATE_EXPLICIT_V1_CLOSURE_FAILS"
        ),
        "claim_boundary": protocol["claim_boundary"],
    }
    summary_path.write_text(
        json.dumps(safe_json(summary), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    make_plot(plot_path, primary_internal, max_z)
    write_report(
        report_path,
        primary_rows,
        covariance_rows,
        release_rows,
        split_rows,
        branches,
        theory,
        budgets,
        primary_internal,
        gates,
        sha256_file(args.protocol),
    )
    manifest = write_manifest(
        args.outdir,
        [
            SCRIPT_PATH,
            args.protocol,
            args.outdir / f"peter_clock_law_v1_protocol_{DATE_TAG}.md",
            args.registry,
            SHARED_CODE / "plamb_clock_distance.py",
            RAWMU_CODE / "run_rawmu_release_grounded_holdouts_2026_07_18.py",
            fit_path,
            release_path,
            split_path,
            theory_path,
            branch_path,
            budget_path,
            summary_path,
            report_path,
            plot_path,
        ],
    )
    print(f"Decision: {summary['decision']}", flush=True)
    print(f"Saved report: {report_path}", flush=True)
    print(f"Saved manifest: {manifest}", flush=True)


if __name__ == "__main__":
    main()
