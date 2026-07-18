#!/usr/bin/env python3
r"""Injection-and-recovery calibration for the locked SPARC group test.

This is a post-result power analysis. It does not test new physical data and it
does not access the five reserved replication galaxies.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[3]
RUN_DATE = "2026-07-18"
SEED = 2026071803
FOLD_SEED = 2026071802
RIDGE_ALPHA = 1.0
OUTER_FOLDS = 10
REPLICATES_PER_EFFECT = 5000
NULL_PERMUTATIONS = 100000
TARGET_ABS_RESIDUAL_CORRELATIONS = (0.0, 0.1, 0.2, 0.3, 0.307254, 0.4, 0.5, 0.6)
PRIOR_EFFECT_ANCHOR = 0.307254

GROUP_RESULT_DIR = (
    ROOT
    / "github_export"
    / "results"
    / "2026-07-18"
    / "fr_group_catalogue_replication"
)
PRIMARY_FRAME = GROUP_RESULT_DIR / "sparc_fr_group_catalogue_primary_frame.csv"
GROUP_PREREGISTRATION = GROUP_RESULT_DIR / "sparc_fr_group_catalogue_preregistration.json"
GROUP_GATE = GROUP_RESULT_DIR / "sparc_fr_group_catalogue_gate.json"
GROUP_SUMMARY = GROUP_RESULT_DIR / "sparc_fr_group_catalogue_test_summary.csv"
GROUP_PROGRAMME = (
    ROOT
    / "github_export"
    / "code"
    / "sparc"
    / "run_sparc_fr_group_catalogue_replication_2026-07-18.py"
)
PRIOR_ENVIRONMENT_SUMMARY = (
    ROOT
    / "github_export"
    / "results"
    / "2026-07-18"
    / "fr_environment_asymmetry"
    / "sparc_fr_environment_test_summary.csv"
)
DEFAULT_OUT = ROOT / "plamb_runs" / "diagnostics" / "sparc_fr_group_power_calibration_20260718"
DEFAULT_EXPORT = (
    ROOT
    / "github_export"
    / "results"
    / "2026-07-18"
    / "fr_group_power_calibration"
)

RESERVED_REPLICATION = {
    "NGC2841",
    "UGC02953",
    "UGC07399",
    "UGC00128",
    "NGC2403",
}
CANDIDATES = {
    "UGC06787",
    "UGC03580",
    "UGC02487",
    "NGC5985",
    "UGC09133",
    "NGC2903",
}
HOST_FEATURES = (
    "T",
    "log_L36",
    "log_SBeff",
    "log_MHI",
    "log_Vflat",
    "Inc_deg",
    "frac_distance_error",
    "bulge_proxy",
    "log_n_points",
    "D_Mpc",
    "abs_gal_b_deg",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def host_matrix(frame: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    matrix = pd.DataFrame(index=frame.index)
    for column in HOST_FEATURES:
        values = pd.to_numeric(frame[column], errors="coerce")
        matrix[column] = values.fillna(float(values.median()))
    distance_dummies = pd.get_dummies(
        pd.to_numeric(frame["f_D"], errors="coerce").fillna(-1).astype(int),
        prefix="f_D",
        dtype=float,
    )
    matrix = pd.concat([matrix, distance_dummies], axis=1)
    return matrix.to_numpy(dtype=float), list(matrix.columns)


def fixed_folds(n_rows: int) -> list[tuple[np.ndarray, np.ndarray]]:
    splitter = KFold(n_splits=OUTER_FOLDS, shuffle=True, random_state=FOLD_SEED)
    return [(train, test) for train, test in splitter.split(np.arange(n_rows))]


def ridge_operators(
    train_design: np.ndarray,
    test_design: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    mean = np.mean(train_design, axis=0)
    scale = np.std(train_design, axis=0, ddof=0)
    scale[scale == 0.0] = 1.0
    train_scaled = (train_design - mean) / scale
    test_scaled = (test_design - mean) / scale
    n_train = len(train_design)
    centring = np.eye(n_train) - np.ones((n_train, n_train)) / n_train
    inverse = np.linalg.pinv(
        train_scaled.T @ train_scaled + RIDGE_ALPHA * np.eye(train_scaled.shape[1])
    )
    coefficient_operator = inverse @ train_scaled.T @ centring
    prediction_operator = (
        test_scaled @ coefficient_operator
        + np.ones((len(test_design), n_train)) / n_train
    )
    return prediction_operator, coefficient_operator


def prepare_operators(
    x: np.ndarray,
    environment: np.ndarray,
) -> tuple[list[dict[str, Any]], np.ndarray]:
    operators: list[dict[str, Any]] = []
    environment_prediction = np.full(len(environment), np.nan)
    for train, test in fixed_folds(len(x)):
        baseline_prediction, _ = ridge_operators(x[train], x[test])
        augmented_design = np.column_stack([x, environment])
        augmented_prediction, augmented_coefficients = ridge_operators(
            augmented_design[train], augmented_design[test]
        )
        environment_prediction[test] = baseline_prediction @ environment[train]
        operators.append(
            {
                "train": train,
                "test": test,
                "baseline_prediction": baseline_prediction,
                "augmented_prediction": augmented_prediction,
                "environment_coefficient": augmented_coefficients[-1],
            }
        )
    return operators, environment - environment_prediction


def full_ridge_signal(design: np.ndarray, response: np.ndarray) -> np.ndarray:
    scaler = StandardScaler().fit(design)
    model = Ridge(alpha=RIDGE_ALPHA).fit(scaler.transform(design), response)
    return model.predict(scaler.transform(design))


def standardise(values: np.ndarray) -> np.ndarray:
    centred = values - np.mean(values)
    scale = np.std(centred, ddof=1)
    if scale <= 0.0:
        raise ValueError("Cannot standardise a constant vector")
    return centred / scale


def spearman_rho(x: np.ndarray, y: np.ndarray) -> float:
    xr = rankdata(x, method="average").astype(float)
    yr = rankdata(y, method="average").astype(float)
    xr -= np.mean(xr)
    yr -= np.mean(yr)
    denominator = math.sqrt(float(np.sum(xr**2) * np.sum(yr**2)))
    return float(np.dot(xr, yr) / denominator)


def null_critical_rho(
    environment_residual: np.ndarray,
    rng: np.random.Generator,
) -> tuple[float, np.ndarray]:
    xr = rankdata(environment_residual, method="average").astype(float)
    xr -= np.mean(xr)
    base = np.arange(len(xr), dtype=float)
    base -= np.mean(base)
    denominator = math.sqrt(float(np.sum(xr**2) * np.sum(base**2)))
    values = np.empty(NULL_PERMUTATIONS, dtype=float)
    for index in range(NULL_PERMUTATIONS):
        values[index] = float(np.dot(xr, rng.permutation(base)) / denominator)
    critical = float(np.quantile(np.abs(values), 0.99, method="higher"))
    return critical, values


def evaluate_response(
    response: np.ndarray,
    environment_residual: np.ndarray,
    operators: list[dict[str, Any]],
    critical_rho: float,
) -> dict[str, Any]:
    baseline_prediction = np.full(len(response), np.nan)
    augmented_prediction = np.full(len(response), np.nan)
    coefficients: list[float] = []
    for item in operators:
        train = item["train"]
        test = item["test"]
        baseline_prediction[test] = item["baseline_prediction"] @ response[train]
        augmented_prediction[test] = item["augmented_prediction"] @ response[train]
        coefficients.append(float(item["environment_coefficient"] @ response[train]))
    outcome_residual = response - baseline_prediction
    rho = spearman_rho(environment_residual, outcome_residual)
    baseline_rmse = math.sqrt(float(np.mean((response - baseline_prediction) ** 2)))
    augmented_rmse = math.sqrt(float(np.mean((response - augmented_prediction) ** 2)))
    improvement = (baseline_rmse - augmented_rmse) / baseline_rmse
    coefficients_array = np.asarray(coefficients, dtype=float)
    sign_fraction = max(
        float(np.mean(coefficients_array > 0.0)),
        float(np.mean(coefficients_array < 0.0)),
    )
    p_gate = abs(rho) >= critical_rho
    cv_gate = improvement >= 0.05
    sign_gate = sign_fraction >= 0.80
    return {
        "partial_spearman_rho": rho,
        "cv_rmse_fractional_improvement": improvement,
        "coefficient_sign_fraction": sign_fraction,
        "p_gate_passed": p_gate,
        "cv_gate_passed": cv_gate,
        "sign_gate_passed": sign_gate,
        "joint_primary_gate_passed": p_gate and cv_gate and sign_gate,
    }


def verify_observed_metric_parity(frame: pd.DataFrame) -> dict[str, float]:
    x, _ = host_matrix(frame)
    environment = frame["group_richness_score"].to_numpy(dtype=float)
    response = frame["signed_log1p_delta_per_point"].to_numpy(dtype=float)
    operators, environment_residual = prepare_operators(x, environment)
    reproduced = evaluate_response(
        response,
        environment_residual,
        operators,
        critical_rho=float("inf"),
    )
    summary = pd.read_csv(GROUP_SUMMARY)
    locked = summary[
        (summary["sample"] == "primary_10_40_mpc")
        & (summary["scenario"] == "combined_conventional")
        & (summary["control"] == "actual")
    ].iloc[0]
    comparisons = {
        "partial_spearman_rho": float(locked["partial_spearman_rho"]),
        "cv_rmse_fractional_improvement": float(
            locked["cv_rmse_fractional_improvement"]
        ),
        "coefficient_sign_fraction": float(
            locked["environment_coefficient_sign_fraction"]
        ),
    }
    for key, expected in comparisons.items():
        if not math.isclose(float(reproduced[key]), expected, rel_tol=0.0, abs_tol=1.0e-12):
            raise RuntimeError(
                f"Calibration does not reproduce locked {key}: "
                f"{reproduced[key]} versus {expected}"
            )
    return {key: float(reproduced[key]) for key in comparisons}


def preregistration_config(programme: Path) -> dict[str, Any]:
    frame = pd.read_csv(PRIMARY_FRAME)
    excluded = set(frame["Galaxy"].astype(str)) & (RESERVED_REPLICATION | CANDIDATES)
    if excluded:
        raise AssertionError(f"Excluded galaxies in calibration frame: {sorted(excluded)}")
    if len(frame) != 65:
        raise ValueError(f"Expected locked N=65 primary frame, found {len(frame)}")
    parity = verify_observed_metric_parity(frame)
    return {
        "date": RUN_DATE,
        "locked_utc": utc_now(),
        "programme": str(programme),
        "programme_sha256": sha256(programme),
        "status": "Post-result injection-and-recovery power calibration; not a physical-data test",
        "input_hashes": {
            "primary_frame": sha256(PRIMARY_FRAME),
            "group_preregistration": sha256(GROUP_PREREGISTRATION),
            "group_gate": sha256(GROUP_GATE),
            "group_summary": sha256(GROUP_SUMMARY),
            "group_programme": sha256(GROUP_PROGRAMME),
            "prior_environment_summary": sha256(PRIOR_ENVIRONMENT_SUMMARY),
        },
        "sample": {
            "N": len(frame),
            "excluded_previous_candidates": sorted(CANDIDATES),
            "excluded_reserved_replication": sorted(RESERVED_REPLICATION),
        },
        "observed_metric_parity": parity,
        "design": {
            "predictor": "group_richness_score = ln(N_group)",
            "host_features": list(HOST_FEATURES) + ["f_D one-hot indicators"],
            "host_signal": "full-sample ridge prediction with alpha=1",
            "empirical_noise": (
                "centred and standardised residuals from the full-sample host-only ridge model, "
                "randomly permuted without replacement in every replicate"
            ),
            "environment_injection": (
                "host-residualised standardised group richness multiplied by "
                "beta=r/sqrt(1-r^2), with the negative sign matching the earlier 2MRS direction"
            ),
            "target_absolute_residual_correlations": list(TARGET_ABS_RESIDUAL_CORRELATIONS),
            "prior_effect_anchor": PRIOR_EFFECT_ANCHOR,
        },
        "monte_carlo": {
            "seed": SEED,
            "locked_fold_seed": FOLD_SEED,
            "replicates_per_effect": REPLICATES_PER_EFFECT,
            "null_permutations_for_rho_critical_value": NULL_PERMUTATIONS,
            "folds": OUTER_FOLDS,
            "ridge_alpha": RIDGE_ALPHA,
        },
        "recovery": {
            "p_component": "absolute partial Spearman rho exceeds the calibrated two-sided 1 per cent null threshold",
            "cv_component": "cross-fitted RMSE fractional improvement >=0.05",
            "sign_component": "one environment-coefficient sign in at least 8/10 folds",
            "joint_primary_gate": "all three primary components pass in the same replicate",
        },
        "locked_interpretation": {
            "adequate_power_at_prior_anchor": "joint primary recovery probability >=0.80 at r=0.307254",
            "underpowered_at_prior_anchor": "joint primary recovery probability <0.80 at r=0.307254",
            "detection_threshold": "smallest tested r with joint primary recovery probability >=0.80",
            "claim_boundary": (
                "Power calibrates the failed group-richness test only. It cannot create evidence for FR, "
                "antimatter or an environmental effect."
            ),
        },
    }


def write_preregistration(path: Path, config: dict[str, Any]) -> None:
    grid = ", ".join(f"{value:.6g}" for value in TARGET_ABS_RESIDUAL_CORRELATIONS)
    lines = [
        "# SPARC FR Group-Richness Injection-and-Recovery Preregistration",
        "",
        f"Date: {RUN_DATE}",
        f"Locked: `{config['locked_utc']}`",
        "",
        "## Status",
        "",
        "This is a post-result power calibration of the failed 65-galaxy published-group test. It is not a new physical-data test and cannot rescue or alter the failed gate.",
        "",
        "## Fixed Design",
        "",
        "The exact locked primary frame, group-richness predictor, eleven intrinsic host controls, distance-method indicators, ten folds and ridge penalty are reused. The six previous candidates and five reserved galaxies are absent.",
        "",
        "A full-sample host-only ridge fit defines the deterministic host signal. Its centred, standardised empirical residuals are permuted without replacement to generate noise. The host-residualised group-richness vector receives a negative injected effect, matching the direction of the earlier 2MRS composite.",
        "",
        f"Target absolute residual correlations are: `{grid}`. The exact earlier-effect anchor is `r={PRIOR_EFFECT_ANCHOR:.6f}`.",
        "",
        "For target r, the injected standardised slope is beta=r/sqrt(1-r^2). There are 5,000 replicates at every grid value.",
        "",
        "## Recovery Gates",
        "",
        "The two-sided p<=0.01 component is represented by an absolute Spearman-rho threshold calibrated from 100,000 null rank permutations of the fixed environment residual. Each simulation also reruns the locked cross-fitted RMSE and coefficient-sign calculations.",
        "",
        "A joint primary recovery requires the p component, RMSE improvement>=5 per cent and one coefficient sign in at least 8/10 folds in the same replicate.",
        "",
        "## Locked Interpretation",
        "",
        "Power is adequate for the earlier-effect anchor only if at least 80 per cent of anchor injections recover the joint primary gate. The detection threshold is the smallest tested absolute residual correlation with at least 80 per cent joint recovery.",
        "",
        "A low recovery rate means the group-richness design is too weak or too coarse to reject an effect of that standardised magnitude. A high recovery rate would make the observed null informative against an effect carried by group richness. Neither outcome is evidence for antimatter.",
        "",
        f"Programme SHA-256: `{config['programme_sha256']}`",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def run_simulations(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    x, feature_names = host_matrix(frame)
    environment = frame["group_richness_score"].to_numpy(dtype=float)
    observed_response = frame["signed_log1p_delta_per_point"].to_numpy(dtype=float)
    operators, crossfit_environment_residual = prepare_operators(x, environment)
    host_signal = full_ridge_signal(x, observed_response)
    empirical_noise = standardise(observed_response - host_signal)
    environment_signal = standardise(environment - full_ridge_signal(x, environment))
    noise_scale = float(np.std(observed_response - host_signal, ddof=1))
    rng = np.random.default_rng(SEED)
    critical_rho, null_rhos = null_critical_rho(crossfit_environment_residual, rng)
    rows: list[dict[str, Any]] = []
    for target in TARGET_ABS_RESIDUAL_CORRELATIONS:
        beta = 0.0 if target == 0.0 else -target / math.sqrt(1.0 - target**2)
        for replicate in range(1, REPLICATES_PER_EFFECT + 1):
            response = host_signal + noise_scale * (
                rng.permutation(empirical_noise) + beta * environment_signal
            )
            result = evaluate_response(
                response,
                crossfit_environment_residual,
                operators,
                critical_rho,
            )
            result.update(
                {
                    "target_abs_residual_correlation": target,
                    "injected_standardised_beta": beta,
                    "replicate": replicate,
                }
            )
            rows.append(result)
    draws = pd.DataFrame(rows)
    summary_rows: list[dict[str, Any]] = []
    for target, selected in draws.groupby("target_abs_residual_correlation", sort=True):
        summary_rows.append(
            {
                "target_abs_residual_correlation": target,
                "injected_standardised_beta": float(selected["injected_standardised_beta"].iloc[0]),
                "replicates": len(selected),
                "median_recovered_rho": float(selected["partial_spearman_rho"].median()),
                "rho_q025": float(selected["partial_spearman_rho"].quantile(0.025)),
                "rho_q975": float(selected["partial_spearman_rho"].quantile(0.975)),
                "p_component_recovery": float(selected["p_gate_passed"].mean()),
                "median_cv_rmse_fractional_improvement": float(
                    selected["cv_rmse_fractional_improvement"].median()
                ),
                "cv_component_recovery": float(selected["cv_gate_passed"].mean()),
                "sign_component_recovery": float(selected["sign_gate_passed"].mean()),
                "joint_primary_recovery": float(selected["joint_primary_gate_passed"].mean()),
            }
        )
    summary = pd.DataFrame(summary_rows)
    observed = evaluate_response(
        observed_response,
        crossfit_environment_residual,
        operators,
        critical_rho,
    )
    diagnostics = {
        "date": RUN_DATE,
        "N": len(frame),
        "host_features": feature_names,
        "critical_absolute_spearman_rho_for_p_le_0p01": critical_rho,
        "null_permutations": NULL_PERMUTATIONS,
        "null_rho_abs_q99": float(np.quantile(np.abs(null_rhos), 0.99)),
        "empirical_noise_sd_outcome_units": noise_scale,
        "observed_locked_primary": observed,
    }
    return draws, summary, diagnostics


def make_plot(summary: pd.DataFrame, path: Path) -> None:
    x = summary["target_abs_residual_correlation"].to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(7.4, 5.2))
    for column, label, colour, marker in [
        ("p_component_recovery", "p<=0.01", "#2f6f6d", "o"),
        ("cv_component_recovery", "CV improvement>=5%", "#b54435", "s"),
        ("sign_component_recovery", "sign>=8/10", "#586994", "^"),
        ("joint_primary_recovery", "joint primary gate", "#222222", "D"),
    ]:
        ax.plot(
            x,
            summary[column].to_numpy(dtype=float),
            marker=marker,
            linewidth=2.0,
            markersize=5,
            label=label,
            color=colour,
        )
    ax.axhline(0.80, color="#777777", linestyle="--", linewidth=1.0)
    ax.axvline(PRIOR_EFFECT_ANCHOR, color="#777777", linestyle=":", linewidth=1.0)
    ax.set_xlim(-0.01, max(x) + 0.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("Injected absolute host-residual correlation")
    ax.set_ylabel("Recovery probability")
    ax.set_title("SPARC Group-Richness Injection-and-Recovery Power")
    ax.grid(alpha=0.18)
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def markdown_table(frame: pd.DataFrame, columns: list[str]) -> str:
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for _, row in frame.iterrows():
        values: list[str] = []
        for column in columns:
            value = row[column]
            if isinstance(value, (float, np.floating)):
                values.append(f"{float(value):.6g}" if math.isfinite(float(value)) else "")
            else:
                values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def write_report(
    path: Path,
    summary: pd.DataFrame,
    diagnostics: dict[str, Any],
) -> None:
    anchor = summary[
        np.isclose(summary["target_abs_residual_correlation"], PRIOR_EFFECT_ANCHOR)
    ].iloc[0]
    adequate = float(anchor["joint_primary_recovery"]) >= 0.80
    qualifying = summary[summary["joint_primary_recovery"] >= 0.80]
    threshold = (
        float(qualifying["target_abs_residual_correlation"].min())
        if len(qualifying)
        else float("nan")
    )
    interpretation = (
        "adequately powered"
        if adequate
        else "underpowered under the locked 80 per cent joint-recovery rule"
    )
    columns = [
        "target_abs_residual_correlation",
        "median_recovered_rho",
        "p_component_recovery",
        "median_cv_rmse_fractional_improvement",
        "cv_component_recovery",
        "sign_component_recovery",
        "joint_primary_recovery",
    ]
    threshold_text = f"{threshold:.6g}" if math.isfinite(threshold) else "not reached"
    lines = [
        "# SPARC Group-Richness Injection-and-Recovery Calibration",
        "",
        f"Date: {RUN_DATE}",
        f"Completed: `{utc_now()}`",
        "",
        "## Bottom Line",
        "",
        f"At the earlier-effect anchor |r|=`{PRIOR_EFFECT_ANCHOR:.6f}`, the p component recovers in `{100.0 * float(anchor['p_component_recovery']):.2f}` per cent of simulations, the 5 per cent CV component in `{100.0 * float(anchor['cv_component_recovery']):.2f}` per cent, and the joint primary gate in `{100.0 * float(anchor['joint_primary_recovery']):.2f}` per cent. The current group-richness design is **{interpretation}** for an effect of that standardised magnitude.",
        "",
        f"The smallest tested effect with at least 80 per cent joint primary recovery is `{threshold_text}`.",
        "",
        "## Simulation Grid",
        "",
        markdown_table(summary, columns),
        "",
        "## Calibration",
        "",
        f"The fixed two-sided p<=0.01 critical value is |rho|>=`{float(diagnostics['critical_absolute_spearman_rho_for_p_le_0p01']):.6g}`, calibrated with `{NULL_PERMUTATIONS}` null rank permutations. Every grid point uses `{REPLICATES_PER_EFFECT}` empirical-residual injections.",
        "",
        "The simulation conditions on the observed 65-galaxy design and empirical residual distribution. It does not include group-catalogue misclassification, alternative host populations, cosmic variance or uncertainty in transferring an effect from the 2MRS composite to group richness.",
        "",
        "## Claim Boundary",
        "",
        "This study calibrates sensitivity only. It cannot turn the observed null into evidence for FR, antimatter or an environmental effect, and the five reserved galaxies remain sealed.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def export_bundle(out_dir: Path, export_dir: Path, programme: Path) -> None:
    export_dir.mkdir(parents=True, exist_ok=True)
    for path in out_dir.iterdir():
        if path.is_file():
            shutil.copy2(path, export_dir / path.name)
    rows: list[dict[str, Any]] = []
    for path in sorted(export_dir.iterdir()):
        if not path.is_file() or path.name == "manifest.csv":
            continue
        rows.append(
            {
                "path": path.name,
                "sha256": sha256(path),
                "bytes": path.stat().st_size,
                "role": "committed power-calibration product",
                "tracked_in_git": True,
            }
        )
    for path, role in [
        (programme, "locked power-calibration programme"),
        (PRIMARY_FRAME, "locked 65-galaxy primary frame"),
        (GROUP_PROGRAMME, "locked group-replication programme"),
        (GROUP_PREREGISTRATION, "locked group-replication preregistration"),
        (GROUP_GATE, "locked failed group-replication gate"),
        (GROUP_SUMMARY, "locked group-replication summary"),
        (PRIOR_ENVIRONMENT_SUMMARY, "earlier 2MRS environment summary"),
    ]:
        rows.append(
            {
                "path": str(path.relative_to(ROOT)).replace("\\", "/"),
                "sha256": sha256(path),
                "bytes": path.stat().st_size,
                "role": role,
                "tracked_in_git": True,
            }
        )
    pd.DataFrame(rows).to_csv(export_dir / "manifest.csv", index=False)


def preregister_phase(out_dir: Path, export_dir: Path, programme: Path) -> None:
    for path in [
        PRIMARY_FRAME,
        GROUP_PREREGISTRATION,
        GROUP_GATE,
        GROUP_SUMMARY,
        GROUP_PROGRAMME,
        PRIOR_ENVIRONMENT_SUMMARY,
    ]:
        if not path.exists():
            raise FileNotFoundError(path)
    out_dir.mkdir(parents=True, exist_ok=True)
    config = preregistration_config(programme)
    json_path = out_dir / "sparc_fr_group_power_preregistration.json"
    markdown_path = out_dir / "sparc_fr_group_power_preregistration.md"
    json_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    write_preregistration(markdown_path, config)
    export_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(json_path, export_dir / json_path.name)
    shutil.copy2(markdown_path, export_dir / markdown_path.name)
    print(f"Saved preregistration: {markdown_path}")
    print(f"Programme SHA-256: {config['programme_sha256']}")
    print(f"Replicates: {len(TARGET_ABS_RESIDUAL_CORRELATIONS)} x {REPLICATES_PER_EFFECT}")
    print("Reserved replication galaxies present: no")


def analyse_phase(out_dir: Path, export_dir: Path, programme: Path) -> None:
    prereg_path = out_dir / "sparc_fr_group_power_preregistration.json"
    if not prereg_path.exists():
        raise FileNotFoundError("Run and commit --phase preregister first")
    config = json.loads(prereg_path.read_text(encoding="utf-8"))
    if sha256(programme) != config["programme_sha256"]:
        raise RuntimeError("Programme hash no longer matches preregistration")
    locked_inputs = {
        "primary_frame": PRIMARY_FRAME,
        "group_preregistration": GROUP_PREREGISTRATION,
        "group_gate": GROUP_GATE,
        "group_summary": GROUP_SUMMARY,
        "group_programme": GROUP_PROGRAMME,
        "prior_environment_summary": PRIOR_ENVIRONMENT_SUMMARY,
    }
    for key, path in locked_inputs.items():
        if sha256(path) != config["input_hashes"][key]:
            raise RuntimeError(f"Input hash changed after preregistration: {key}")
    frame = pd.read_csv(PRIMARY_FRAME)
    excluded = set(frame["Galaxy"].astype(str)) & (RESERVED_REPLICATION | CANDIDATES)
    if excluded:
        raise AssertionError(f"Excluded galaxies entered simulation: {sorted(excluded)}")
    draws, summary, diagnostics = run_simulations(frame)
    draws.to_csv(
        out_dir / "sparc_fr_group_power_simulation_draws.csv.gz",
        index=False,
        compression="gzip",
    )
    summary.to_csv(out_dir / "sparc_fr_group_power_summary.csv", index=False)
    (out_dir / "sparc_fr_group_power_diagnostics.json").write_text(
        json.dumps(diagnostics, indent=2) + "\n", encoding="utf-8"
    )
    make_plot(summary, out_dir / "sparc_fr_group_power_curve.png")
    write_report(out_dir / "sparc_fr_group_power_report.md", summary, diagnostics)
    export_bundle(out_dir, export_dir, programme)
    anchor = summary[
        np.isclose(summary["target_abs_residual_correlation"], PRIOR_EFFECT_ANCHOR)
    ].iloc[0]
    print(f"Saved report: {out_dir / 'sparc_fr_group_power_report.md'}")
    print(f"Anchor joint recovery: {float(anchor['joint_primary_recovery']):.6f}")
    print("Reserved replication galaxies opened: no")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", choices=["preregister", "analyse"], required=True)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--export-dir", type=Path, default=DEFAULT_EXPORT)
    args = parser.parse_args()
    programme = Path(__file__).resolve()
    if args.phase == "preregister":
        preregister_phase(args.out_dir, args.export_dir, programme)
    else:
        analyse_phase(args.out_dir, args.export_dir, programme)


if __name__ == "__main__":
    main()
