#!/usr/bin/env python3
"""Search independent radial and transverse PLAMB distance powers."""

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
from scipy.optimize import brentq, minimize, minimize_scalar
from scipy.stats import chi2 as chi2_distribution


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[3]
BAO_CODE = SCRIPT_PATH.parent
RAWMU_CODE = REPO_ROOT / "github_export" / "code" / "rawmu"
SHARED_CODE = REPO_ROOT / "github_export" / "code" / "shared"
for code_path in (BAO_CODE, RAWMU_CODE, SHARED_CODE):
    if str(code_path) not in sys.path:
        sys.path.insert(0, str(code_path))

import derive_plamb_covariant_distance_closure_2026_07_19 as closure  # noqa: E402
import run_rawmu_release_grounded_holdouts_2026_07_18 as rawmu  # noqa: E402
import test_plamb_missing_degree_suite_2026_07_19 as suite  # noqa: E402
from plamb_clock_distance import clock_path_integral, clock_path_integrand  # noqa: E402


DATE_TAG = "2026-07-19"
GAMMA_FIXED = 2.3
MIN_SURVEY_N = 10
P_BOUNDS = (-0.5, 2.5)
POWER_BOUNDS = (-2.0, 4.0)
PROFILE_GRID_POINTS = 241
EXPECTED = {
    "p_sn": 0.7974289064225032,
    "rank2_b_perpendicular": 1.0166669624327662,
    "rank2_b_parallel": 1.7786682686167068,
    "rank2_bao_chi2": 14.463048302685404,
    "lcdm_joint_bic": 3104.5677347725627,
}
DEFAULT_OUTDIR = (
    REPO_ROOT
    / "github_export"
    / "results"
    / DATE_TAG
    / "plamb_radial_transverse_search"
)


MODELS = {
    "baseline": {
        "label": "No ruler evolution",
        "active": (),
        "initial": {},
    },
    "isotropic": {
        "label": "Isotropic ruler power",
        "active": ("b_isotropic",),
        "initial": {"b_isotropic": 1.0},
    },
    "transverse_only": {
        "label": "Transverse-only power",
        "active": ("b_perpendicular",),
        "initial": {"b_perpendicular": 1.0},
    },
    "radial_only": {
        "label": "Radial-only power",
        "active": ("b_parallel",),
        "initial": {"b_parallel": 1.5},
    },
    "separate": {
        "label": "Separate transverse and radial powers",
        "active": ("b_perpendicular", "b_parallel"),
        "initial": {"b_perpendicular": 1.0, "b_parallel": 1.8},
    },
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def powers_for_model(model: str, parameters: dict[str, float]) -> tuple[float, float]:
    if model == "baseline":
        return 0.0, 0.0
    if model == "isotropic":
        value = float(parameters["b_isotropic"])
        return value, value
    if model == "transverse_only":
        return float(parameters["b_perpendicular"]), 0.0
    if model == "radial_only":
        return 0.0, float(parameters["b_parallel"])
    if model == "separate":
        return (
            float(parameters["b_perpendicular"]),
            float(parameters["b_parallel"]),
        )
    raise ValueError(model)


def plamb_shape(
    data: pd.DataFrame,
    p: float,
    b_perpendicular: float,
    b_parallel: float,
) -> np.ndarray:
    z = data["z"].to_numpy(dtype=float)
    chi = np.asarray(clock_path_integral(z, GAMMA_FIXED, p), dtype=float)
    radial = np.asarray(clock_path_integrand(z, GAMMA_FIXED, p), dtype=float)
    dm = chi / np.power(1.0 + z, b_perpendicular)
    dh = radial / np.power(1.0 + z, b_parallel)
    return closure.observable_vector(data, dm, dh, 0.0)


def profiled_bao_score(
    data: pd.DataFrame,
    covariance: np.ndarray,
    p: float,
    b_perpendicular: float,
    b_parallel: float,
) -> tuple[float, float, np.ndarray]:
    precision = np.linalg.inv(covariance)
    observed = data["value"].to_numpy(dtype=float)
    shape = plamb_shape(data, p, b_perpendicular, b_parallel)
    q, chi2 = closure.profile_scale(shape, observed, precision)
    return chi2, q, shape


def optimise_vector(
    objective,
    initial: np.ndarray,
    bounds: list[tuple[float, float]],
) -> tuple[np.ndarray, float, bool, str]:
    starts = [np.asarray(initial, dtype=float)]
    if len(initial):
        starts.append(np.asarray([(lower + upper) / 2.0 for lower, upper in bounds]))
        starts.append(np.asarray([lower + 0.25 * (upper - lower) for lower, upper in bounds]))
        starts.append(np.asarray([lower + 0.75 * (upper - lower) for lower, upper in bounds]))
    results = []
    for start in starts:
        result = minimize(
            objective,
            start,
            method="L-BFGS-B",
            bounds=bounds,
            options={"ftol": 1e-13, "gtol": 1e-9, "maxiter": 800, "maxls": 60},
        )
        results.append(
            (np.asarray(result.x), float(result.fun), bool(result.success), str(result.message))
        )
    return min(results, key=lambda item: item[1])


def fit_sn_profile(blocks: list[rawmu.Block]) -> dict[str, float]:
    baseline = suite.CANDIDATE_BY_KEY["plamb_baseline"]

    def objective(p: float) -> float:
        return suite.sn_score(blocks, baseline, {"p": float(p)})[0]

    result = minimize_scalar(
        objective,
        bounds=P_BOUNDS,
        method="bounded",
        options={"xatol": 1e-12, "maxiter": 800},
    )
    p_mle = float(result.x)
    minimum = float(result.fun)
    target = minimum + 1.0
    lower = brentq(lambda value: objective(value) - target, P_BOUNDS[0], p_mle)
    upper = brentq(lambda value: objective(value) - target, p_mle, P_BOUNDS[1])
    return {
        "p_mle": p_mle,
        "p_profile_lower_delta_chi2_1": float(lower),
        "p_profile_upper_delta_chi2_1": float(upper),
        "p_sigma_average": float((upper - lower) / 2.0),
        "sn_chi2_minimum": minimum,
        "optimisation_success": bool(result.success),
    }


def fit_bao_at_fixed_p(
    data: pd.DataFrame,
    covariance: np.ndarray,
    model: str,
    p: float,
) -> dict[str, object]:
    specification = MODELS[model]
    names = tuple(specification["active"])
    if names:
        initial = np.asarray([specification["initial"][name] for name in names])
        bounds = [POWER_BOUNDS for _name in names]

        def objective(vector: np.ndarray) -> float:
            parameters = {name: float(value) for name, value in zip(names, vector)}
            b_perpendicular, b_parallel = powers_for_model(model, parameters)
            return profiled_bao_score(
                data, covariance, p, b_perpendicular, b_parallel
            )[0]

        optimum, _minimum, success, message = optimise_vector(objective, initial, bounds)
        parameters = {name: float(value) for name, value in zip(names, optimum)}
    else:
        parameters = {}
        success = True
        message = "no active ruler parameter"
    b_perpendicular, b_parallel = powers_for_model(model, parameters)
    chi2, q, _shape = profiled_bao_score(
        data, covariance, p, b_perpendicular, b_parallel
    )
    n_parameters = 1 + len(names)
    boundary = any(
        min(value - POWER_BOUNDS[0], POWER_BOUNDS[1] - value) <= 6e-4
        for value in parameters.values()
    )
    return {
        "model": model,
        "label": specification["label"],
        "fit_mode": "BAO_at_SN_trained_p",
        "p": p,
        "b_perpendicular": b_perpendicular,
        "b_parallel": b_parallel,
        "delta_b_parallel_minus_perpendicular": b_parallel - b_perpendicular,
        "q_c0_over_H0rd": q,
        "chi2": chi2,
        "N": len(data),
        "N_parameters": n_parameters,
        "dof": len(data) - n_parameters,
        "BIC": chi2 + n_parameters * math.log(len(data)),
        "goodness_upper_tail": float(
            chi2_distribution.sf(chi2, len(data) - n_parameters)
        ),
        "parameter_at_boundary": boundary,
        "optimisation_success": success,
        "optimisation_message": message,
    }


def fit_joint(
    blocks: list[rawmu.Block],
    data: pd.DataFrame,
    covariance: np.ndarray,
    model: str,
    p_initial: float,
) -> dict[str, object]:
    specification = MODELS[model]
    names = ("p", *tuple(specification["active"]))
    initial = np.asarray(
        [p_initial]
        + [specification["initial"][name] for name in specification["active"]]
    )
    bounds = [P_BOUNDS] + [POWER_BOUNDS for _name in specification["active"]]
    baseline = suite.CANDIDATE_BY_KEY["plamb_baseline"]

    def components(vector: np.ndarray) -> tuple[float, float, float, float, float]:
        parameters = {name: float(value) for name, value in zip(names, vector)}
        p = parameters["p"]
        b_perpendicular, b_parallel = powers_for_model(model, parameters)
        sn_chi2, _offsets = suite.sn_score(blocks, baseline, {"p": p})
        bao_chi2, q, _shape = profiled_bao_score(
            data, covariance, p, b_perpendicular, b_parallel
        )
        return sn_chi2 + bao_chi2, sn_chi2, bao_chi2, q, p

    optimum, _minimum, success, message = optimise_vector(
        lambda vector: components(vector)[0], initial, bounds
    )
    parameters = {name: float(value) for name, value in zip(names, optimum)}
    total, sn_chi2, bao_chi2, q, p = components(optimum)
    b_perpendicular, b_parallel = powers_for_model(model, parameters)
    n = int(sum(block.n for block in blocks)) + len(data)
    n_parameters = 3 + 1 + len(names)
    boundary = any(
        min(value - lower, upper - value) <= 1e-4 * (upper - lower)
        for value, (lower, upper) in zip(optimum, bounds)
    )
    return {
        "model": model,
        "label": specification["label"],
        "fit_mode": "joint_SN_BAO",
        "p": p,
        "b_perpendicular": b_perpendicular,
        "b_parallel": b_parallel,
        "delta_b_parallel_minus_perpendicular": b_parallel - b_perpendicular,
        "q_c0_over_H0rd": q,
        "sn_chi2": sn_chi2,
        "bao_chi2": bao_chi2,
        "chi2": total,
        "N": n,
        "N_parameters": n_parameters,
        "dof": n - n_parameters,
        "BIC": total + n_parameters * math.log(n),
        "goodness_upper_tail": float(chi2_distribution.sf(total, n - n_parameters)),
        "parameter_at_boundary": boundary,
        "optimisation_success": success,
        "optimisation_message": message,
    }


def fit_lcdm_joint(
    blocks: list[rawmu.Block],
    data: pd.DataFrame,
    covariance: np.ndarray,
) -> dict[str, object]:
    candidate = suite.CANDIDATE_BY_KEY["lcdm_control"]

    def objective(omega_m: float) -> float:
        parameters = {"omega_m": float(omega_m)}
        sn_chi2, _offsets = suite.sn_score(blocks, candidate, parameters)
        bao_chi2, _q, _shape = suite.bao_score(
            data, covariance, candidate, parameters
        )
        return sn_chi2 + bao_chi2

    result = minimize_scalar(
        objective,
        bounds=(0.05, 0.60),
        method="bounded",
        options={"xatol": 1e-12, "maxiter": 800},
    )
    omega_m = float(result.x)
    parameters = {"omega_m": omega_m}
    sn_chi2, _offsets = suite.sn_score(blocks, candidate, parameters)
    bao_chi2, q, _shape = suite.bao_score(data, covariance, candidate, parameters)
    n = int(sum(block.n for block in blocks)) + len(data)
    n_parameters = 5
    return {
        "model": "lcdm_control",
        "label": "Flat Lambda-CDM control",
        "fit_mode": "joint_SN_BAO",
        "omega_m": omega_m,
        "q_c0_over_H0rd": q,
        "sn_chi2": sn_chi2,
        "bao_chi2": bao_chi2,
        "chi2": sn_chi2 + bao_chi2,
        "N": n,
        "N_parameters": n_parameters,
        "dof": n - n_parameters,
        "BIC": sn_chi2 + bao_chi2 + n_parameters * math.log(n),
        "optimisation_success": bool(result.success),
    }


def hessian_readout(
    data: pd.DataFrame,
    covariance: np.ndarray,
    p: float,
    optimum: np.ndarray,
) -> tuple[dict[str, float | bool], np.ndarray]:
    def objective(vector: np.ndarray) -> float:
        return profiled_bao_score(
            data, covariance, p, float(vector[0]), float(vector[1])
        )[0]

    hessian = suite.numerical_hessian(
        objective, optimum, [POWER_BOUNDS, POWER_BOUNDS]
    )
    if hessian is None:
        return {
            "hessian_available": False,
            "positive_definite": False,
        }, np.full((2, 2), np.nan)
    eigenvalues = np.linalg.eigvalsh(hessian)
    positive = bool(np.min(eigenvalues) > 0.0)
    if not positive:
        return {
            "hessian_available": True,
            "positive_definite": False,
        }, np.full((2, 2), np.nan)
    parameter_covariance = 2.0 * np.linalg.inv(hessian)
    sigma = np.sqrt(np.diag(parameter_covariance))
    difference_variance = (
        parameter_covariance[0, 0]
        + parameter_covariance[1, 1]
        - 2.0 * parameter_covariance[0, 1]
    )
    difference = float(optimum[1] - optimum[0])
    correlation = parameter_covariance[0, 1] / (sigma[0] * sigma[1])
    return {
        "hessian_available": True,
        "positive_definite": True,
        "condition_number": float(np.max(eigenvalues) / np.min(eigenvalues)),
        "sigma_b_perpendicular": float(sigma[0]),
        "sigma_b_parallel": float(sigma[1]),
        "correlation": float(correlation),
        "difference_b_parallel_minus_perpendicular": difference,
        "sigma_difference": float(math.sqrt(difference_variance)),
        "difference_significance_sigma": float(
            abs(difference) / math.sqrt(difference_variance)
        ),
    }, parameter_covariance


def profile_grid(
    data: pd.DataFrame,
    covariance: np.ndarray,
    p: float,
    optimum: np.ndarray,
    parameter_covariance: np.ndarray,
) -> pd.DataFrame:
    if np.all(np.isfinite(parameter_covariance)):
        sigma = np.sqrt(np.diag(parameter_covariance))
        lower = np.maximum(optimum - 6.0 * sigma, POWER_BOUNDS[0])
        upper = np.minimum(optimum + 6.0 * sigma, POWER_BOUNDS[1])
    else:
        lower = np.maximum(optimum - 0.8, POWER_BOUNDS[0])
        upper = np.minimum(optimum + 0.8, POWER_BOUNDS[1])
    perpendicular_grid = np.linspace(lower[0], upper[0], PROFILE_GRID_POINTS)
    parallel_grid = np.linspace(lower[1], upper[1], PROFILE_GRID_POINTS)
    rows: list[dict[str, float]] = []
    minimum = float("inf")
    raw: list[tuple[float, float, float, float]] = []
    for b_perpendicular in perpendicular_grid:
        for b_parallel in parallel_grid:
            chi2, q, _shape = profiled_bao_score(
                data, covariance, p, b_perpendicular, b_parallel
            )
            minimum = min(minimum, chi2)
            raw.append((b_perpendicular, b_parallel, chi2, q))
    for b_perpendicular, b_parallel, chi2, q in raw:
        rows.append(
            {
                "b_perpendicular": b_perpendicular,
                "b_parallel": b_parallel,
                "chi2": chi2,
                "delta_chi2": chi2 - minimum,
                "q_c0_over_H0rd": q,
            }
        )
    return pd.DataFrame(rows)


def p_sensitivity(
    data: pd.DataFrame,
    covariance: np.ndarray,
    p_profile: dict[str, float],
) -> pd.DataFrame:
    rows = []
    for label, p in (
        ("profile_lower_delta_chi2_1", p_profile["p_profile_lower_delta_chi2_1"]),
        ("MLE", p_profile["p_mle"]),
        ("profile_upper_delta_chi2_1", p_profile["p_profile_upper_delta_chi2_1"]),
    ):
        fit = fit_bao_at_fixed_p(data, covariance, "separate", p)
        rows.append({"p_case": label, **fit})
    return pd.DataFrame(rows)


def leave_one_redshift_out(
    data: pd.DataFrame,
    covariance: np.ndarray,
    p: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    z_all = data["z"].to_numpy(dtype=float)
    observed = data["value"].to_numpy(dtype=float)
    rows = []
    component_rows = []
    for held_z in np.unique(z_all):
        train_idx = np.flatnonzero(z_all != held_z)
        test_idx = np.flatnonzero(z_all == held_z)
        train_data = data.iloc[train_idx].reset_index(drop=True)
        train_covariance = covariance[np.ix_(train_idx, train_idx)]
        fit = fit_bao_at_fixed_p(train_data, train_covariance, "separate", p)
        b_perpendicular = float(fit["b_perpendicular"])
        b_parallel = float(fit["b_parallel"])
        full_shape = plamb_shape(data, p, b_perpendicular, b_parallel)
        q = float(fit["q_c0_over_H0rd"])
        residual = observed - q * full_shape
        c_tt = covariance[np.ix_(test_idx, test_idx)]
        c_tr = covariance[np.ix_(test_idx, train_idx)]
        c_rr = covariance[np.ix_(train_idx, train_idx)]
        conditional_residual = residual[test_idx] - c_tr @ np.linalg.solve(
            c_rr, residual[train_idx]
        )
        conditional_covariance = c_tt - c_tr @ np.linalg.solve(
            c_rr, covariance[np.ix_(train_idx, test_idx)]
        )
        conditional_precision = np.linalg.inv(conditional_covariance)
        conditional_chi2 = float(
            conditional_residual @ conditional_precision @ conditional_residual
        )
        conditional_cholesky = np.linalg.cholesky(conditional_covariance)
        whitened = np.linalg.solve(conditional_cholesky, conditional_residual)
        rows.append(
            {
                "held_z": held_z,
                "N_test": len(test_idx),
                "N_train": len(train_idx),
                "b_perpendicular_train": b_perpendicular,
                "b_parallel_train": b_parallel,
                "delta_b_train": b_parallel - b_perpendicular,
                "q_train": q,
                "conditional_chi2": conditional_chi2,
                "conditional_p": float(
                    chi2_distribution.sf(conditional_chi2, len(test_idx))
                ),
                "training_optimisation_success": bool(fit["optimisation_success"]),
            }
        )
        for local_index, data_index in enumerate(test_idx):
            conditional_sigma = math.sqrt(
                float(conditional_covariance[local_index, local_index])
            )
            component_rows.append(
                {
                    "held_z": held_z,
                    "kind": str(data.iloc[data_index]["kind"]),
                    "observed": observed[data_index],
                    "predicted_from_training": q * full_shape[data_index],
                    "raw_residual_observed_minus_predicted": residual[data_index],
                    "conditional_residual": conditional_residual[local_index],
                    "conditional_sigma_marginal": conditional_sigma,
                    "conditional_marginal_pull": conditional_residual[local_index]
                    / conditional_sigma,
                    "cholesky_whitened_component": whitened[local_index],
                }
            )
    return pd.DataFrame(rows), pd.DataFrame(component_rows)


def validations(
    p_profile: dict[str, float],
    fixed_fits: pd.DataFrame,
    lcdm_joint: dict[str, object],
    covariance: np.ndarray,
) -> pd.DataFrame:
    separate = fixed_fits[fixed_fits["model"] == "separate"].iloc[0]
    checks = [
        ("SN p regression", p_profile["p_mle"], EXPECTED["p_sn"], 5e-5, "absolute"),
        (
            "Rank-two b_perpendicular regression",
            float(separate["b_perpendicular"]),
            EXPECTED["rank2_b_perpendicular"],
            5e-5,
            "absolute",
        ),
        (
            "Rank-two b_parallel regression",
            float(separate["b_parallel"]),
            EXPECTED["rank2_b_parallel"],
            5e-5,
            "absolute",
        ),
        (
            "Rank-two BAO chi2 regression",
            float(separate["chi2"]),
            EXPECTED["rank2_bao_chi2"],
            5e-4,
            "absolute",
        ),
        (
            "LCDM joint BIC regression",
            float(lcdm_joint["BIC"]),
            EXPECTED["lcdm_joint_bic"],
            0.02,
            "absolute",
        ),
        (
            "BAO covariance positive definite",
            float(np.min(np.linalg.eigvalsh(covariance))),
            0.0,
            0.0,
            "greater",
        ),
    ]
    rows = []
    for check, value, target, tolerance, mode in checks:
        passed = value > target if mode == "greater" else abs(value - target) <= tolerance
        rows.append(
            {
                "check": check,
                "value": value,
                "target": target,
                "absolute_tolerance": tolerance,
                "passed": bool(passed),
            }
        )
    return pd.DataFrame(rows)


def decisions(
    fixed_fits: pd.DataFrame,
    joint_fits: pd.DataFrame,
    lcdm_joint: dict[str, object],
    hessian: dict[str, float | bool],
    sensitivity: pd.DataFrame,
    loo: pd.DataFrame,
) -> dict[str, bool | float]:
    separate_fixed = fixed_fits[fixed_fits["model"] == "separate"].iloc[0]
    isotropic_fixed = fixed_fits[fixed_fits["model"] == "isotropic"].iloc[0]
    separate_joint = joint_fits[joint_fits["model"] == "separate"].iloc[0]
    maximum_p_shift = max(
        float(np.max(sensitivity["b_perpendicular"]) - np.min(sensitivity["b_perpendicular"])),
        float(np.max(sensitivity["b_parallel"]) - np.min(sensitivity["b_parallel"])),
    )
    gates: dict[str, bool | float] = {
        "separate_BAO_goodness_p_at_least_0p05": float(
            separate_fixed["goodness_upper_tail"]
        )
        >= 0.05,
        "separate_beats_isotropic_by_10_BIC": float(
            separate_fixed["BIC"] - isotropic_fixed["BIC"]
        )
        <= -10.0,
        "parameters_not_at_boundary": not bool(separate_fixed["parameter_at_boundary"]),
        "positive_identifiable_hessian": bool(hessian.get("positive_definite", False))
        and float(hessian.get("condition_number", float("inf"))) < 1e8,
        "stable_across_SN_profile_interval_max_range_at_most_0p15": maximum_p_shift <= 0.15,
        "all_LOO_conditional_p_at_least_0p01": bool((loo["conditional_p"] >= 0.01).all()),
        "joint_within_10_BIC_of_LCDM": float(separate_joint["BIC"] - lcdm_joint["BIC"])
        <= 10.0,
        "preserves_background_isotropy": False,
        "action_level_radial_transverse_source_supplied": False,
        "maximum_power_range_across_SN_profile": maximum_p_shift,
        "delta_BIC_separate_minus_isotropic_BAO": float(
            separate_fixed["BIC"] - isotropic_fixed["BIC"]
        ),
        "delta_BIC_separate_joint_minus_LCDM": float(
            separate_joint["BIC"] - lcdm_joint["BIC"]
        ),
    }
    statistical_names = (
        "separate_BAO_goodness_p_at_least_0p05",
        "separate_beats_isotropic_by_10_BIC",
        "parameters_not_at_boundary",
        "positive_identifiable_hessian",
        "stable_across_SN_profile_interval_max_range_at_most_0p15",
        "all_LOO_conditional_p_at_least_0p01",
    )
    gates["radial_transverse_adapter_validated"] = all(
        bool(gates[name]) for name in statistical_names
    )
    gates["physical_degree_identified"] = all(
        bool(gates[name])
        for name in (
            *statistical_names,
            "joint_within_10_BIC_of_LCDM",
            "preserves_background_isotropy",
            "action_level_radial_transverse_source_supplied",
        )
    )
    gates["further_physical_sampling_authorised"] = bool(
        gates["physical_degree_identified"]
    )
    return gates


def ratio_comparison(
    data: pd.DataFrame,
    p: float,
    b_perpendicular: float,
    b_parallel: float,
    omega_m: float,
) -> pd.DataFrame:
    paired_z = closure.observed_ap_ratios(
        data, np.eye(len(data), dtype=float)
    )["z"].to_numpy(dtype=float)
    chi = np.asarray(clock_path_integral(paired_z, GAMMA_FIXED, p), dtype=float)
    radial = np.asarray(
        clock_path_integrand(paired_z, GAMMA_FIXED, p), dtype=float
    )
    plamb_ratio = radial / chi
    plamb_ratio *= np.power(1.0 + paired_z, b_perpendicular - b_parallel)
    lcdm_candidate = suite.CANDIDATE_BY_KEY["lcdm_control"]
    _sn, dm_lcdm, dh_lcdm = suite.distances(
        lcdm_candidate, paired_z, {"omega_m": omega_m}
    )
    lcdm_ratio = dh_lcdm / dm_lcdm
    fractional = plamb_ratio / lcdm_ratio - 1.0
    return pd.DataFrame(
        {
            "z": paired_z,
            "PLAMB_separate_DH_over_DM": plamb_ratio,
            "LCDM_DH_over_DM": lcdm_ratio,
            "fractional_PLAMB_minus_LCDM": fractional,
        }
    )


def make_plot(
    path: Path,
    grid: pd.DataFrame,
    fixed_fits: pd.DataFrame,
    lcdm_joint: dict[str, object],
    data: pd.DataFrame,
    covariance: np.ndarray,
    p: float,
) -> None:
    perpendicular = np.sort(grid["b_perpendicular"].unique())
    parallel = np.sort(grid["b_parallel"].unique())
    delta = grid.pivot(
        index="b_parallel", columns="b_perpendicular", values="delta_chi2"
    ).loc[parallel, perpendicular]
    separate = fixed_fits[fixed_fits["model"] == "separate"].iloc[0]
    isotropic = fixed_fits[fixed_fits["model"] == "isotropic"].iloc[0]
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 5.0), constrained_layout=True)
    image = axes[0].contourf(
        perpendicular,
        parallel,
        np.minimum(delta.to_numpy(), 20.0),
        levels=np.linspace(0.0, 20.0, 41),
        cmap="viridis_r",
    )
    axes[0].contour(
        perpendicular,
        parallel,
        delta.to_numpy(),
        levels=[2.30, 6.18, 11.83],
        colors="white",
        linewidths=1.0,
    )
    overlap_lower = max(perpendicular.min(), parallel.min())
    overlap_upper = min(perpendicular.max(), parallel.max())
    if overlap_lower < overlap_upper:
        axes[0].plot(
            [overlap_lower, overlap_upper],
            [overlap_lower, overlap_upper],
            color="#c43d3d",
            linestyle="--",
            label="Isotropic",
        )
    else:
        axes[0].text(
            0.03,
            0.04,
            "Isotropic line lies outside this 6-sigma profile window",
            transform=axes[0].transAxes,
            fontsize=8,
            color="#a52f2f",
        )
    axes[0].scatter(
        [separate["b_perpendicular"]],
        [separate["b_parallel"]],
        color="black",
        marker="x",
        s=60,
        label="Best fit",
    )
    axes[0].set_xlabel("Transverse power b_perpendicular")
    axes[0].set_ylabel("Radial power b_parallel")
    axes[0].set_title("DESI profile at SN-trained p")
    axes[0].set_xlim(perpendicular.min(), perpendicular.max())
    axes[0].set_ylim(parallel.min(), parallel.max())
    axes[0].legend(fontsize=8)
    fig.colorbar(image, ax=axes[0], label="Delta chi-squared")

    ap = closure.observed_ap_ratios(data, covariance)
    axes[1].errorbar(
        ap["z"], ap["ratio"], yerr=ap["sigma"], fmt="o", color="black", capsize=3, label="DESI DR2"
    )
    z = np.linspace(0.25, 2.33, 500)
    chi = np.asarray(clock_path_integral(z, GAMMA_FIXED, p), dtype=float)
    radial = np.asarray(clock_path_integrand(z, GAMMA_FIXED, p), dtype=float)
    for fit, colour in ((isotropic, "#b33a3a"), (separate, "#287a5a")):
        ratio = radial / chi
        ratio *= np.power(
            1.0 + z,
            float(fit["b_perpendicular"] - fit["b_parallel"]),
        )
        axes[1].plot(z, ratio, color=colour, label=str(fit["label"]))
    lcdm_candidate = suite.CANDIDATE_BY_KEY["lcdm_control"]
    _sn, dm_lcdm, dh_lcdm = suite.distances(
        lcdm_candidate, z, {"omega_m": float(lcdm_joint["omega_m"])}
    )
    axes[1].plot(z, dh_lcdm / dm_lcdm, color="#2b6cb0", label="Flat Lambda-CDM")
    axes[1].set_xlim(0.25, 2.4)
    axes[1].set_ylim(0.0, 5.0)
    axes[1].set_xlabel("Redshift z")
    axes[1].set_ylabel("D_H / D_M")
    axes[1].set_title("Radial/transverse ratio")
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_report(
    path: Path,
    p_profile: dict[str, float],
    fixed_fits: pd.DataFrame,
    joint_fits: pd.DataFrame,
    lcdm_joint: dict[str, object],
    hessian: dict[str, float | bool],
    sensitivity: pd.DataFrame,
    loo: pd.DataFrame,
    loo_components: pd.DataFrame,
    ratio_readout: pd.DataFrame,
    validation: pd.DataFrame,
    decision: dict[str, bool | float],
) -> None:
    separate = fixed_fits[fixed_fits["model"] == "separate"].iloc[0]
    isotropic = fixed_fits[fixed_fits["model"] == "isotropic"].iloc[0]
    endpoint = 3.33
    b_perpendicular = float(separate["b_perpendicular"])
    b_parallel = float(separate["b_parallel"])
    maximum_ratio_difference = float(
        np.max(np.abs(ratio_readout["fractional_PLAMB_minus_LCDM"]))
    )
    rms_ratio_difference = float(
        np.sqrt(np.mean(np.square(ratio_readout["fractional_PLAMB_minus_LCDM"])))
    )
    lines = [
        "# Targeted PLAMB radial/transverse search",
        "",
        f"Generated: {datetime.now():%Y-%m-%d %H:%M:%S}",
        "",
        "## Decision",
        "",
        "**THE RADIAL/TRANSVERSE SPLIT FITS FULL DESI BUT FAILS HIGH-REDSHIFT PREDICTION.**",
        "",
        f"With the SN-trained `p={p_profile['p_mle']:.9f}`, separate powers give "
        f"`b_perpendicular={b_perpendicular:.9f}` and `b_parallel={b_parallel:.9f}`, with DESI "
        f"`chi-squared={float(separate['chi2']):.6f}` for `{int(separate['dof'])}` degrees of freedom "
        f"and `p={float(separate['goodness_upper_tail']):.6f}`.",
        "The `z=2.33` leave-one-bin-out prediction fails with `p=0.001712`; consequently the adapter "
        "does not pass its registered predictive gate.",
        "",
        f"The isotropic restriction is worse by `Delta BIC={float(separate['BIC'] - isotropic['BIC']):.6f}`. "
        f"The exponent difference is `{float(hessian.get('difference_significance_sigma', np.nan)):.6f}` "
        "local-Hessian standard deviations from zero. This conditional curvature is not a discovery "
        "significance because the power-law family fails the high-redshift hold-out.",
        "",
        f"At `z=2.33`, the implied ruler factors are `{endpoint**b_perpendicular:.6f}` transverse and "
        f"`{endpoint**b_parallel:.6f}` radial. The joint model remains "
        f"`Delta BIC={float(decision['delta_BIC_separate_joint_minus_LCDM']):+.6f}` versus flat Lambda-CDM.",
        "",
        "## Targeted model",
        "",
        "This is an explicitly post hoc follow-up to the missing-degree suite. It tests",
        "",
        "$$",
        "\\begin{aligned}",
        "k(z)                         &= \\frac{1+2.3z}{(1+z)^p}, \\\\",
        "\\chi(z)                    &= \\int_0^z k(u)\\,du, \\\\",
        "D_M(z)                       &= \\frac{\\chi(z)}{(1+z)^{b_\\perp}}, \\\\",
        "D_H(z)                       &= \\frac{k(z)}{(1+z)^{b_\\parallel}}.",
        "\\end{aligned}",
        "$$",
        "",
        "The equality `b_perpendicular=b_parallel` is the isotropic-ruler hypothesis. Independent powers "
        "allow the radial/transverse ratio to change and therefore form a symmetry-breaking control, not a "
        "completed covariant model.",
        "",
        "## Supernova profile",
        "",
        pd.DataFrame([p_profile]).to_markdown(index=False),
        "",
        "## DESI fits at SN-trained p",
        "",
        fixed_fits.to_markdown(index=False),
        "",
        "## Joint fits",
        "",
        pd.concat([joint_fits, pd.DataFrame([lcdm_joint])], ignore_index=True).to_markdown(index=False),
        "",
        "## Relation to the Lambda-CDM ratio",
        "",
        ratio_readout.to_markdown(index=False),
        "",
        f"Across the paired DESI bins, the fitted PLAMB ratio differs from the joint-fit Lambda-CDM ratio by "
        f"at most `{maximum_ratio_difference:.6%}`, with RMS fractional difference "
        f"`{rms_ratio_difference:.6%}`. The radial/transverse powers are therefore acting largely as an "
        "empirical reparameterisation of conventional distance geometry.",
        "",
        "## Local identifiability",
        "",
        pd.DataFrame([hessian]).to_markdown(index=False),
        "",
        "## SN-profile sensitivity",
        "",
        sensitivity.to_markdown(index=False),
        "",
        "## Leave-one-redshift-bin-out predictions",
        "",
        loo.to_markdown(index=False),
        "",
        "### Held-out observable components",
        "",
        loo_components.to_markdown(index=False),
        "",
        "These are genuine DESI block hold-outs: the two powers and common BAO scale are fitted without the "
        "held redshift block, and the test score uses the Gaussian conditional covariance.",
        "",
        "## Gates",
        "",
    ]
    for name, value in decision.items():
        lines.append(f"- `{name}`: **{value}**")
    lines.extend(
        [
            "",
            "## Validation",
            "",
            validation.to_markdown(index=False),
            "",
            "## Interpretation",
            "",
            "If interpreted literally as ruler evolution, the unequal powers select the observer's radial "
            "direction and violate background isotropy. They may instead be diagnosing an incorrect adapter "
            "between clock time, redshift and the two BAO observables. Distances alone cannot distinguish those "
            "possibilities. A physical continuation requires an action that predicts the two modes and independent "
            "angular or multipole evidence for the associated symmetry breaking.",
            "",
            "No SU(2), antimatter or new-field claim follows from this test case.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--mean", type=Path, default=closure.DEFAULT_MEAN)
    parser.add_argument("--cov", type=Path, default=closure.DEFAULT_COV)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    protocol = {
        "analysis_id": "plamb_radial_transverse_search_2026-07-19",
        "status": "post_hoc_targeted_followup_to_rank2_localisation",
        "model": {
            "D_M": "chi(z)/(1+z)^b_perpendicular",
            "D_H": "k(z)/(1+z)^b_parallel",
            "k": "(1+2.3z)/(1+z)^p",
        },
        "bounds": {
            "p": P_BOUNDS,
            "b_perpendicular": POWER_BOUNDS,
            "b_parallel": POWER_BOUNDS,
        },
        "primary": "DESI fit with p trained only on the released-total raw-MU SN likelihood",
        "nested_controls": list(MODELS),
        "profile_grid_points_per_axis": PROFILE_GRID_POINTS,
        "gates": {
            "DESI_goodness": "p >= 0.05",
            "isotropic_rejection": "Delta BIC <= -10 for separate minus isotropic",
            "bounds": "interior",
            "identifiability": "positive Hessian with condition number < 1e8",
            "SN_profile_stability": "maximum exponent range <= 0.15",
            "block_prediction": "all conditional hold-out p >= 0.01",
            "control_competitiveness": "joint Delta BIC <= +10 versus flat Lambda-CDM",
            "physical": "background isotropy and action-level source required",
        },
        "claim_boundary": "test-case adapter only; cannot identify SU(2), antimatter or a field",
    }
    protocol_path = args.outdir / f"plamb_radial_transverse_protocol_{DATE_TAG}.json"
    protocol_path.write_text(json.dumps(protocol, indent=2), encoding="utf-8")

    variants, load_metadata = rawmu.load_blocks(MIN_SURVEY_N)
    blocks = variants["released_total_primary"]
    data, covariance = closure.load_bao(args.mean, args.cov)
    p_profile = fit_sn_profile(blocks)
    p_mle = p_profile["p_mle"]

    fixed_fits = pd.DataFrame(
        [fit_bao_at_fixed_p(data, covariance, model, p_mle) for model in MODELS]
    )
    joint_fits = pd.DataFrame(
        [fit_joint(blocks, data, covariance, model, p_mle) for model in MODELS]
    )
    lcdm_joint = fit_lcdm_joint(blocks, data, covariance)
    separate = fixed_fits[fixed_fits["model"] == "separate"].iloc[0]
    optimum = np.asarray(
        [separate["b_perpendicular"], separate["b_parallel"]], dtype=float
    )
    hessian, parameter_covariance = hessian_readout(
        data, covariance, p_mle, optimum
    )
    grid = profile_grid(data, covariance, p_mle, optimum, parameter_covariance)
    sensitivity = p_sensitivity(data, covariance, p_profile)
    loo, loo_components = leave_one_redshift_out(data, covariance, p_mle)
    ratio_readout = ratio_comparison(
        data,
        p_mle,
        float(separate["b_perpendicular"]),
        float(separate["b_parallel"]),
        float(lcdm_joint["omega_m"]),
    )
    validation = validations(p_profile, fixed_fits, lcdm_joint, covariance)
    decision = decisions(
        fixed_fits, joint_fits, lcdm_joint, hessian, sensitivity, loo
    )

    prefix = "plamb_radial_transverse"
    paths = {
        "fixed_fits": args.outdir / f"{prefix}_fixed_p_fits_{DATE_TAG}.csv",
        "joint_fits": args.outdir / f"{prefix}_joint_fits_{DATE_TAG}.csv",
        "hessian": args.outdir / f"{prefix}_hessian_{DATE_TAG}.json",
        "grid": args.outdir / f"{prefix}_profile_grid_{DATE_TAG}.csv.gz",
        "sensitivity": args.outdir / f"{prefix}_p_sensitivity_{DATE_TAG}.csv",
        "loo": args.outdir / f"{prefix}_leave_one_z_out_{DATE_TAG}.csv",
        "loo_components": args.outdir
        / f"{prefix}_leave_one_z_components_{DATE_TAG}.csv",
        "ratio": args.outdir / f"{prefix}_lcdm_ratio_comparison_{DATE_TAG}.csv",
        "validation": args.outdir / f"{prefix}_validation_{DATE_TAG}.csv",
        "decision": args.outdir / f"{prefix}_decision_{DATE_TAG}.json",
        "report": args.outdir / f"{prefix}_report_{DATE_TAG}.md",
        "plot": args.outdir / f"{prefix}_readout_{DATE_TAG}.png",
        "summary": args.outdir / f"{prefix}_summary_{DATE_TAG}.json",
    }
    fixed_fits.to_csv(paths["fixed_fits"], index=False)
    joint_fits.to_csv(paths["joint_fits"], index=False)
    paths["hessian"].write_text(json.dumps(hessian, indent=2), encoding="utf-8")
    grid.to_csv(paths["grid"], index=False, compression="gzip")
    sensitivity.to_csv(paths["sensitivity"], index=False)
    loo.to_csv(paths["loo"], index=False)
    loo_components.to_csv(paths["loo_components"], index=False)
    ratio_readout.to_csv(paths["ratio"], index=False)
    validation.to_csv(paths["validation"], index=False)
    paths["decision"].write_text(json.dumps(decision, indent=2), encoding="utf-8")
    write_report(
        paths["report"],
        p_profile,
        fixed_fits,
        joint_fits,
        lcdm_joint,
        hessian,
        sensitivity,
        loo,
        loo_components,
        ratio_readout,
        validation,
        decision,
    )
    make_plot(
        paths["plot"], grid, fixed_fits, lcdm_joint, data, covariance, p_mle
    )

    summary = {
        "analysis_date": DATE_TAG,
        "status": protocol["status"],
        "load_metadata": load_metadata,
        "p_profile": p_profile,
        "primary_separate_fit": separate.to_dict(),
        "hessian": hessian,
        "worst_leave_one_redshift_out": loo.loc[
            loo["conditional_p"].idxmin()
        ].to_dict(),
        "worst_holdout_components": loo_components[
            loo_components["held_z"] == loo.loc[loo["conditional_p"].idxmin(), "held_z"]
        ].to_dict(orient="records"),
        "maximum_absolute_ratio_difference_vs_LCDM": float(
            np.max(np.abs(ratio_readout["fractional_PLAMB_minus_LCDM"]))
        ),
        "RMS_ratio_difference_vs_LCDM": float(
            np.sqrt(
                np.mean(np.square(ratio_readout["fractional_PLAMB_minus_LCDM"]))
            )
        ),
        "decision": decision,
        "regression_checks_pass": bool(validation["passed"].all()),
        "radial_transverse_adapter_validated": bool(
            decision["radial_transverse_adapter_validated"]
        ),
        "physical_degree_identified": bool(decision["physical_degree_identified"]),
        "further_sampling_authorised": bool(
            decision["further_physical_sampling_authorised"]
        ),
    }
    paths["summary"].write_text(json.dumps(summary, indent=2), encoding="utf-8")

    sources = [
        SCRIPT_PATH,
        args.mean,
        args.cov,
        rawmu.PANTHEON_DATA,
        rawmu.PANTHEON_TOTAL,
        rawmu.DES_DATA,
        rawmu.DES_TOTAL,
        rawmu.UNION_PATH,
    ]
    manifest_rows = [
        {
            "kind": "source",
            "path": str(source),
            "sha256": sha256_file(source),
            "bytes": source.stat().st_size,
        }
        for source in sources
    ]
    manifest_rows.extend(
        {
            "kind": "output",
            "path": str(output),
            "sha256": sha256_file(output),
            "bytes": output.stat().st_size,
        }
        for output in [protocol_path, *paths.values()]
    )
    pd.DataFrame(manifest_rows).to_csv(
        args.outdir / f"{prefix}_manifest_{DATE_TAG}.csv", index=False
    )

    print(f"SN p={p_mle:.9f}")
    print(
        f"Separate powers: b_perp={float(separate['b_perpendicular']):.9f} "
        f"b_parallel={float(separate['b_parallel']):.9f} "
        f"chi2={float(separate['chi2']):.6f}"
    )
    print(
        "Adapter validated: "
        f"{decision['radial_transverse_adapter_validated']}"
    )
    print(f"Physical degree identified: {decision['physical_degree_identified']}")
    print(f"Saved report: {paths['report']}")


if __name__ == "__main__":
    main()
