#!/usr/bin/env python3
"""Test observable proxies for a missing PLAMB cosmological degree of freedom."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize
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
from diagnose_pantheon_rawmu_fr import C_KMS  # noqa: E402
from plamb_clock_distance import clock_path_integral, clock_path_integrand  # noqa: E402


DATE_TAG = "2026-07-19"
GAMMA_FIXED = 2.3
H0_FIXED = 67.5
MIN_SURVEY_N = 10
EXPECTED = {
    "plamb_sn_p": 0.7974289064225032,
    "plamb_sn_chi2": 3047.674475661129,
    "lcdm_sn_omega_m": 0.33067964339262307,
    "lcdm_sn_chi2": 3047.6625171859296,
    "lcdm_bao_omega_m": 0.29746182078985123,
    "lcdm_bao_chi2": 10.271041425445219,
    "fixed_plamb_bao_chi2": 33458.051670260524,
}
DEFAULT_OUTDIR = (
    REPO_ROOT
    / "github_export"
    / "results"
    / DATE_TAG
    / "plamb_missing_degree_suite"
)


@dataclass(frozen=True)
class Candidate:
    key: str
    label: str
    parameters: tuple[str, ...]
    bounds: dict[str, tuple[float, float]]
    initial: dict[str, float]
    sn_active: tuple[str, ...]
    symmetry_class: str
    interpretation: str
    action_complete: bool = False


CANDIDATES = (
    Candidate(
        "plamb_baseline",
        "PLAMB static closure",
        ("p",),
        {"p": (-0.5, 2.5)},
        {"p": 0.8004694772058842},
        ("p",),
        "homogeneous_isotropic",
        "Registered gamma=2.3 static clock-path closure.",
    ),
    Candidate(
        "common_conformal_power",
        "Common conformal power",
        ("p", "s_conformal"),
        {"p": (-0.5, 2.5), "s_conformal": (-3.0, 3.0)},
        {"p": 0.8004694772058842, "s_conformal": 0.0},
        ("p", "s_conformal"),
        "homogeneous_isotropic",
        "Observable proxy for an omitted A(T,X) or common unit-scaling factor.",
    ),
    Candidate(
        "scalar_ruler_power",
        "Scalar ruler evolution",
        ("p", "b_ruler"),
        {"p": (-0.5, 2.5), "b_ruler": (-5.0, 5.0)},
        {"p": 0.8004694772058842, "b_ruler": 0.0},
        ("p",),
        "homogeneous_isotropic",
        "Comoving ruler B(z)=(1+z)^b; invisible to supernovae.",
    ),
    Candidate(
        "spatial_curvature",
        "Spatial-curvature completion",
        ("p", "omega_k"),
        {"p": (-0.5, 2.5), "omega_k": (-2.0, 2.0)},
        {"p": 0.8004694772058842, "omega_k": 0.0},
        ("p", "omega_k"),
        "homogeneous_isotropic",
        "Constant-curvature transverse map S_K(chi).",
    ),
    Candidate(
        "redshift_remap",
        "Atomic redshift remapping",
        ("p", "nu_redshift"),
        {"p": (-0.5, 2.5), "nu_redshift": (0.2, 2.5)},
        {"p": 0.8004694772058842, "nu_redshift": 1.0},
        ("p", "nu_redshift"),
        "homogeneous_isotropic",
        "Clock/atomic proxy y=[(1+z)^nu-1]/nu with dy/dz explicit.",
    ),
    Candidate(
        "clock_quadratic",
        "Clock-law curvature",
        ("p", "gamma_2"),
        {"p": (-0.5, 2.5), "gamma_2": (-3.0, 3.0)},
        {"p": 0.8004694772058842, "gamma_2": 0.0},
        ("p", "gamma_2"),
        "homogeneous_isotropic",
        "One smooth nonlinear term in C(z), retaining C(0)=1.",
    ),
    Candidate(
        "transverse_focusing",
        "Transverse optical focusing",
        ("p", "f_focus"),
        {"p": (-0.5, 2.5), "f_focus": (-4.0, 4.0)},
        {"p": 0.8004694772058842, "f_focus": 0.0},
        ("p", "f_focus"),
        "homogeneous_isotropic_effective",
        "Jacobi-map proxy that changes transverse and luminosity distance, not radial path.",
    ),
    Candidate(
        "anisotropic_ruler",
        "Volume-preserving anisotropic ruler",
        ("p", "e_anisotropy"),
        {"p": (-0.5, 2.5), "e_anisotropy": (-5.0, 5.0)},
        {"p": 0.8004694772058842, "e_anisotropy": 0.0},
        ("p",),
        "breaks_spatial_isotropy",
        "B_perp=(1+z)^(-e/3), B_parallel=(1+z)^(2e/3).",
    ),
    Candidate(
        "radial_closure_slip",
        "Independent radial closure slip",
        ("p", "h_radial"),
        {"p": (-0.5, 2.5), "h_radial": (-5.0, 5.0)},
        {"p": 0.8004694772058842, "h_radial": 0.0},
        ("p",),
        "breaks_distance_closure",
        "Control allowing D_H to move independently of D_M.",
    ),
    Candidate(
        "lcdm_control",
        "Flat Lambda-CDM control",
        ("omega_m",),
        {"omega_m": (0.05, 0.60)},
        {"omega_m": 0.30},
        ("omega_m",),
        "homogeneous_isotropic_control",
        "Standard flat Lambda-CDM distance relation.",
        action_complete=True,
    ),
)
CANDIDATE_BY_KEY = {candidate.key: candidate for candidate in CANDIDATES}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def redshift_mapping(z: np.ndarray, nu: float) -> tuple[np.ndarray, np.ndarray]:
    z = np.asarray(z, dtype=float)
    mapped = np.expm1(nu * np.log1p(z)) / nu
    derivative = np.power(1.0 + z, nu - 1.0)
    return mapped, derivative


def transverse_curvature(chi: np.ndarray, omega_k: float) -> np.ndarray:
    chi = np.asarray(chi, dtype=float)
    if abs(omega_k) < 1e-9:
        return chi * (1.0 + omega_k * np.square(chi) / 6.0)
    root = math.sqrt(abs(omega_k))
    argument = root * chi
    if omega_k > 0.0:
        return np.sinh(argument) / root
    return np.sin(argument) / root


def quadratic_clock_integral(z: np.ndarray, p: float, gamma_2: float) -> np.ndarray:
    base = np.asarray(clock_path_integral(z, GAMMA_FIXED, p), dtype=float)
    extra = closure.power_integral(z, 2.0 - p)
    extra -= 2.0 * closure.power_integral(z, 1.0 - p)
    extra += closure.power_integral(z, -p)
    return base + gamma_2 * extra


def distances(
    candidate: Candidate,
    z: np.ndarray,
    parameters: dict[str, float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    z = np.asarray(z, dtype=float)
    if candidate.key == "lcdm_control":
        omega_m = float(parameters["omega_m"])
        dm = rawmu.lcdm_integral(z, omega_m)
        dh = 1.0 / np.sqrt(omega_m * np.power(1.0 + z, 3.0) + 1.0 - omega_m)
        return (1.0 + z) * dm, dm, dh

    p = float(parameters["p"])
    chi = np.asarray(clock_path_integral(z, GAMMA_FIXED, p), dtype=float)
    dh = np.asarray(clock_path_integrand(z, GAMMA_FIXED, p), dtype=float)
    sn_distance = chi.copy()
    dm = chi.copy()

    if candidate.key == "common_conformal_power":
        factor = np.power(1.0 + z, float(parameters["s_conformal"]))
        sn_distance *= factor
        dm *= factor
        dh *= factor
    elif candidate.key == "scalar_ruler_power":
        ruler = np.power(1.0 + z, float(parameters["b_ruler"]))
        dm /= ruler
        dh /= ruler
    elif candidate.key == "spatial_curvature":
        dm = transverse_curvature(chi, float(parameters["omega_k"]))
        sn_distance = dm.copy()
    elif candidate.key == "redshift_remap":
        mapped, derivative = redshift_mapping(z, float(parameters["nu_redshift"]))
        chi = np.asarray(clock_path_integral(mapped, GAMMA_FIXED, p), dtype=float)
        dh = np.asarray(clock_path_integrand(mapped, GAMMA_FIXED, p), dtype=float)
        dh *= derivative
        sn_distance = chi.copy()
        dm = chi.copy()
    elif candidate.key == "clock_quadratic":
        gamma_2 = float(parameters["gamma_2"])
        c_ratio = 1.0 + GAMMA_FIXED * z + gamma_2 * np.square(z) / (1.0 + z)
        if np.any(c_ratio <= 0.0):
            invalid = np.full_like(z, np.nan)
            return invalid, invalid, invalid
        chi = quadratic_clock_integral(z, p, gamma_2)
        dh = c_ratio / np.power(1.0 + z, p)
        sn_distance = chi.copy()
        dm = chi.copy()
    elif candidate.key == "transverse_focusing":
        factor = np.exp(float(parameters["f_focus"]) * z / (1.0 + z))
        sn_distance *= factor
        dm *= factor
    elif candidate.key == "anisotropic_ruler":
        exponent = float(parameters["e_anisotropy"])
        b_perpendicular = np.power(1.0 + z, -exponent / 3.0)
        b_parallel = np.power(1.0 + z, 2.0 * exponent / 3.0)
        dm /= b_perpendicular
        dh /= b_parallel
    elif candidate.key == "radial_closure_slip":
        dh *= np.power(1.0 + z, float(parameters["h_radial"]))
    elif candidate.key != "plamb_baseline":
        raise ValueError(candidate.key)

    return sn_distance, dm, dh


def sn_score(
    blocks: list[rawmu.Block],
    candidate: Candidate,
    parameters: dict[str, float],
) -> tuple[float, dict[str, float]]:
    total = 0.0
    offsets: dict[str, float] = {}
    for block in blocks:
        distance, _dm, _dh = distances(candidate, block.z, parameters)
        if np.any(~np.isfinite(distance)) or np.any(distance <= 0.0):
            return 1e100, {}
        model = 5.0 * np.log10((C_KMS / H0_FIXED) * distance) + 25.0
        residual = block.mu - model
        ones = np.ones(block.n)
        denominator = float(ones @ block.precision @ ones)
        offset = float(ones @ block.precision @ residual / denominator)
        profiled = residual - offset
        total += float(profiled @ block.precision @ profiled)
        offsets[block.label] = offset
    return total, offsets


def bao_shape(
    data: pd.DataFrame,
    candidate: Candidate,
    parameters: dict[str, float],
) -> np.ndarray:
    _sn, dm, dh = distances(candidate, data["z"].to_numpy(dtype=float), parameters)
    if np.any(~np.isfinite(dm)) or np.any(~np.isfinite(dh)):
        return np.full(len(data), np.nan)
    return closure.observable_vector(data, dm, dh, 0.0)


def bao_score(
    data: pd.DataFrame,
    covariance: np.ndarray,
    candidate: Candidate,
    parameters: dict[str, float],
) -> tuple[float, float, np.ndarray]:
    shape = bao_shape(data, candidate, parameters)
    if np.any(~np.isfinite(shape)) or np.any(shape <= 0.0):
        return 1e100, float("nan"), shape
    precision = np.linalg.inv(covariance)
    observed = data["value"].to_numpy(dtype=float)
    q, chi2 = closure.profile_scale(shape, observed, precision)
    if not np.isfinite(q) or q <= 0.0:
        return 1e100, q, shape
    return chi2, q, shape


def active_parameters(candidate: Candidate, mode: str) -> tuple[str, ...]:
    return candidate.sn_active if mode == "sn" else candidate.parameters


def parameters_from_vector(
    candidate: Candidate,
    names: tuple[str, ...],
    vector: np.ndarray,
) -> dict[str, float]:
    result = dict(candidate.initial)
    result.update({name: float(value) for name, value in zip(names, vector)})
    return result


def numerical_hessian(
    objective,
    optimum: np.ndarray,
    bounds: list[tuple[float, float]],
) -> np.ndarray | None:
    n = len(optimum)
    steps = np.asarray(
        [max(2e-4 * (upper - lower), 1e-5) for lower, upper in bounds]
    )
    for value, step, (lower, upper) in zip(optimum, steps, bounds):
        if value - step <= lower or value + step >= upper:
            return None
    centre = float(objective(optimum))
    hessian = np.zeros((n, n), dtype=float)
    for i in range(n):
        plus = optimum.copy()
        minus = optimum.copy()
        plus[i] += steps[i]
        minus[i] -= steps[i]
        hessian[i, i] = (objective(plus) - 2.0 * centre + objective(minus)) / steps[i] ** 2
        for j in range(i):
            pp = optimum.copy()
            pm = optimum.copy()
            mp = optimum.copy()
            mm = optimum.copy()
            pp[i] += steps[i]
            pp[j] += steps[j]
            pm[i] += steps[i]
            pm[j] -= steps[j]
            mp[i] -= steps[i]
            mp[j] += steps[j]
            mm[i] -= steps[i]
            mm[j] -= steps[j]
            value = (objective(pp) - objective(pm) - objective(mp) + objective(mm))
            value /= 4.0 * steps[i] * steps[j]
            hessian[i, j] = value
            hessian[j, i] = value
    return hessian


def fit_candidate(
    blocks: list[rawmu.Block],
    bao_data: pd.DataFrame,
    bao_covariance: np.ndarray,
    candidate: Candidate,
    mode: str,
) -> tuple[dict[str, object], dict[str, object]]:
    names = active_parameters(candidate, mode)
    bounds = [candidate.bounds[name] for name in names]

    def objective(vector: np.ndarray) -> float:
        parameters = parameters_from_vector(candidate, names, vector)
        sn_chi2, _offsets = sn_score(blocks, candidate, parameters)
        bao_chi2, _q, _shape = bao_score(bao_data, bao_covariance, candidate, parameters)
        if mode == "sn":
            return sn_chi2
        if mode == "bao":
            return bao_chi2
        if mode == "joint":
            return sn_chi2 + bao_chi2
        raise ValueError(mode)

    initial = np.asarray([candidate.initial[name] for name in names], dtype=float)
    starts = [initial]
    if len(names) > 1:
        lower, upper = bounds[-1]
        for fraction in (0.2, 0.8):
            trial = initial.copy()
            trial[-1] = lower + fraction * (upper - lower)
            starts.append(trial)
    candidates: list[tuple[float, np.ndarray, bool, str]] = []
    for start in starts:
        result = minimize(
            objective,
            start,
            method="L-BFGS-B",
            bounds=bounds,
            options={"ftol": 1e-12, "gtol": 1e-8, "maxiter": 500, "maxls": 50},
        )
        candidates.append((float(result.fun), np.asarray(result.x), bool(result.success), str(result.message)))
    minimum, optimum, success, message = min(candidates, key=lambda item: item[0])
    parameters = parameters_from_vector(candidate, names, optimum)
    sn_chi2, offsets = sn_score(blocks, candidate, parameters)
    bao_chi2, q, _shape = bao_score(bao_data, bao_covariance, candidate, parameters)
    n_sn = int(sum(block.n for block in blocks))
    n_bao = len(bao_data)
    nuisance = {"sn": len(blocks), "bao": 1, "joint": len(blocks) + 1}[mode]
    n_data = {"sn": n_sn, "bao": n_bao, "joint": n_sn + n_bao}[mode]
    n_parameters = len(names) + nuisance
    selected_chi2 = {"sn": sn_chi2, "bao": bao_chi2, "joint": sn_chi2 + bao_chi2}[mode]
    boundary = False
    for value, (lower, upper) in zip(optimum, bounds):
        tolerance = 1e-4 * max(upper - lower, 1.0)
        boundary = boundary or value - lower <= tolerance or upper - value <= tolerance

    hessian = numerical_hessian(objective, optimum, bounds)
    identifiability: dict[str, object] = {
        "candidate": candidate.key,
        "mode": mode,
        "active_parameters": json.dumps(names),
        "hessian_available": hessian is not None,
        "positive_definite": False,
        "condition_number": float("nan"),
        "maximum_absolute_correlation": float("nan"),
    }
    if hessian is not None and np.all(np.isfinite(hessian)):
        eigenvalues = np.linalg.eigvalsh(hessian)
        positive = bool(np.min(eigenvalues) > 1e-8)
        identifiability["positive_definite"] = positive
        if positive:
            condition = float(np.max(eigenvalues) / np.min(eigenvalues))
            covariance = 2.0 * np.linalg.inv(hessian)
            sigma = np.sqrt(np.maximum(np.diag(covariance), 0.0))
            correlation = covariance / np.outer(sigma, sigma)
            off_diagonal = correlation - np.eye(len(names))
            identifiability["condition_number"] = condition
            identifiability["maximum_absolute_correlation"] = float(
                np.max(np.abs(off_diagonal))
            )
            for name, value in zip(names, sigma):
                identifiability[f"sigma_{name}"] = float(value)

    row: dict[str, object] = {
        "candidate": candidate.key,
        "label": candidate.label,
        "mode": mode,
        "symmetry_class": candidate.symmetry_class,
        "active_parameters": json.dumps(names),
        "parameters": json.dumps(parameters, sort_keys=True),
        "p": parameters.get("p", np.nan),
        "theta_name": next((name for name in candidate.parameters if name != "p"), ""),
        "theta_value": next((parameters[name] for name in candidate.parameters if name != "p"), np.nan),
        "sn_chi2_at_fit": sn_chi2,
        "bao_chi2_at_fit": bao_chi2,
        "selected_chi2": selected_chi2,
        "N": n_data,
        "N_parameters": n_parameters,
        "BIC": selected_chi2 + n_parameters * math.log(n_data),
        "goodness_upper_tail": float(
            chi2_distribution.sf(selected_chi2, max(n_data - n_parameters, 1))
        ),
        "bao_q_c0_over_H0rd": q,
        "offsets": json.dumps(offsets, sort_keys=True),
        "optimisation_success": success and np.isfinite(minimum),
        "optimisation_message": message,
        "parameter_at_boundary": boundary,
        "all_candidate_parameters_active": set(names) == set(candidate.parameters),
    }
    return row, identifiability


def parameter_dict(row: pd.Series) -> dict[str, float]:
    return {key: float(value) for key, value in json.loads(str(row["parameters"])).items()}


def make_transfer_table(
    fits: pd.DataFrame,
    identifiability: pd.DataFrame,
    blocks: list[rawmu.Block],
    bao_data: pd.DataFrame,
    bao_covariance: np.ndarray,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for candidate in CANDIDATES:
        model_fits = fits[fits["candidate"] == candidate.key].set_index("mode")
        sn_fit = model_fits.loc["sn"]
        bao_fit = model_fits.loc["bao"]
        sn_parameters = parameter_dict(sn_fit)
        bao_parameters = parameter_dict(bao_fit)
        inactive = tuple(
            name for name in candidate.parameters if name not in candidate.sn_active
        )
        conditional_parameters = dict(sn_parameters)
        conditional_success = True
        if inactive:
            bounds = [candidate.bounds[name] for name in inactive]

            def conditional_objective(vector: np.ndarray) -> float:
                trial = dict(sn_parameters)
                trial.update(
                    {name: float(value) for name, value in zip(inactive, vector)}
                )
                return bao_score(
                    bao_data, bao_covariance, candidate, trial
                )[0]

            result = minimize(
                conditional_objective,
                np.asarray([candidate.initial[name] for name in inactive]),
                method="L-BFGS-B",
                bounds=bounds,
                options={"ftol": 1e-12, "gtol": 1e-9, "maxiter": 500},
            )
            conditional_parameters.update(
                {name: float(value) for name, value in zip(inactive, result.x)}
            )
            conditional_success = bool(result.success)
        bao_transfer_chi2, _q, _shape = bao_score(
            bao_data, bao_covariance, candidate, conditional_parameters
        )
        dof = max(len(bao_data) - 1 - len(inactive), 1)
        bao_transfer_p = float(chi2_distribution.sf(bao_transfer_chi2, dof))
        sn_transfer_chi2, _offsets = sn_score(blocks, candidate, bao_parameters)
        sn_best = float(sn_fit["sn_chi2_at_fit"])

        theta_name = str(sn_fit["theta_name"])
        tension = float("nan")
        if theta_name and theta_name in candidate.sn_active:
            sn_id = identifiability[
                (identifiability["candidate"] == candidate.key)
                & (identifiability["mode"] == "sn")
            ].iloc[0]
            bao_id = identifiability[
                (identifiability["candidate"] == candidate.key)
                & (identifiability["mode"] == "bao")
            ].iloc[0]
            sigma_column = f"sigma_{theta_name}"
            if sigma_column in identifiability.columns:
                sigma_sn = float(sn_id.get(sigma_column, np.nan))
                sigma_bao = float(bao_id.get(sigma_column, np.nan))
                if np.isfinite(sigma_sn) and np.isfinite(sigma_bao):
                    difference = float(sn_fit["theta_value"] - bao_fit["theta_value"])
                    tension = abs(difference) / math.sqrt(sigma_sn**2 + sigma_bao**2)
        rows.append(
            {
                "candidate": candidate.key,
                "theta_name": theta_name,
                "SN_invisible_BAO_parameters": json.dumps(inactive),
                "conditional_parameters": json.dumps(
                    conditional_parameters, sort_keys=True
                ),
                "conditional_optimisation_success": conditional_success,
                "bao_chi2_after_SN_training": bao_transfer_chi2,
                "bao_p_after_SN_training": bao_transfer_p,
                "sn_chi2_at_bao_fit": sn_transfer_chi2,
                "delta_sn_chi2_at_bao_fit": sn_transfer_chi2 - sn_best,
                "sn_bao_theta_tension_sigma": tension,
            }
        )
    return pd.DataFrame(rows)


def exploratory_rank2_localisation(
    fits: pd.DataFrame,
    data: pd.DataFrame,
    covariance: np.ndarray,
) -> tuple[pd.DataFrame, dict[str, float]]:
    """Post-registered two-mode control after every rank-one candidate failed."""

    sn_fit = fits[
        (fits["candidate"] == "plamb_baseline") & (fits["mode"] == "sn")
    ].iloc[0]
    p = float(sn_fit["p"])
    z = data["z"].to_numpy(dtype=float)
    observed = data["value"].to_numpy(dtype=float)
    precision = np.linalg.inv(covariance)
    chi = np.asarray(clock_path_integral(z, GAMMA_FIXED, p), dtype=float)
    dh = np.asarray(clock_path_integrand(z, GAMMA_FIXED, p), dtype=float)

    def score(vector: np.ndarray) -> tuple[float, float, np.ndarray]:
        b_ruler, h_radial = np.asarray(vector, dtype=float)
        ruler = np.power(1.0 + z, b_ruler)
        dm_modified = chi / ruler
        dh_modified = dh / ruler * np.power(1.0 + z, h_radial)
        shape = closure.observable_vector(data, dm_modified, dh_modified, 0.0)
        q, chi2 = closure.profile_scale(shape, observed, precision)
        return chi2, q, shape

    result = minimize(
        lambda vector: score(vector)[0],
        np.asarray([1.0, -0.75]),
        method="L-BFGS-B",
        bounds=[(-5.0, 5.0), (-5.0, 5.0)],
        options={"ftol": 1e-13, "gtol": 1e-10, "maxiter": 800},
    )
    b_ruler, h_radial = np.asarray(result.x, dtype=float)
    chi2, q, _shape = score(result.x)
    n_bao_parameters = 3
    bao_dof = len(data) - n_bao_parameters
    sn_chi2 = float(sn_fit["sn_chi2_at_fit"])
    n_joint = int(sn_fit["N"]) + len(data)
    n_joint_parameters = 3 + 1 + 3
    joint_chi2 = sn_chi2 + chi2
    lcdm_bao_bic = float(
        fits[
            (fits["candidate"] == "lcdm_control") & (fits["mode"] == "bao")
        ].iloc[0]["BIC"]
    )
    lcdm_joint_bic = float(
        fits[
            (fits["candidate"] == "lcdm_control") & (fits["mode"] == "joint")
        ].iloc[0]["BIC"]
    )
    bao_bic = chi2 + n_bao_parameters * math.log(len(data))
    joint_bic = joint_chi2 + n_joint_parameters * math.log(n_joint)
    endpoint_one_plus_z = 3.33
    row = {
        "analysis_status": "post_hoc_rank2_localisation_after_registered_rank1_failure",
        "p_fixed_from_SN": p,
        "b_common_ruler": float(b_ruler),
        "h_independent_radial": float(h_radial),
        "B_parallel_power_b_minus_h": float(b_ruler - h_radial),
        "B_perpendicular_z2p33": float(endpoint_one_plus_z**b_ruler),
        "B_parallel_z2p33": float(endpoint_one_plus_z ** (b_ruler - h_radial)),
        "AP_ratio_multiplier_z2p33": float(endpoint_one_plus_z**h_radial),
        "bao_q_c0_over_H0rd": q,
        "bao_chi2": chi2,
        "bao_dof": bao_dof,
        "bao_goodness_upper_tail": float(chi2_distribution.sf(chi2, bao_dof)),
        "bao_BIC": bao_bic,
        "delta_bao_BIC_vs_LCDM": bao_bic - lcdm_bao_bic,
        "sn_chi2": sn_chi2,
        "joint_chi2": joint_chi2,
        "joint_N_parameters": n_joint_parameters,
        "joint_BIC": joint_bic,
        "delta_joint_BIC_vs_LCDM": joint_bic - lcdm_joint_bic,
        "optimisation_success": bool(result.success),
        "parameter_at_boundary": bool(
            abs(b_ruler) >= 4.999 or abs(h_radial) >= 4.999
        ),
        "preserves_background_isotropy": False,
        "preserves_covariant_distance_closure": False,
        "physical_degree_identified": False,
    }
    return pd.DataFrame([row]), {
        "p": p,
        "b_ruler": float(b_ruler),
        "h_radial": float(h_radial),
        "q": q,
    }


def conditional_bao_stress(
    fits: pd.DataFrame,
    data: pd.DataFrame,
    covariance: np.ndarray,
) -> pd.DataFrame:
    observed = data["value"].to_numpy(dtype=float)
    z_all = data["z"].to_numpy(dtype=float)
    rows: list[dict[str, object]] = []
    for candidate in CANDIDATES:
        fit = fits[(fits["candidate"] == candidate.key) & (fits["mode"] == "joint")].iloc[0]
        parameters = parameter_dict(fit)
        shape = bao_shape(data, candidate, parameters)
        q = float(fit["bao_q_c0_over_H0rd"])
        residual = observed - q * shape
        for z in np.unique(z_all):
            test = np.flatnonzero(z_all == z)
            train = np.flatnonzero(z_all != z)
            c_tt = covariance[np.ix_(test, test)]
            if len(train):
                c_tr = covariance[np.ix_(test, train)]
                c_rr = covariance[np.ix_(train, train)]
                projection = c_tr @ np.linalg.solve(c_rr, residual[train])
                conditional_residual = residual[test] - projection
                conditional_covariance = c_tt - c_tr @ np.linalg.solve(
                    c_rr, covariance[np.ix_(train, test)]
                )
            else:
                conditional_residual = residual[test]
                conditional_covariance = c_tt
            precision = np.linalg.inv(conditional_covariance)
            score = float(conditional_residual @ precision @ conditional_residual)
            rows.append(
                {
                    "candidate": candidate.key,
                    "z": z,
                    "N": len(test),
                    "conditional_chi2": score,
                    "conditional_p": float(chi2_distribution.sf(score, len(test))),
                }
            )
    return pd.DataFrame(rows)


def validation_table(
    fits: pd.DataFrame,
    blocks: list[rawmu.Block],
    bao_data: pd.DataFrame,
    bao_covariance: np.ndarray,
) -> pd.DataFrame:
    def fit_row(candidate: str, mode: str) -> pd.Series:
        return fits[(fits["candidate"] == candidate) & (fits["mode"] == mode)].iloc[0]

    baseline_sn = fit_row("plamb_baseline", "sn")
    lcdm_sn = fit_row("lcdm_control", "sn")
    lcdm_bao = fit_row("lcdm_control", "bao")
    fixed_parameters = {"p": 0.8004694772058842}
    fixed_chi2, _q, _shape = bao_score(
        bao_data,
        bao_covariance,
        CANDIDATE_BY_KEY["plamb_baseline"],
        fixed_parameters,
    )
    checks = [
        ("SN sample size", sum(block.n for block in blocks), 3422, 0.0),
        ("PLAMB SN p regression", float(baseline_sn["p"]), EXPECTED["plamb_sn_p"], 5e-5),
        ("PLAMB SN chi2 regression", float(baseline_sn["sn_chi2_at_fit"]), EXPECTED["plamb_sn_chi2"], 0.02),
        ("LCDM SN Omega_m regression", float(json.loads(lcdm_sn["parameters"])["omega_m"]), EXPECTED["lcdm_sn_omega_m"], 5e-5),
        ("LCDM SN chi2 regression", float(lcdm_sn["sn_chi2_at_fit"]), EXPECTED["lcdm_sn_chi2"], 0.02),
        ("LCDM BAO Omega_m regression", float(json.loads(lcdm_bao["parameters"])["omega_m"]), EXPECTED["lcdm_bao_omega_m"], 5e-5),
        ("LCDM BAO chi2 regression", float(lcdm_bao["bao_chi2_at_fit"]), EXPECTED["lcdm_bao_chi2"], 0.02),
        ("Fixed PLAMB BAO chi2 regression", fixed_chi2, EXPECTED["fixed_plamb_bao_chi2"], 0.05),
        ("BAO covariance minimum eigenvalue positive", float(np.min(np.linalg.eigvalsh(bao_covariance))), 0.0, 0.0),
    ]
    rows = []
    for name, value, target, tolerance in checks:
        if name.endswith("positive"):
            passed = value > target
        else:
            passed = abs(value - target) <= tolerance
        rows.append(
            {
                "check": name,
                "value": value,
                "target": target,
                "absolute_tolerance": tolerance,
                "passed": bool(passed),
            }
        )
    return pd.DataFrame(rows)


def decision_table(
    fits: pd.DataFrame,
    transfers: pd.DataFrame,
    identifiability: pd.DataFrame,
) -> pd.DataFrame:
    joint = fits[fits["mode"] == "joint"].set_index("candidate")
    baseline_bic = float(joint.loc["plamb_baseline", "BIC"])
    lcdm_bic = float(joint.loc["lcdm_control", "BIC"])
    rows = []
    for candidate in CANDIDATES:
        fit = joint.loc[candidate.key]
        transfer = transfers[transfers["candidate"] == candidate.key].iloc[0]
        identity = identifiability[
            (identifiability["candidate"] == candidate.key)
            & (identifiability["mode"] == "joint")
        ].iloc[0]
        identifiable = bool(identity["positive_definite"]) and float(
            identity["condition_number"]
        ) < 1e8
        bao_dof = max(13 - (len(candidate.parameters) + 1), 1)
        bao_p_joint = float(chi2_distribution.sf(float(fit["bao_chi2_at_fit"]), bao_dof))
        transfer_pass = bool(
            bool(transfer["conditional_optimisation_success"])
            and float(transfer["bao_p_after_SN_training"]) >= 0.05
        )
        symmetry_pass = not candidate.symmetry_class.startswith("breaks_")
        gates = {
            "delta_BIC_joint_at_most_minus_10": float(fit["BIC"] - baseline_bic) <= -10.0,
            "joint_BAO_goodness_p_at_least_0p05": bao_p_joint >= 0.05,
            "SN_trained_BAO_transfer": transfer_pass,
            "joint_hessian_identifiable": identifiable,
            "parameter_not_at_boundary": not bool(fit["parameter_at_boundary"]),
            "preserves_background_symmetry": symmetry_pass,
            "within_10_BIC_of_LCDM": float(fit["BIC"] - lcdm_bic) <= 10.0,
            "action_level_derivation_supplied": candidate.action_complete,
        }
        observational = all(
            gates[name]
            for name in (
                "delta_BIC_joint_at_most_minus_10",
                "joint_BAO_goodness_p_at_least_0p05",
                "SN_trained_BAO_transfer",
                "joint_hessian_identifiable",
                "parameter_not_at_boundary",
                "preserves_background_symmetry",
                "within_10_BIC_of_LCDM",
            )
        )
        rows.append(
            {
                "candidate": candidate.key,
                "symmetry_class": candidate.symmetry_class,
                "delta_BIC_joint_vs_PLAMB": float(fit["BIC"] - baseline_bic),
                "delta_BIC_joint_vs_LCDM": float(fit["BIC"] - lcdm_bic),
                "bao_chi2_at_joint": float(fit["bao_chi2_at_fit"]),
                "bao_p_at_joint": bao_p_joint,
                **gates,
                "observational_signature_promoted": observational,
                "physical_degree_identified": observational and candidate.action_complete,
            }
        )
    return pd.DataFrame(rows)


def observed_ap_ratios(data: pd.DataFrame, covariance: np.ndarray) -> pd.DataFrame:
    return closure.observed_ap_ratios(data, covariance)


def make_plot(
    path: Path,
    fits: pd.DataFrame,
    decisions: pd.DataFrame,
    data: pd.DataFrame,
    covariance: np.ndarray,
    rank2_parameters: dict[str, float],
) -> None:
    joint = fits[fits["mode"] == "joint"].copy()
    joint = joint.merge(decisions[["candidate", "delta_BIC_joint_vs_PLAMB"]], on="candidate")
    joint = joint.sort_values("delta_BIC_joint_vs_PLAMB")
    colours = [
        "#b33a3a" if str(value).startswith("breaks_") else "#287a5a"
        for value in joint["symmetry_class"]
    ]
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.2), constrained_layout=True)
    axes[0].barh(joint["label"], joint["delta_BIC_joint_vs_PLAMB"], color=colours)
    axes[0].axvline(0.0, color="black", linewidth=1.0)
    axes[0].axvline(-10.0, color="black", linewidth=1.0, linestyle="--")
    axes[0].set_xlabel("Delta BIC relative to PLAMB baseline")
    axes[0].set_title("Joint SN plus DESI model ranking")
    axes[0].grid(axis="x", alpha=0.25)

    ap = observed_ap_ratios(data, covariance)
    axes[1].errorbar(
        ap["z"], ap["ratio"], yerr=ap["sigma"], fmt="o", color="black", capsize=3, label="DESI DR2"
    )
    selected = list(joint["candidate"].head(3))
    if "lcdm_control" not in selected:
        selected.append("lcdm_control")
    z = np.linspace(0.25, 2.33, 500)
    palette = ("#287a5a", "#b33a3a", "#5d55a5", "#2672b8")
    for colour, key in zip(palette, selected):
        candidate = CANDIDATE_BY_KEY[key]
        fit = joint[joint["candidate"] == key].iloc[0]
        parameters = parameter_dict(fit)
        _sn, dm, dh = distances(candidate, z, parameters)
        axes[1].plot(z, dh / dm, color=colour, label=candidate.label)
    p_rank2 = rank2_parameters["p"]
    chi_rank2 = np.asarray(
        clock_path_integral(z, GAMMA_FIXED, p_rank2), dtype=float
    )
    dh_rank2 = np.asarray(
        clock_path_integrand(z, GAMMA_FIXED, p_rank2), dtype=float
    )
    rank2_ratio = dh_rank2 / chi_rank2
    rank2_ratio *= np.power(1.0 + z, rank2_parameters["h_radial"])
    axes[1].plot(
        z,
        rank2_ratio,
        color="#d07a21",
        linestyle="--",
        label="Post hoc rank-two closure",
    )
    axes[1].set_xlim(0.25, 2.4)
    axes[1].set_ylim(0.0, 5.0)
    axes[1].set_xlabel("Redshift z")
    axes[1].set_ylabel("D_H / D_M")
    axes[1].set_title("Ruler-independent ratio at joint fit")
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_report(
    path: Path,
    fits: pd.DataFrame,
    transfers: pd.DataFrame,
    identifiability: pd.DataFrame,
    stresses: pd.DataFrame,
    validations: pd.DataFrame,
    decisions: pd.DataFrame,
    rank2: pd.DataFrame,
) -> None:
    joint = fits[fits["mode"] == "joint"].copy()
    joint = joint.merge(
        decisions[
            [
                "candidate",
                "delta_BIC_joint_vs_PLAMB",
                "delta_BIC_joint_vs_LCDM",
                "bao_p_at_joint",
                "observational_signature_promoted",
                "physical_degree_identified",
            ]
        ],
        on="candidate",
    ).sort_values("BIC")
    top_noncontrol = joint[
        ~joint["candidate"].isin(["plamb_baseline", "lcdm_control"])
    ].iloc[0]
    rank2_row = rank2.iloc[0]
    worst_stress = (
        stresses.groupby("candidate", as_index=False)
        .apply(lambda frame: frame.loc[frame["conditional_chi2"].idxmax()], include_groups=False)
        .reset_index(drop=True)
    )
    lines = [
        "# PLAMB missing-degree-of-freedom diagnostic suite",
        "",
        f"Generated: {datetime.now():%Y-%m-%d %H:%M:%S}",
        "",
        "## Decision",
        "",
        "This suite searches for an **observable signature** of an omitted degree of freedom. It cannot identify a "
        "fundamental field without an action and independent coupling law.",
        "",
        f"The largest numerical one-degree compromise is `{top_noncontrol['candidate']}` with joint "
        f"`Delta BIC={float(top_noncontrol['delta_BIC_joint_vs_PLAMB']):.6f}` relative to the static PLAMB "
        f"baseline and DESI `chi-squared={float(top_noncontrol['bao_chi2_at_fit']):.6f}` at its joint fit. "
        "It is boundary-seeking and fails the DESI goodness gate.",
        "",
        f"Promoted observational signatures: `{int(decisions['observational_signature_promoted'].sum())}`. "
        f"Physically identified new degrees of freedom: `{int(decisions['physical_degree_identified'].sum())}`.",
        "The registered result therefore disfavours the hypothesis that one omitted scalar background variable "
        "repairs the existing SN--BAO closure. The post hoc rank test requires two independent observable modes: "
        "a common scale evolution and a separate radial/transverse distortion.",
        "",
        "## Registered observable basis",
        "",
        "The common static path is",
        "",
        "$$",
        "\\begin{aligned}",
        "C(z)                         &= 1+2.3z, \\\\",
        "k(z)                         &= \\frac{C(z)}{(1+z)^p}, \\\\",
        "\\chi(z)                    &= \\int_0^z k(u)\\,du, \\\\",
        "D_M                          &= \\chi, \\\\",
        "D_H                          &= k.",
        "\\end{aligned}",
        "$$",
        "",
        "Each extension adds at most one shape variable beyond `p`:",
        "",
        pd.DataFrame(
            [
                {
                    "candidate": candidate.key,
                    "parameters": ", ".join(candidate.parameters),
                    "symmetry": candidate.symmetry_class,
                    "interpretation": candidate.interpretation,
                }
                for candidate in CANDIDATES
            ]
        ).to_markdown(index=False),
        "",
        "The anisotropic-ruler and radial-slip rows are controls. A good fit from either would locate the missing "
        "observable freedom in the radial/transverse closure, but would not establish a covariant or isotropic theory.",
        "",
        "## Joint ranking",
        "",
        joint[
            [
                "candidate",
                "p",
                "theta_name",
                "theta_value",
                "sn_chi2_at_fit",
                "bao_chi2_at_fit",
                "BIC",
                "delta_BIC_joint_vs_PLAMB",
                "delta_BIC_joint_vs_LCDM",
                "bao_p_at_joint",
                "parameter_at_boundary",
            ]
        ].to_markdown(index=False),
        "",
        "## Cross-probe transfer",
        "",
        transfers.to_markdown(index=False),
        "",
        "A candidate that exists only in the BAO likelihood cannot be validated by supernova-to-BAO transfer. "
        "For those rows, the suite fixes the supernova-trained `p` and fits only the SN-invisible BAO variable. "
        "A scalar candidate visible to supernovae must predict DESI without retuning its extra variable.",
        "",
        "## Exploratory rank-two localisation",
        "",
        "After every registered one-degree extension failed, a deliberately post hoc rank test combined a common "
        "ruler power with an independent radial slip while fixing `p` to the supernova-only result:",
        "",
        "$$",
        "\\begin{aligned}",
        "D_M(z)                      &= \\frac{\\chi(z)}{(1+z)^b}, \\\\",
        "D_H(z)                      &= \\frac{k(z)}{(1+z)^b}(1+z)^h, \\\\",
        "B_\\perp(z)                 &= (1+z)^b, \\\\",
        "B_\\parallel(z)             &= (1+z)^{b-h}.",
        "\\end{aligned}",
        "$$",
        "",
        rank2.to_markdown(index=False),
        "",
        f"This reaches DESI `chi-squared={float(rank2_row['bao_chi2']):.6f}` with "
        f"`p={float(rank2_row['p_fixed_from_SN']):.6f}`, `b={float(rank2_row['b_common_ruler']):.6f}` "
        f"and `h={float(rank2_row['h_independent_radial']):.6f}`. It shows that two observational modes can "
        "span most of the discrepancy. At `z=2.33` they require transverse and radial ruler factors of "
        f"`{float(rank2_row['B_perpendicular_z2p33']):.6f}` and "
        f"`{float(rank2_row['B_parallel_z2p33']):.6f}`. The joint result remains "
        f"`Delta BIC={float(rank2_row['delta_joint_BIC_vs_LCDM']):+.6f}` versus Lambda-CDM. Because the "
        "second mode breaks the covariant radial/transverse closure and has no action-level source, this is "
        "localisation evidence, not a new physical model.",
        "",
        "## Identifiability",
        "",
        identifiability[
            [
                "candidate",
                "mode",
                "active_parameters",
                "positive_definite",
                "condition_number",
                "maximum_absolute_correlation",
            ]
        ].to_markdown(index=False),
        "",
        "## Largest conditional DESI stress",
        "",
        worst_stress.to_markdown(index=False),
        "",
        "## Promotion gates",
        "",
        decisions.to_markdown(index=False),
        "",
        "A physical identification additionally requires an action-level derivation of the candidate and its "
        "clock, photon, atomic and ruler couplings. None of the PLAMB proxy extensions supplies that derivation.",
        "",
        "## Regression checks",
        "",
        validations.to_markdown(index=False),
        "",
        "## Claim boundary",
        "",
        "A lower BIC identifies where the present observable closure is deficient; it does not identify SU(2), "
        "antimatter, a new field or broken Lorentz symmetry. Candidates that break isotropy must also predict "
        "direction-dependent BAO and sky multipoles before they can be interpreted physically.",
        "",
        "## Primary references",
        "",
        "- [General disformal scalar-tensor transformations](https://arxiv.org/abs/1412.6210)",
        "- [Homogeneous and isotropic SU(2) cosmology](https://arxiv.org/abs/1012.2861)",
        "- [Etherington distance-duality assumptions](https://arxiv.org/abs/1612.08784)",
        "- [DESI DR2 publications and products](https://data.desi.lbl.gov/doc/papers/dr2/)",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--mean", type=Path, default=closure.DEFAULT_MEAN)
    parser.add_argument("--cov", type=Path, default=closure.DEFAULT_COV)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    protocol = {
        "analysis_id": "plamb_missing_degree_suite_2026-07-19",
        "status": "preregistered_observable_proxy_diagnostic",
        "gamma_c_fixed": GAMMA_FIXED,
        "candidate_keys": [candidate.key for candidate in CANDIDATES],
        "fit_modes": ["sn", "bao", "joint"],
        "primary_SN_covariance": "released_total_primary",
        "primary_redshift_frame": "zHD for Pantheon+ and DES; release z for Union3.1",
        "gates": {
            "delta_BIC_joint": "<= -10 versus PLAMB baseline",
            "DESI_goodness": "p >= 0.05 at joint fit",
            "cross_probe_transfer": (
                "SN shape held fixed; SN-visible variables are not retuned, while any "
                "SN-invisible BAO variable may be fitted once; DESI p >= 0.05"
            ),
            "identifiability": "positive Hessian and condition number < 1e8",
            "boundary": "no active parameter at registered bound",
            "symmetry": "background isotropy retained",
            "control_competitiveness": "within 10 BIC of flat Lambda-CDM",
            "physical_identification": "all statistical gates plus action-level derivation",
        },
        "exploratory_followup": (
            "A rank-two common-ruler plus radial-slip control may be run only after all "
            "registered rank-one candidates fail; it is labelled post hoc and cannot promote."
        ),
        "claim_boundary": "locates observable closure deficiency; does not identify a fundamental field",
    }
    protocol_path = args.outdir / f"plamb_missing_degree_protocol_{DATE_TAG}.json"
    protocol_path.write_text(json.dumps(protocol, indent=2), encoding="utf-8")

    variants, load_metadata = rawmu.load_blocks(MIN_SURVEY_N)
    blocks = variants["released_total_primary"]
    bao_data, bao_covariance = closure.load_bao(args.mean, args.cov)

    fit_rows: list[dict[str, object]] = []
    identity_rows: list[dict[str, object]] = []
    for candidate in CANDIDATES:
        for mode in ("sn", "bao", "joint"):
            row, identity = fit_candidate(
                blocks, bao_data, bao_covariance, candidate, mode
            )
            fit_rows.append(row)
            identity_rows.append(identity)
            print(
                f"{candidate.key} {mode}: chi2={float(row['selected_chi2']):.6f} "
                f"SN={float(row['sn_chi2_at_fit']):.6f} BAO={float(row['bao_chi2_at_fit']):.6f}"
            )
    fits = pd.DataFrame(fit_rows)
    identifiability = pd.DataFrame(identity_rows)
    transfers = make_transfer_table(
        fits, identifiability, blocks, bao_data, bao_covariance
    )
    rank2, rank2_parameters = exploratory_rank2_localisation(
        fits, bao_data, bao_covariance
    )
    stresses = conditional_bao_stress(fits, bao_data, bao_covariance)
    validations = validation_table(fits, blocks, bao_data, bao_covariance)
    decisions = decision_table(fits, transfers, identifiability)

    prefix = "plamb_missing_degree"
    paths = {
        "fits": args.outdir / f"{prefix}_fits_{DATE_TAG}.csv",
        "transfers": args.outdir / f"{prefix}_cross_probe_transfer_{DATE_TAG}.csv",
        "identifiability": args.outdir / f"{prefix}_identifiability_{DATE_TAG}.csv",
        "stress": args.outdir / f"{prefix}_conditional_bao_stress_{DATE_TAG}.csv",
        "validation": args.outdir / f"{prefix}_validation_{DATE_TAG}.csv",
        "decision": args.outdir / f"{prefix}_decisions_{DATE_TAG}.csv",
        "rank2": args.outdir / f"{prefix}_exploratory_rank2_{DATE_TAG}.csv",
        "report": args.outdir / f"{prefix}_report_{DATE_TAG}.md",
        "plot": args.outdir / f"{prefix}_readout_{DATE_TAG}.png",
        "summary": args.outdir / f"{prefix}_summary_{DATE_TAG}.json",
    }
    fits.to_csv(paths["fits"], index=False)
    transfers.to_csv(paths["transfers"], index=False)
    identifiability.to_csv(paths["identifiability"], index=False)
    stresses.to_csv(paths["stress"], index=False)
    validations.to_csv(paths["validation"], index=False)
    decisions.to_csv(paths["decision"], index=False)
    rank2.to_csv(paths["rank2"], index=False)
    write_report(
        paths["report"],
        fits,
        transfers,
        identifiability,
        stresses,
        validations,
        decisions,
        rank2,
    )
    make_plot(
        paths["plot"],
        fits,
        decisions,
        bao_data,
        bao_covariance,
        rank2_parameters,
    )

    joint = fits[fits["mode"] == "joint"].sort_values("BIC")
    summary = {
        "analysis_date": DATE_TAG,
        "N_SN": int(sum(block.n for block in blocks)),
        "N_BAO": len(bao_data),
        "load_metadata": load_metadata,
        "regression_checks_pass": bool(validations["passed"].all()),
        "best_joint_candidate": str(joint.iloc[0]["candidate"]),
        "best_joint_BIC": float(joint.iloc[0]["BIC"]),
        "observational_signatures_promoted": decisions.loc[
            decisions["observational_signature_promoted"], "candidate"
        ].tolist(),
        "physical_degrees_identified": decisions.loc[
            decisions["physical_degree_identified"], "candidate"
        ].tolist(),
        "exploratory_rank2": rank2.iloc[0].to_dict(),
        "further_sampling_authorised": bool(decisions["physical_degree_identified"].any()),
        "joint_ranking": joint[
            ["candidate", "BIC", "sn_chi2_at_fit", "bao_chi2_at_fit", "parameters"]
        ].to_dict(orient="records"),
    }
    paths["summary"].write_text(json.dumps(summary, indent=2), encoding="utf-8")

    source_paths = [
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
        {"kind": "source", "path": str(path), "sha256": sha256_file(path), "bytes": path.stat().st_size}
        for path in source_paths
    ]
    manifest_rows.extend(
        {"kind": "output", "path": str(path), "sha256": sha256_file(path), "bytes": path.stat().st_size}
        for path in [protocol_path, *paths.values()]
    )
    manifest_path = args.outdir / f"{prefix}_manifest_{DATE_TAG}.csv"
    pd.DataFrame(manifest_rows).to_csv(manifest_path, index=False)

    print(f"Regression checks: {'PASS' if validations['passed'].all() else 'FAIL'}")
    print(
        "Physical degrees identified: "
        + (", ".join(summary["physical_degrees_identified"]) or "none")
    )
    print(f"Saved report: {paths['report']}")


if __name__ == "__main__":
    main()
