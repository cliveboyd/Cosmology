#!/usr/bin/env python3
"""Audit PLAMB/Lambda-CDM equivalence and high-redshift Lyman-alpha BAO tests."""

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
from matplotlib.patches import Ellipse
import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
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
import search_plamb_radial_transverse_2026_07_19 as radial  # noqa: E402
import test_plamb_missing_degree_suite_2026_07_19 as suite  # noqa: E402
from plamb_clock_distance import clock_path_integral, clock_path_integrand  # noqa: E402


DATE_TAG = "2026-07-19"
GAMMA_FIXED = 2.3
MIN_SURVEY_N = 10
POWER_BOUNDS = radial.POWER_BOUNDS
DATA_DIR = (
    REPO_ROOT
    / "github_export"
    / "data"
    / "bao"
    / "lya_validation_2026-07-19"
)
MEASUREMENT_CONFIG = DATA_DIR / "lya_published_measurements_2026-07-19.json"
DEFAULT_OUTDIR = (
    REPO_ROOT
    / "github_export"
    / "results"
    / DATE_TAG
    / "plamb_lcdm_lya_validation"
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def covariance_2d(
    sigma_dm: float,
    sigma_dh: float,
    rho: float,
    systematic_fraction: float = 0.0,
    mean_dm: float | None = None,
    mean_dh: float | None = None,
) -> np.ndarray:
    covariance = np.asarray(
        [
            [sigma_dm**2, rho * sigma_dm * sigma_dh],
            [rho * sigma_dm * sigma_dh, sigma_dh**2],
        ],
        dtype=float,
    )
    if systematic_fraction:
        if mean_dm is None or mean_dh is None:
            raise ValueError("Means are required when adding a fractional systematic")
        covariance += np.diag(
            [
                (systematic_fraction * mean_dm) ** 2,
                (systematic_fraction * mean_dh) ** 2,
            ]
        )
    if np.min(np.linalg.eigvalsh(covariance)) <= 0.0:
        raise ValueError("Measurement covariance is not positive definite")
    return covariance


def predict_plamb(
    z: float | np.ndarray,
    p: float,
    b_perpendicular: float,
    b_parallel: float,
    q: float,
) -> tuple[np.ndarray, np.ndarray]:
    z_array = np.atleast_1d(np.asarray(z, dtype=float))
    chi = np.asarray(clock_path_integral(z_array, GAMMA_FIXED, p), dtype=float)
    radial_distance = np.asarray(
        clock_path_integrand(z_array, GAMMA_FIXED, p), dtype=float
    )
    dm = q * chi / np.power(1.0 + z_array, b_perpendicular)
    dh = q * radial_distance / np.power(1.0 + z_array, b_parallel)
    return dm, dh


def predict_lcdm(
    z: float | np.ndarray,
    omega_m: float,
    q: float,
) -> tuple[np.ndarray, np.ndarray]:
    z_array = np.atleast_1d(np.asarray(z, dtype=float))
    candidate = suite.CANDIDATE_BY_KEY["lcdm_control"]
    _sn_distance, dm, dh = suite.distances(
        candidate, z_array, {"omega_m": omega_m}
    )
    return q * np.asarray(dm, dtype=float), q * np.asarray(dh, dtype=float)


def model_prediction(
    model: str,
    z: float,
    plamb_fit: dict[str, object],
    lcdm_fit: dict[str, object],
) -> np.ndarray:
    if model == "PLAMB_frozen_lower_z":
        dm, dh = predict_plamb(
            z,
            float(plamb_fit["p"]),
            float(plamb_fit["b_perpendicular"]),
            float(plamb_fit["b_parallel"]),
            float(plamb_fit["q_c0_over_H0rd"]),
        )
    elif model == "LambdaCDM_frozen_lower_z":
        dm, dh = predict_lcdm(
            z,
            float(lcdm_fit["omega_m"]),
            float(lcdm_fit["q_c0_over_H0rd"]),
        )
    else:
        raise ValueError(model)
    return np.asarray([dm[0], dh[0]], dtype=float)


def fit_equivalence(
    data: pd.DataFrame,
    covariance: np.ndarray,
    p: float,
    lcdm_fit: dict[str, object],
    full_data_plamb_fit: dict[str, object],
    lower_z_plamb_fit: dict[str, object],
) -> tuple[dict[str, object], pd.DataFrame, pd.DataFrame]:
    z = data["z"].to_numpy(dtype=float)
    lcdm_candidate = suite.CANDIDATE_BY_KEY["lcdm_control"]
    _sn, dm_lcdm_base, dh_lcdm_base = suite.distances(
        lcdm_candidate,
        z,
        {"omega_m": float(lcdm_fit["omega_m"])},
    )
    target_shape = closure.observable_vector(
        data,
        np.asarray(dm_lcdm_base, dtype=float),
        np.asarray(dh_lcdm_base, dtype=float),
        0.0,
    )
    target = float(lcdm_fit["q_c0_over_H0rd"]) * target_shape
    precision = np.linalg.inv(covariance)

    def score(vector: np.ndarray) -> tuple[float, float, np.ndarray]:
        shape = radial.plamb_shape(data, p, float(vector[0]), float(vector[1]))
        q, chi2 = closure.profile_scale(shape, target, precision)
        return chi2, q, q * shape

    optimum, _minimum, success, message = radial.optimise_vector(
        lambda vector: score(vector)[0],
        np.asarray(
            [
                float(full_data_plamb_fit["b_perpendicular"]),
                float(full_data_plamb_fit["b_parallel"]),
            ]
        ),
        [POWER_BOUNDS, POWER_BOUNDS],
    )
    chi2_best, q_best, prediction_best = score(optimum)
    full_vector = np.asarray(
        [
            full_data_plamb_fit["b_perpendicular"],
            full_data_plamb_fit["b_parallel"],
        ],
        dtype=float,
    )
    lower_vector = np.asarray(
        [
            lower_z_plamb_fit["b_perpendicular"],
            lower_z_plamb_fit["b_parallel"],
        ],
        dtype=float,
    )
    chi2_full_data_powers, q_full_data_powers, _ = score(full_vector)
    chi2_lower_z_powers, q_lower_z_powers, _ = score(lower_vector)

    hessian, parameter_covariance = radial.hessian_readout(
        data, covariance, p, full_vector
    )
    difference = optimum - full_vector
    if np.all(np.isfinite(parameter_covariance)):
        parameter_distance = float(
            math.sqrt(difference @ np.linalg.solve(parameter_covariance, difference))
        )
    else:
        parameter_distance = float("nan")

    fractional = prediction_best / target - 1.0
    summary: dict[str, object] = {
        "p_SN_trained": p,
        "omega_m_LambdaCDM_lower_z_joint": float(lcdm_fit["omega_m"]),
        "q_LambdaCDM_lower_z_joint": float(lcdm_fit["q_c0_over_H0rd"]),
        "b_perpendicular_equivalent": float(optimum[0]),
        "b_parallel_equivalent": float(optimum[1]),
        "delta_b_equivalent": float(optimum[1] - optimum[0]),
        "q_PLAMB_equivalent": q_best,
        "chi2_PLAMB_to_LambdaCDM_geometry_best": chi2_best,
        "dof_geometry_vector_minus_profiled_parameters": int(len(data) - 3),
        "goodness_p_geometry": float(
            chi2_distribution.sf(chi2_best, max(len(data) - 3, 1))
        ),
        "chi2_at_full_data_fitted_powers": chi2_full_data_powers,
        "q_at_full_data_fitted_powers": q_full_data_powers,
        "chi2_at_lower_z_fitted_powers": chi2_lower_z_powers,
        "q_at_lower_z_fitted_powers": q_lower_z_powers,
        "full_data_b_perpendicular": float(full_vector[0]),
        "full_data_b_parallel": float(full_vector[1]),
        "equivalent_minus_data_parameter_distance_sigma": parameter_distance,
        "maximum_absolute_fractional_geometry_difference": float(
            np.max(np.abs(fractional))
        ),
        "RMS_fractional_geometry_difference": float(
            np.sqrt(np.mean(np.square(fractional)))
        ),
        "optimisation_success": success,
        "optimisation_message": message,
        "data_hessian_positive_definite": bool(hessian.get("positive_definite", False)),
    }

    vector_rows = []
    for index, row in data.iterrows():
        vector_rows.append(
            {
                "z": float(row["z"]),
                "kind": str(row["kind"]),
                "LambdaCDM_lower_z_prediction": target[index],
                "PLAMB_equivalent_prediction": prediction_best[index],
                "fractional_PLAMB_minus_LambdaCDM": fractional[index],
            }
        )

    paired_z = []
    for z_value in np.unique(z):
        kinds = set(data.loc[np.isclose(data["z"], z_value), "kind"].astype(str))
        if any(kind.startswith("DM") for kind in kinds) and any(
            kind.startswith("DH") for kind in kinds
        ):
            paired_z.append(float(z_value))
    paired = np.asarray(paired_z, dtype=float)
    chi = np.asarray(clock_path_integral(paired, GAMMA_FIXED, p), dtype=float)
    k = np.asarray(clock_path_integrand(paired, GAMMA_FIXED, p), dtype=float)
    _sn, dm_lcdm, dh_lcdm = suite.distances(
        lcdm_candidate,
        paired,
        {"omega_m": float(lcdm_fit["omega_m"])},
    )
    q_lcdm = float(lcdm_fit["q_c0_over_H0rd"])
    effective_perpendicular = np.log(q_best * chi / (q_lcdm * dm_lcdm)) / np.log1p(
        paired
    )
    effective_parallel = np.log(q_best * k / (q_lcdm * dh_lcdm)) / np.log1p(
        paired
    )

    dense_z = np.linspace(0.05, 2.50, 2001)
    dense_log_z = np.log1p(dense_z)
    dense_chi = np.asarray(
        clock_path_integral(dense_z, GAMMA_FIXED, p), dtype=float
    )
    dense_k = np.asarray(
        clock_path_integrand(dense_z, GAMMA_FIXED, p), dtype=float
    )
    _sn, dense_dm_lcdm, dense_dh_lcdm = suite.distances(
        lcdm_candidate,
        dense_z,
        {"omega_m": float(lcdm_fit["omega_m"])},
    )
    local_perpendicular = np.gradient(np.log(dense_chi), dense_log_z) - np.gradient(
        np.log(dense_dm_lcdm), dense_log_z
    )
    local_parallel = np.gradient(np.log(dense_k), dense_log_z) - np.gradient(
        np.log(dense_dh_lcdm), dense_log_z
    )
    effective_table = pd.DataFrame(
        {
            "z": paired,
            "pointwise_b_perpendicular": effective_perpendicular,
            "pointwise_b_parallel": effective_parallel,
            "pointwise_delta_b": effective_parallel - effective_perpendicular,
            "local_derivative_b_perpendicular": np.interp(
                paired, dense_z, local_perpendicular
            ),
            "local_derivative_b_parallel": np.interp(
                paired, dense_z, local_parallel
            ),
        }
    )
    return summary, effective_table, pd.DataFrame(vector_rows)


def alpha_measurement(
    measurement_id: str,
    classification: str,
    z: float,
    alpha: dict[str, float],
    fiducial_dm: float,
    fiducial_dh: float,
    systematic_fraction: float = 0.0,
) -> dict[str, object]:
    mean_dm = fiducial_dm * float(alpha["alpha_perpendicular"])
    mean_dh = fiducial_dh * float(alpha["alpha_parallel"])
    sigma_dm = fiducial_dm * float(alpha["sigma_alpha_perpendicular"])
    sigma_dh = fiducial_dh * float(alpha["sigma_alpha_parallel"])
    covariance = covariance_2d(
        sigma_dm,
        sigma_dh,
        float(alpha["rho_parallel_perpendicular"]),
        systematic_fraction,
        mean_dm,
        mean_dh,
    )
    return {
        "measurement_id": measurement_id,
        "classification": classification,
        "z": z,
        "mean": np.asarray([mean_dm, mean_dh], dtype=float),
        "covariance": covariance,
        "systematic_fraction_added": systematic_fraction,
    }


def direct_measurement(
    measurement_id: str,
    classification: str,
    z: float,
    values: dict[str, float],
    systematic_fraction: float = 0.0,
) -> dict[str, object]:
    mean_dm = float(values["DM_over_rd"])
    mean_dh = float(values["DH_over_rd"])
    covariance = covariance_2d(
        float(values["sigma_DM"]),
        float(values["sigma_DH"]),
        float(values["rho_DM_DH"]),
        systematic_fraction,
        mean_dm,
        mean_dh,
    )
    return {
        "measurement_id": measurement_id,
        "classification": classification,
        "z": z,
        "mean": np.asarray([mean_dm, mean_dh], dtype=float),
        "covariance": covariance,
        "systematic_fraction_added": systematic_fraction,
    }


def build_gaussian_measurements(config: dict[str, object]) -> list[dict[str, object]]:
    measurements: list[dict[str, object]] = []
    dr2 = config["desi_dr2"]
    dr2_class = str(dr2["classification"])
    measurements.append(
        direct_measurement(
            "DESI_DR2_combined_stat_sys",
            dr2_class,
            float(dr2["z"]),
            dr2["combined_stat_sys"],
        )
    )
    combined_alpha = dr2["combined_alpha_stat"]
    combined_direct = dr2["combined_stat_sys"]
    fiducial_dm_dr2 = float(combined_direct["DM_over_rd"]) / float(
        combined_alpha["alpha_perpendicular"]
    )
    fiducial_dh_dr2 = float(combined_direct["DH_over_rd"]) / float(
        combined_alpha["alpha_parallel"]
    )
    systematic_fraction = float(dr2["theory_systematic_fraction_per_axis"])
    for label, key in (("auto", "auto_alpha_stat"), ("cross", "cross_alpha_stat")):
        measurements.append(
            alpha_measurement(
                f"DESI_DR2_{label}_plus_0p3pct_theory_budget",
                dr2_class,
                float(dr2["z"]),
                dr2[key],
                fiducial_dm_dr2,
                fiducial_dh_dr2,
                systematic_fraction,
            )
        )

    multipoles = config["desi_dr2_multipoles"]
    measurements.append(
        direct_measurement(
            "DESI_DR2_multipoles_stat_only",
            str(multipoles["classification"]),
            float(multipoles["z"]),
            multipoles,
        )
    )
    measurements.append(
        direct_measurement(
            "DESI_DR2_multipoles_plus_0p3pct_theory_budget",
            str(multipoles["classification"]),
            float(multipoles["z"]),
            multipoles,
            systematic_fraction,
        )
    )

    dr1 = config["desi_dr1"]
    dr1_class = str(dr1["classification"])
    measurements.append(
        direct_measurement(
            "DESI_DR1_combined",
            dr1_class,
            float(dr1["z"]),
            dr1["combined"],
        )
    )
    fiducial_dm_dr1 = float(dr1["combined"]["DM_over_rd"]) / float(
        dr1["combined"]["alpha_perpendicular"]
    )
    fiducial_dh_dr1 = float(dr1["combined"]["DH_over_rd"]) / float(
        dr1["combined"]["alpha_parallel"]
    )
    for label, key in (("auto", "auto_alpha"), ("cross", "cross_alpha")):
        measurements.append(
            alpha_measurement(
                f"DESI_DR1_{label}",
                dr1_class,
                float(dr1["z"]),
                dr1[key],
                fiducial_dm_dr1,
                fiducial_dh_dr1,
            )
        )
    return measurements


def gaussian_validation(
    measurements: list[dict[str, object]],
    plamb_fit: dict[str, object],
    lcdm_fit: dict[str, object],
) -> pd.DataFrame:
    rows = []
    for measurement in measurements:
        observed = np.asarray(measurement["mean"], dtype=float)
        covariance = np.asarray(measurement["covariance"], dtype=float)
        precision = np.linalg.inv(covariance)
        eigenvalues, eigenvectors = np.linalg.eigh(covariance)
        for model in ("PLAMB_frozen_lower_z", "LambdaCDM_frozen_lower_z"):
            prediction = model_prediction(
                model,
                float(measurement["z"]),
                plamb_fit,
                lcdm_fit,
            )
            residual = prediction - observed
            chi2 = float(residual @ precision @ residual)
            eigen_pulls = eigenvectors.T @ residual / np.sqrt(eigenvalues)
            rows.append(
                {
                    "measurement_id": measurement["measurement_id"],
                    "classification": measurement["classification"],
                    "z": measurement["z"],
                    "model": model,
                    "observed_DM_over_rd": observed[0],
                    "observed_DH_over_rd": observed[1],
                    "predicted_DM_over_rd": prediction[0],
                    "predicted_DH_over_rd": prediction[1],
                    "residual_DM": residual[0],
                    "residual_DH": residual[1],
                    "marginal_pull_DM": residual[0] / math.sqrt(covariance[0, 0]),
                    "marginal_pull_DH": residual[1] / math.sqrt(covariance[1, 1]),
                    "principal_mode_pull_1": eigen_pulls[0],
                    "principal_mode_pull_2": eigen_pulls[1],
                    "chi2_2d": chi2,
                    "p_2d": float(chi2_distribution.sf(chi2, 2)),
                    "systematic_fraction_added": measurement[
                        "systematic_fraction_added"
                    ],
                    "covariance_minimum_eigenvalue": float(np.min(eigenvalues)),
                }
            )
    return pd.DataFrame(rows)


def load_native_grid(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    raw = np.loadtxt(path, dtype=float)
    dm_axis = np.unique(raw[:, 0])
    dh_axis = np.unique(raw[:, 1])
    frame = pd.DataFrame(raw, columns=["dm", "dh", "likelihood"])
    likelihood = (
        frame.pivot(index="dm", columns="dh", values="likelihood")
        .loc[dm_axis, dh_axis]
        .to_numpy(dtype=float)
    )
    if not np.all(np.isfinite(likelihood)) or np.max(likelihood) <= 0.0:
        raise ValueError(f"Invalid native likelihood grid: {path}")
    return dm_axis, dh_axis, likelihood / np.max(likelihood)


def likelihood_grid_score(
    dm_axis: np.ndarray,
    dh_axis: np.ndarray,
    likelihood: np.ndarray,
    prediction: np.ndarray,
) -> dict[str, float | bool]:
    log_likelihood = np.log(np.maximum(likelihood, np.finfo(float).tiny))
    interpolator = RegularGridInterpolator(
        (dm_axis, dh_axis),
        log_likelihood,
        bounds_error=False,
        fill_value=np.nan,
    )
    point_log_likelihood = float(interpolator(prediction).item())
    inside = bool(np.isfinite(point_log_likelihood))
    if not inside:
        return {
            "inside_grid": False,
            "delta_chi2_native": float("inf"),
            "wilks_p_2d": 0.0,
            "HPD_tail_mass": 0.0,
        }
    maximum = float(np.max(log_likelihood))
    delta_chi2 = -2.0 * (point_log_likelihood - maximum)
    point_likelihood = math.exp(point_log_likelihood - maximum)
    relative = np.exp(log_likelihood - maximum)
    tail_mass = float(np.sum(relative[relative <= point_likelihood]) / np.sum(relative))
    return {
        "inside_grid": True,
        "delta_chi2_native": delta_chi2,
        "wilks_p_2d": float(chi2_distribution.sf(delta_chi2, 2)),
        "HPD_tail_mass": tail_mass,
    }


def native_eboss_validation(
    config: dict[str, object],
    plamb_fit: dict[str, object],
    lcdm_fit: dict[str, object],
) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    eboss = config["eboss_dr16"]
    auto_path = DATA_DIR / str(eboss["auto_grid"])
    cross_path = DATA_DIR / str(eboss["cross_grid"])
    dm_auto, dh_auto, likelihood_auto = load_native_grid(auto_path)
    dm_cross, dh_cross, likelihood_cross = load_native_grid(cross_path)
    if not np.array_equal(dm_auto, dm_cross) or not np.array_equal(dh_auto, dh_cross):
        raise ValueError("eBOSS native auto and cross grids use different axes")
    combined_log = np.log(np.maximum(likelihood_auto, np.finfo(float).tiny))
    combined_log += np.log(np.maximum(likelihood_cross, np.finfo(float).tiny))
    combined = np.exp(combined_log - np.max(combined_log))
    grids = {
        "auto": likelihood_auto,
        "cross": likelihood_cross,
        "combined": combined,
    }
    rows = []
    z = float(eboss["z"])
    for component, likelihood in grids.items():
        best_index = np.unravel_index(np.argmax(likelihood), likelihood.shape)
        best_dm = float(dm_auto[best_index[0]])
        best_dh = float(dh_auto[best_index[1]])
        for model in ("PLAMB_frozen_lower_z", "LambdaCDM_frozen_lower_z"):
            prediction = model_prediction(model, z, plamb_fit, lcdm_fit)
            score = likelihood_grid_score(
                dm_auto, dh_auto, likelihood, prediction
            )
            rows.append(
                {
                    "measurement_id": f"eBOSS_DR16_native_{component}",
                    "classification": eboss["classification"],
                    "z": z,
                    "model": model,
                    "predicted_DM_over_rd": prediction[0],
                    "predicted_DH_over_rd": prediction[1],
                    "native_grid_best_DM": best_dm,
                    "native_grid_best_DH": best_dh,
                    **score,
                }
            )
    return pd.DataFrame(rows), {
        "dm_axis": dm_auto,
        "dh_axis": dh_auto,
        "auto": likelihood_auto,
        "cross": likelihood_cross,
        "combined": combined,
    }


def ellipse_for_covariance(
    axis: plt.Axes,
    mean: np.ndarray,
    covariance: np.ndarray,
    colour: str,
    label: str,
) -> None:
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    angle = math.degrees(math.atan2(eigenvectors[1, 0], eigenvectors[0, 0]))
    ellipse = Ellipse(
        xy=mean,
        width=2.0 * math.sqrt(eigenvalues[0]),
        height=2.0 * math.sqrt(eigenvalues[1]),
        angle=angle,
        facecolor="none",
        edgecolor=colour,
        linewidth=1.4,
        label=label,
    )
    axis.add_patch(ellipse)


def make_plot(
    path: Path,
    equivalence_table: pd.DataFrame,
    equivalence_summary: dict[str, object],
    measurements: list[dict[str, object]],
    plamb_fit: dict[str, object],
    lcdm_fit: dict[str, object],
    native_grids: dict[str, np.ndarray],
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(16.0, 4.9), constrained_layout=True)

    axes[0].plot(
        equivalence_table["z"],
        equivalence_table["pointwise_b_perpendicular"],
        "o-",
        color="#2b6cb0",
        label="Pointwise transverse",
    )
    axes[0].plot(
        equivalence_table["z"],
        equivalence_table["pointwise_b_parallel"],
        "s-",
        color="#b33a3a",
        label="Pointwise radial",
    )
    axes[0].axhline(
        float(equivalence_summary["b_perpendicular_equivalent"]),
        color="#2b6cb0",
        linestyle="--",
        label="Best constant transverse",
    )
    axes[0].axhline(
        float(equivalence_summary["b_parallel_equivalent"]),
        color="#b33a3a",
        linestyle="--",
        label="Best constant radial",
    )
    axes[0].set_xlabel("Redshift z")
    axes[0].set_ylabel("Effective power")
    axes[0].set_title("Lambda-CDM-equivalent PLAMB powers")
    axes[0].grid(alpha=0.25)
    axes[0].legend(fontsize=7)

    selected = {
        "DESI_DR2_combined_stat_sys": "#b33a3a",
        "DESI_DR2_auto_plus_0p3pct_theory_budget": "#7b4f9d",
        "DESI_DR2_cross_plus_0p3pct_theory_budget": "#d47a20",
        "DESI_DR2_multipoles_plus_0p3pct_theory_budget": "#287a5a",
        "DESI_DR1_combined": "#555555",
    }
    for measurement in measurements:
        measurement_id = str(measurement["measurement_id"])
        if measurement_id not in selected:
            continue
        ellipse_for_covariance(
            axes[1],
            np.asarray(measurement["mean"], dtype=float),
            np.asarray(measurement["covariance"], dtype=float),
            selected[measurement_id],
            measurement_id.replace("_", " "),
        )
    for model, marker, colour in (
        ("PLAMB_frozen_lower_z", "x", "black"),
        ("LambdaCDM_frozen_lower_z", "+", "#2b6cb0"),
    ):
        prediction = model_prediction(model, 2.33, plamb_fit, lcdm_fit)
        axes[1].scatter(
            prediction[0],
            prediction[1],
            marker=marker,
            s=90,
            linewidths=2.0,
            color=colour,
            label=model.replace("_", " "),
            zorder=5,
        )
    axes[1].set_xlabel("D_M / r_d")
    axes[1].set_ylabel("D_H / r_d")
    axes[1].set_title("Published 68% Gaussian summaries")
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=6, loc="best")

    dm_axis = native_grids["dm_axis"]
    dh_axis = native_grids["dh_axis"]
    combined = native_grids["combined"]
    contour_levels = np.sort(np.exp(-0.5 * np.asarray([6.18, 2.30])))
    axes[2].contour(
        dm_axis,
        dh_axis,
        combined.T,
        levels=contour_levels,
        colors=["#7393b3", "#24557a"],
        linewidths=[1.0, 1.6],
    )
    for model, marker, colour in (
        ("PLAMB_frozen_lower_z", "x", "black"),
        ("LambdaCDM_frozen_lower_z", "+", "#b33a3a"),
    ):
        prediction = model_prediction(model, 2.334, plamb_fit, lcdm_fit)
        axes[2].scatter(
            prediction[0],
            prediction[1],
            marker=marker,
            s=90,
            linewidths=2.0,
            color=colour,
            label=model.replace("_", " "),
        )
    axes[2].set_xlim(34.0, 42.5)
    axes[2].set_ylim(7.6, 9.8)
    axes[2].set_xlabel("D_M / r_d")
    axes[2].set_ylabel("D_H / r_d")
    axes[2].set_title("Native eBOSS DR16 auto x cross")
    axes[2].grid(alpha=0.25)
    axes[2].legend(fontsize=7)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_report(
    path: Path,
    p_profile: dict[str, float],
    lower_plamb_fit: dict[str, object],
    lower_lcdm_fit: dict[str, object],
    equivalence_summary: dict[str, object],
    effective_table: pd.DataFrame,
    gaussian_results: pd.DataFrame,
    native_results: pd.DataFrame,
    availability: pd.DataFrame,
    validations: pd.DataFrame,
    decisions: dict[str, bool | float | str],
) -> None:
    plamb_gaussian = gaussian_results[
        gaussian_results["model"] == "PLAMB_frozen_lower_z"
    ].copy()
    lcdm_gaussian = gaussian_results[
        gaussian_results["model"] == "LambdaCDM_frozen_lower_z"
    ].copy()
    plamb_native = native_results[
        native_results["model"] == "PLAMB_frozen_lower_z"
    ].copy()
    lines = [
        "# PLAMB / Lambda-CDM equivalence and Lyman-alpha validation audit",
        "",
        f"Generated: {datetime.now():%Y-%m-%d %H:%M:%S}",
        "",
        "## Decision",
        "",
        f"**{decisions['headline']}**",
        "",
        f"The analytic geometry audit gives `chi-squared={float(equivalence_summary['chi2_PLAMB_to_LambdaCDM_geometry_best']):.6f}` "
        f"for the best constant-power representation of the lower-redshift-trained Lambda-CDM curve, "
        f"with `b_perpendicular={float(equivalence_summary['b_perpendicular_equivalent']):.9f}` and "
        f"`b_parallel={float(equivalence_summary['b_parallel_equivalent']):.9f}`. The maximum absolute "
        f"observable difference is `{float(equivalence_summary['maximum_absolute_fractional_geometry_difference']):.6%}`.",
        "",
        f"This is a close approximation with geometry goodness `p={float(equivalence_summary['goodness_p_geometry']):.6f}`, "
        f"not an exact equivalence: the Lambda-CDM-equivalent powers are "
        f"`{float(equivalence_summary['equivalent_minus_data_parameter_distance_sigma']):.6f}` fitted-parameter "
        "standard deviations from the powers selected by the observed full DESI vector. The strict pre-run "
        "equivalence gate therefore fails.",
        "",
        f"The frozen PLAMB prediction remains rejected by the official compressed DESI DR2 combined "
        f"posterior at `p={float(decisions['PLAMB_DESI_DR2_combined_p']):.6g}`. The auto/cross split "
        f"localises the stronger stress to `{decisions['worst_DR2_split']}` with "
        f"`p={float(decisions['worst_DR2_split_p']):.6g}`.",
        "",
        f"The alternative DR2 multipole-method control gives `p={float(decisions['PLAMB_DR2_multipole_p']):.6g}`; "
        f"DESI DR1 gives `p={float(decisions['PLAMB_DESI_DR1_combined_p']):.6g}`; and the native "
        f"eBOSS combined grid gives HPD tail mass `{float(decisions['PLAMB_eBOSS_combined_HPD_tail']):.6g}` "
        f"(`Wilks p={float(decisions['PLAMB_eBOSS_combined_Wilks_p']):.6g}`). These are controls, "
        "not independent confirmations of a new radial/transverse degree of freedom.",
        "",
        "The multipole result shows that the DR2 rejection is stable to a substantially different compression "
        "and covariance method. It is not an independent-sky replication because it uses the same DR2 spectra. "
        "Conversely, the 1% rejection threshold is not crossed by DESI DR1 or native eBOSS DR16.",
        "",
        f"Relative to the equally frozen Lambda-CDM control, PLAMB has predictive "
        f"`Delta chi-squared={float(decisions['delta_chi2_PLAMB_minus_LambdaCDM_DESI_DR2_combined']):+.6f}` "
        f"for the DR2 combined summary and "
        f"`{float(decisions['delta_chi2_PLAMB_minus_LambdaCDM_DR2_multipoles']):+.6f}` for the multipole "
        f"summary. The external eBOSS native grid instead gives "
        f"`{float(decisions['delta_native_chi2_PLAMB_minus_LambdaCDM_eBOSS']):+.6f}`. The sign change "
        "across surveys is further evidence against interpreting the DR2 residual as a universal new degree "
        "of freedom.",
        "",
        "The native DESI DR2 correlation likelihood is not in the current public release. The DR2 rows "
        "below therefore use collaboration-published posterior compressions. The unavailable native rerun "
        "is retained as an explicit failed availability gate.",
        "",
        "## Frozen protocol",
        "",
        "The raw-MU supernova likelihood fixes `p`; all DESI points below `z=2` then fix the two powers "
        "and common BAO scale. No Lyman-alpha measurement is used to retune either model before it is scored.",
        "",
        "$$",
        "\\begin{aligned}",
        "k(z)                         &= \\frac{1+2.3z}{(1+z)^p}, \\\\",
        "\\chi(z)                    &= \\int_0^z k(u)\\,du, \\\\",
        "D_M(z)                       &= q\\,\\frac{\\chi(z)}{(1+z)^{b_\\perp}}, \\\\",
        "D_H(z)                       &= q\\,\\frac{k(z)}{(1+z)^{b_\\parallel}}.",
        "\\end{aligned}",
        "$$",
        "",
        "The comparator is flat Lambda-CDM with `Omega_m` fitted jointly to the same raw-MU supernova "
        "blocks and lower-redshift DESI BAO data, while its common BAO scale is profiled.",
        "",
        "## Training fits",
        "",
        pd.DataFrame(
            [
                {
                    "model": "PLAMB radial/transverse",
                    "p": lower_plamb_fit["p"],
                    "b_perpendicular": lower_plamb_fit["b_perpendicular"],
                    "b_parallel": lower_plamb_fit["b_parallel"],
                    "q": lower_plamb_fit["q_c0_over_H0rd"],
                    "omega_m": np.nan,
                },
                {
                    "model": "flat Lambda-CDM",
                    "p": np.nan,
                    "b_perpendicular": np.nan,
                    "b_parallel": np.nan,
                    "q": lower_lcdm_fit["q_c0_over_H0rd"],
                    "omega_m": lower_lcdm_fit["omega_m"],
                },
            ]
        ).to_markdown(index=False),
        "",
        "## Analytic equivalence",
        "",
        pd.DataFrame([equivalence_summary]).to_markdown(index=False),
        "",
        "The pointwise powers include the globally profiled scale. The local-derivative powers remove "
        "that scale and show whether a constant power is genuinely implied by the two background curves.",
        "",
        effective_table.to_markdown(index=False),
        "",
        "## DESI compressed validation",
        "",
        "### Frozen PLAMB",
        "",
        plamb_gaussian.to_markdown(index=False),
        "",
        "### Frozen Lambda-CDM control",
        "",
        lcdm_gaussian.to_markdown(index=False),
        "",
        "The 0.3% diagonal theoretical BAO-shift budget is added to the DR2 auto/cross and multipole "
        "controls where stated. The published DR2 combined covariance already includes it.",
        "",
        "## Native eBOSS DR16 likelihood",
        "",
        native_results.to_markdown(index=False),
        "",
        "The eBOSS release states that its native auto- and cross-correlation likelihood grids may be "
        "treated as independent. The combined row uses their product. `HPD_tail_mass` is the integrated "
        "grid likelihood at density no greater than the frozen prediction; `Wilks_p_2d` is supplied only "
        "as an asymptotic likelihood-ratio readout.",
        "",
        "## Product availability",
        "",
        availability.to_markdown(index=False),
        "",
        "## Gates",
        "",
    ]
    for name, value in decisions.items():
        if name != "headline":
            lines.append(f"- `{name}`: **{value}**")
    lines.extend(
        [
            "",
            "## Validation checks",
            "",
            validations.to_markdown(index=False),
            "",
            "## Interpretation",
            "",
            "A successful constant-power reconstruction of Lambda-CDM does not identify a new symmetry. "
            "It shows that the two fitted powers can encode an already familiar distance geometry. The "
            "DESI DR2 split then tests whether the remaining high-redshift residual is common to the forest "
            "auto-correlation and quasar cross-correlation or concentrated in one tracer construction.",
            "",
            "The alternative multipole analysis is a valuable method and covariance control, but it uses "
            "the same DR2 spectra. DESI DR1 is an earlier, partly overlapping survey control; eBOSS DR16 is "
            "the external spectroscopic control. These distinctions are preserved in every output table.",
            "",
            "The auto/cross readout is a localisation diagnostic, not a formal test of the difference between "
            "the two split estimates: their shared spectra and cross-covariance prevent treating those p-values "
            "as independent. It nevertheless identifies the Ly-alpha x quasar construction as the immediate "
            "place to inspect redshift errors, proximity modelling and quasar-related nuisance priors.",
            "",
            "No SU(2), anisotropic metric or evolving-ruler claim is authorised by this audit. A physical "
            "claim would additionally require the native DR2 likelihood, independent high-redshift "
            "replication, and an action-level source for the radial/transverse split.",
            "",
            "## Primary sources",
            "",
            "- DESI DR2 baseline: https://arxiv.org/abs/2503.14739",
            "- DESI DR2 figure data: https://zenodo.org/records/15690869",
            "- DESI DR2 multipoles: https://arxiv.org/abs/2603.04281",
            "- DESI DR1: https://arxiv.org/abs/2404.03001",
            "- eBOSS DR16: https://arxiv.org/abs/2007.08995",
            "- SDSS DR16 likelihood release: https://svn.sdss.org/public/data/eboss/DR16cosmo/tags/v1_0_1/likelihoods/BAO-only/",
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
        "analysis_id": "plamb_lcdm_lya_validation_2026-07-19",
        "status": "outcome_aware_followup_with_frozen_lower_redshift_training",
        "training": {
            "p": "released-total raw-MU supernova profile only",
            "radial_transverse_powers_and_q": "DESI DR2 BAO rows with z < 2 only",
            "LambdaCDM_Omega_m_and_q": "same raw-MU supernova blocks plus DESI rows with z < 2",
        },
        "tests": [
            "analytic constant-power representation of lower-z-trained Lambda-CDM geometry",
            "DESI DR2 compressed combined and auto/cross posterior summaries",
            "same-spectra alternative DR2 multipole-method summary",
            "earlier DESI DR1 combined and auto/cross summaries",
            "external eBOSS DR16 native auto/cross likelihood grids",
        ],
        "gates": {
            "practical_geometry_equivalence": "chi2 <= 2.30, parameter distance <= 2 sigma, and maximum fractional mismatch <= 1.5%",
            "individual_validation": "two-dimensional p or HPD tail >= 0.01",
            "promotion_validation": "two-dimensional p or HPD tail >= 0.05",
            "native_DR2": "native public likelihood required for native-likelihood claim",
        },
        "claim_boundary": "method and external-data validation only; no symmetry, SU(2), metric or ruler identification",
    }
    protocol_path = args.outdir / f"plamb_lcdm_lya_protocol_{DATE_TAG}.json"
    protocol_path.write_text(json.dumps(protocol, indent=2), encoding="utf-8")

    config = json.loads(MEASUREMENT_CONFIG.read_text(encoding="utf-8"))
    variants, load_metadata = rawmu.load_blocks(MIN_SURVEY_N)
    blocks = variants["released_total_primary"]
    data, covariance = closure.load_bao(args.mean, args.cov)
    lower_mask = data["z"].to_numpy(dtype=float) < 2.0
    lower_indices = np.flatnonzero(lower_mask)
    lower_data = data.iloc[lower_indices].reset_index(drop=True)
    lower_covariance = covariance[np.ix_(lower_indices, lower_indices)]

    p_profile = radial.fit_sn_profile(blocks)
    p_mle = float(p_profile["p_mle"])
    lower_plamb_fit = radial.fit_bao_at_fixed_p(
        lower_data, lower_covariance, "separate", p_mle
    )
    full_plamb_fit = radial.fit_bao_at_fixed_p(
        data, covariance, "separate", p_mle
    )
    lower_lcdm_fit = radial.fit_lcdm_joint(blocks, lower_data, lower_covariance)

    equivalence_summary, effective_table, equivalence_vector = fit_equivalence(
        data,
        covariance,
        p_mle,
        lower_lcdm_fit,
        full_plamb_fit,
        lower_plamb_fit,
    )
    measurements = build_gaussian_measurements(config)
    gaussian_results = gaussian_validation(
        measurements, lower_plamb_fit, lower_lcdm_fit
    )
    native_results, native_grids = native_eboss_validation(
        config, lower_plamb_fit, lower_lcdm_fit
    )

    availability = pd.DataFrame(
        [
            {
                "product": "DESI DR2 baseline native 9306-vector likelihood",
                "publicly_available": False,
                "used": False,
                "role": "pending native rerun",
            },
            {
                "product": "DESI DR2 combined/auto/cross posterior compression",
                "publicly_available": True,
                "used": True,
                "role": "primary compressed validation",
            },
            {
                "product": "DESI DR2 alternative multipole likelihood/data",
                "publicly_available": False,
                "used": False,
                "role": "published summary method control",
            },
            {
                "product": "DESI DR1 correlation data and covariance",
                "publicly_available": True,
                "used": False,
                "role": "published summary; native Vega refit not performed",
            },
            {
                "product": "eBOSS DR16 auto/cross likelihood grids",
                "publicly_available": True,
                "used": True,
                "role": "external native non-Gaussian validation",
            },
        ]
    )

    def gaussian_row(measurement_id: str, model: str = "PLAMB_frozen_lower_z") -> pd.Series:
        selected = gaussian_results[
            (gaussian_results["measurement_id"] == measurement_id)
            & (gaussian_results["model"] == model)
        ]
        if len(selected) != 1:
            raise ValueError(f"Expected one row for {measurement_id} and {model}")
        return selected.iloc[0]

    def native_row(measurement_id: str, model: str = "PLAMB_frozen_lower_z") -> pd.Series:
        selected = native_results[
            (native_results["measurement_id"] == measurement_id)
            & (native_results["model"] == model)
        ]
        if len(selected) != 1:
            raise ValueError(f"Expected one row for {measurement_id} and {model}")
        return selected.iloc[0]

    dr2_combined = gaussian_row("DESI_DR2_combined_stat_sys")
    dr2_combined_lcdm = gaussian_row(
        "DESI_DR2_combined_stat_sys", "LambdaCDM_frozen_lower_z"
    )
    dr2_auto = gaussian_row("DESI_DR2_auto_plus_0p3pct_theory_budget")
    dr2_cross = gaussian_row("DESI_DR2_cross_plus_0p3pct_theory_budget")
    split_rows = pd.DataFrame([dr2_auto, dr2_cross])
    worst_split = split_rows.loc[split_rows["p_2d"].idxmin()]
    multipole = gaussian_row("DESI_DR2_multipoles_plus_0p3pct_theory_budget")
    multipole_lcdm = gaussian_row(
        "DESI_DR2_multipoles_plus_0p3pct_theory_budget",
        "LambdaCDM_frozen_lower_z",
    )
    dr1_combined = gaussian_row("DESI_DR1_combined")
    dr1_combined_lcdm = gaussian_row(
        "DESI_DR1_combined", "LambdaCDM_frozen_lower_z"
    )
    eboss_combined = native_row("eBOSS_DR16_native_combined")
    eboss_combined_lcdm = native_row(
        "eBOSS_DR16_native_combined", "LambdaCDM_frozen_lower_z"
    )

    practical_equivalence = bool(
        float(equivalence_summary["chi2_PLAMB_to_LambdaCDM_geometry_best"]) <= 2.30
        and float(
            equivalence_summary["equivalent_minus_data_parameter_distance_sigma"]
        )
        <= 2.0
        and float(
            equivalence_summary["maximum_absolute_fractional_geometry_difference"]
        )
        <= 0.015
    )
    independent_controls_pass_0p01 = bool(
        float(dr1_combined["p_2d"]) >= 0.01
        and float(eboss_combined["HPD_tail_mass"]) >= 0.01
    )
    all_compressed_controls_pass_0p01 = bool(
        float(dr2_combined["p_2d"]) >= 0.01
        and float(dr2_auto["p_2d"]) >= 0.01
        and float(dr2_cross["p_2d"]) >= 0.01
        and float(multipole["p_2d"]) >= 0.01
        and independent_controls_pass_0p01
    )
    decisions: dict[str, bool | float | str] = {
        "headline": (
            "THE POWERS CLOSELY APPROXIMATE LAMBDA-CDM BUT MISS THE STRICT EQUIVALENCE GATE; THE DR2 FAILURE IS CROSS-CORRELATION-LED AND SURVIVES THE MULTIPOLE CONTROL."
            if not practical_equivalence
            and float(dr2_combined["p_2d"]) < 0.01
            and float(multipole["p_2d"]) < 0.01
            else "THE HIGH-REDSHIFT VALIDATION DOES NOT AUTHORISE A PHYSICAL RADIAL/TRANSVERSE CLAIM."
        ),
        "analytic_practical_LambdaCDM_equivalence": practical_equivalence,
        "analytic_geometry_goodness_p": float(
            equivalence_summary["goodness_p_geometry"]
        ),
        "analytic_maximum_fractional_geometry_difference": float(
            equivalence_summary["maximum_absolute_fractional_geometry_difference"]
        ),
        "analytic_equivalent_minus_data_parameter_distance_sigma": float(
            equivalence_summary["equivalent_minus_data_parameter_distance_sigma"]
        ),
        "PLAMB_DESI_DR2_combined_p": float(dr2_combined["p_2d"]),
        "delta_chi2_PLAMB_minus_LambdaCDM_DESI_DR2_combined": float(
            dr2_combined["chi2_2d"] - dr2_combined_lcdm["chi2_2d"]
        ),
        "PLAMB_DESI_DR2_auto_p": float(dr2_auto["p_2d"]),
        "PLAMB_DESI_DR2_cross_p": float(dr2_cross["p_2d"]),
        "worst_DR2_split": str(worst_split["measurement_id"]),
        "worst_DR2_split_p": float(worst_split["p_2d"]),
        "PLAMB_DR2_multipole_p": float(multipole["p_2d"]),
        "delta_chi2_PLAMB_minus_LambdaCDM_DR2_multipoles": float(
            multipole["chi2_2d"] - multipole_lcdm["chi2_2d"]
        ),
        "PLAMB_DESI_DR1_combined_p": float(dr1_combined["p_2d"]),
        "delta_chi2_PLAMB_minus_LambdaCDM_DESI_DR1": float(
            dr1_combined["chi2_2d"] - dr1_combined_lcdm["chi2_2d"]
        ),
        "PLAMB_eBOSS_combined_HPD_tail": float(eboss_combined["HPD_tail_mass"]),
        "PLAMB_eBOSS_combined_Wilks_p": float(eboss_combined["wilks_p_2d"]),
        "delta_native_chi2_PLAMB_minus_LambdaCDM_eBOSS": float(
            eboss_combined["delta_chi2_native"]
            - eboss_combined_lcdm["delta_chi2_native"]
        ),
        "independent_DR1_and_eBOSS_controls_pass_0p01": independent_controls_pass_0p01,
        "DR2_failure_survives_same_spectra_multipole_method": bool(
            float(dr2_combined["p_2d"]) < 0.01
            and float(multipole["p_2d"]) < 0.01
        ),
        "DR2_cross_fails_while_auto_passes_0p01": bool(
            float(dr2_cross["p_2d"]) < 0.01
            and float(dr2_auto["p_2d"]) >= 0.01
        ),
        "DR2_failure_replicated_by_independent_survey_at_0p01": bool(
            float(dr1_combined["p_2d"]) < 0.01
            or float(eboss_combined["HPD_tail_mass"]) < 0.01
        ),
        "all_available_compressed_or_native_controls_pass_0p01": all_compressed_controls_pass_0p01,
        "native_DESI_DR2_likelihood_publicly_available": False,
        "native_DESI_DR2_likelihood_rerun_complete": False,
        "radial_transverse_adapter_validated": False,
        "physical_degree_identified": False,
        "further_physical_sampling_authorised": False,
    }

    plamb_z233 = model_prediction(
        "PLAMB_frozen_lower_z", 2.33, lower_plamb_fit, lower_lcdm_fit
    )
    validation_rows = [
        {
            "check": "SN p regression",
            "value": p_mle,
            "target": radial.EXPECTED["p_sn"],
            "tolerance": 5e-5,
            "passed": abs(p_mle - radial.EXPECTED["p_sn"]) <= 5e-5,
        },
        {
            "check": "lower-z b_perpendicular regression",
            "value": float(lower_plamb_fit["b_perpendicular"]),
            "target": 0.992066,
            "tolerance": 5e-5,
            "passed": abs(float(lower_plamb_fit["b_perpendicular"]) - 0.992066)
            <= 5e-5,
        },
        {
            "check": "lower-z b_parallel regression",
            "value": float(lower_plamb_fit["b_parallel"]),
            "target": 1.752423,
            "tolerance": 5e-5,
            "passed": abs(float(lower_plamb_fit["b_parallel"]) - 1.752423)
            <= 5e-5,
        },
        {
            "check": "frozen z2.33 DM regression",
            "value": float(plamb_z233[0]),
            "target": 40.0853,
            "tolerance": 5e-4,
            "passed": abs(float(plamb_z233[0]) - 40.0853) <= 5e-4,
        },
        {
            "check": "frozen z2.33 DH regression",
            "value": float(plamb_z233[1]),
            "target": 8.80755,
            "tolerance": 5e-4,
            "passed": abs(float(plamb_z233[1]) - 8.80755) <= 5e-4,
        },
        {
            "check": "all Gaussian covariances positive definite",
            "value": float(gaussian_results["covariance_minimum_eigenvalue"].min()),
            "target": 0.0,
            "tolerance": 0.0,
            "passed": bool(
                (gaussian_results["covariance_minimum_eigenvalue"] > 0.0).all()
            ),
        },
        {
            "check": "native eBOSS predictions inside both grids",
            "value": float(native_results["inside_grid"].astype(int).mean()),
            "target": 1.0,
            "tolerance": 0.0,
            "passed": bool(native_results["inside_grid"].all()),
        },
        {
            "check": "equivalence optimisation succeeded",
            "value": float(bool(equivalence_summary["optimisation_success"])),
            "target": 1.0,
            "tolerance": 0.0,
            "passed": bool(equivalence_summary["optimisation_success"]),
        },
    ]
    validations = pd.DataFrame(validation_rows)

    prefix = "plamb_lcdm_lya"
    paths = {
        "equivalence_summary": args.outdir
        / f"{prefix}_equivalence_summary_{DATE_TAG}.json",
        "effective_powers": args.outdir
        / f"{prefix}_effective_powers_{DATE_TAG}.csv",
        "equivalence_vector": args.outdir
        / f"{prefix}_equivalence_vector_{DATE_TAG}.csv",
        "gaussian_validation": args.outdir
        / f"{prefix}_gaussian_validation_{DATE_TAG}.csv",
        "native_eboss": args.outdir
        / f"{prefix}_native_eboss_validation_{DATE_TAG}.csv",
        "availability": args.outdir
        / f"{prefix}_product_availability_{DATE_TAG}.csv",
        "validation": args.outdir / f"{prefix}_validation_checks_{DATE_TAG}.csv",
        "decision": args.outdir / f"{prefix}_decision_{DATE_TAG}.json",
        "plot": args.outdir / f"{prefix}_readout_{DATE_TAG}.png",
        "report": args.outdir / f"{prefix}_report_{DATE_TAG}.md",
        "summary": args.outdir / f"{prefix}_summary_{DATE_TAG}.json",
    }
    paths["equivalence_summary"].write_text(
        json.dumps(equivalence_summary, indent=2), encoding="utf-8"
    )
    effective_table.to_csv(paths["effective_powers"], index=False)
    equivalence_vector.to_csv(paths["equivalence_vector"], index=False)
    gaussian_results.to_csv(paths["gaussian_validation"], index=False)
    native_results.to_csv(paths["native_eboss"], index=False)
    availability.to_csv(paths["availability"], index=False)
    validations.to_csv(paths["validation"], index=False)
    paths["decision"].write_text(json.dumps(decisions, indent=2), encoding="utf-8")
    make_plot(
        paths["plot"],
        effective_table,
        equivalence_summary,
        measurements,
        lower_plamb_fit,
        lower_lcdm_fit,
        native_grids,
    )
    write_report(
        paths["report"],
        p_profile,
        lower_plamb_fit,
        lower_lcdm_fit,
        equivalence_summary,
        effective_table,
        gaussian_results,
        native_results,
        availability,
        validations,
        decisions,
    )
    summary = {
        "analysis_date": DATE_TAG,
        "status": protocol["status"],
        "load_metadata": load_metadata,
        "p_profile": p_profile,
        "lower_z_PLAMB_fit": lower_plamb_fit,
        "lower_z_LambdaCDM_fit": lower_lcdm_fit,
        "equivalence": equivalence_summary,
        "decisions": decisions,
        "validation_checks_pass": bool(validations["passed"].all()),
    }
    paths["summary"].write_text(json.dumps(summary, indent=2), encoding="utf-8")

    sources = [
        SCRIPT_PATH,
        MEASUREMENT_CONFIG,
        DATA_DIR / "desi_dr2_lya_alphas.txt",
        DATA_DIR / "desi_dr2_lya_data_splits.txt",
        DATA_DIR / str(config["eboss_dr16"]["auto_grid"]),
        DATA_DIR / str(config["eboss_dr16"]["cross_grid"]),
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

    print(
        "Equivalence: "
        f"chi2={float(equivalence_summary['chi2_PLAMB_to_LambdaCDM_geometry_best']):.6f} "
        f"b_perp={float(equivalence_summary['b_perpendicular_equivalent']):.9f} "
        f"b_parallel={float(equivalence_summary['b_parallel_equivalent']):.9f}"
    )
    print(
        f"DR2 combined PLAMB p={float(dr2_combined['p_2d']):.8g}; "
        f"worst split={worst_split['measurement_id']} p={float(worst_split['p_2d']):.8g}"
    )
    print(
        f"DR1 PLAMB p={float(dr1_combined['p_2d']):.8g}; "
        f"eBOSS combined HPD={float(eboss_combined['HPD_tail_mass']):.8g}"
    )
    print(f"Saved report: {paths['report']}")


if __name__ == "__main__":
    main()
