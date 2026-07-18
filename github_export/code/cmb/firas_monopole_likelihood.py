#!/usr/bin/env python3
"""COBE/FIRAS monopole ingestion and covariance-enforced linear likelihood."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray


FloatArray = NDArray[np.float64]
H_PLANCK = 6.62607015e-34
K_BOLTZMANN = 1.380649e-23
C_M_S = 299792458.0
T_CMB_K = 2.725
BETA_MU = 2.1923


@dataclass(frozen=True)
class FirasMonopoleData:
    frequency_ghz: FloatArray
    monopole_mjy_sr: FloatArray
    residual_kjy_sr: FloatArray
    sigma_kjy_sr: FloatArray
    galaxy_kjy_sr: FloatArray

    @classmethod
    def load(cls, path: Path) -> "FirasMonopoleData":
        values = np.loadtxt(path, comments="#", dtype=float)
        if values.shape != (43, 5):
            raise ValueError(f"Expected FIRAS shape (43, 5), found {values.shape}")
        # The released coordinate is a wavenumber in cm^-1; convert cm^-1 to m^-1.
        frequency_ghz = values[:, 0] * 100.0 * C_M_S / 1.0e9
        result = cls(
            frequency_ghz=frequency_ghz,
            monopole_mjy_sr=values[:, 1],
            residual_kjy_sr=values[:, 2],
            sigma_kjy_sr=values[:, 3],
            galaxy_kjy_sr=values[:, 4],
        )
        if not np.all(np.diff(result.frequency_ghz) > 0):
            raise ValueError("FIRAS frequencies are not strictly increasing")
        if not np.all(np.isfinite(values)) or np.any(result.sigma_kjy_sr <= 0):
            raise ValueError("FIRAS table contains invalid values")
        return result


def spectral_templates_kjy_sr(frequency_ghz: FloatArray, temperature_k: float = T_CMB_K) -> dict[str, FloatArray]:
    nu = np.asarray(frequency_ghz, dtype=float) * 1.0e9
    x = H_PLANCK * nu / (K_BOLTZMANN * temperature_k)
    exponential = np.exp(x)
    occupation_factor = exponential / np.square(exponential - 1.0)
    intensity_scale_kjy_sr = 2.0 * np.power(K_BOLTZMANN * temperature_k, 3) / np.square(H_PLANCK * C_M_S) * 1.0e23
    temperature = intensity_scale_kjy_sr * np.power(x, 4) * occupation_factor / temperature_k
    mu = intensity_scale_kjy_sr * np.power(x, 3) * occupation_factor * (x / BETA_MU - 1.0)
    y = intensity_scale_kjy_sr * np.power(x, 4) * occupation_factor * (
        x * (exponential + 1.0) / (exponential - 1.0) - 4.0
    )
    return {"delta_temperature_K": temperature, "mu": mu, "y": y}


def load_covariance(
    data: FirasMonopoleData,
    covariance_path: Path | None,
    *,
    allow_diagonal_diagnostic: bool = False,
) -> tuple[FloatArray, str]:
    if covariance_path is None:
        if not allow_diagonal_diagnostic:
            raise RuntimeError(
                "A registered FIRAS covariance is required for a primary fit; "
                "set allow_diagonal_diagnostic only for non-gating software checks"
            )
        return np.diag(np.square(data.sigma_kjy_sr)), "diagonal_diagnostic_only"
    covariance = np.load(covariance_path)
    if covariance.shape != (43, 43):
        raise ValueError(f"Expected FIRAS covariance shape (43, 43), found {covariance.shape}")
    covariance = 0.5 * (covariance + covariance.T)
    if not np.all(np.isfinite(covariance)) or not np.all(np.linalg.eigvalsh(covariance) > 0):
        raise ValueError("FIRAS covariance must be finite and positive definite")
    return covariance, "registered_full_covariance"


def fit_linear_templates(
    data: FirasMonopoleData,
    covariance: FloatArray,
    *,
    injection_template_kjy_sr: FloatArray | None = None,
) -> dict[str, object]:
    templates = spectral_templates_kjy_sr(data.frequency_ghz)
    names = ["delta_temperature_K", "mu", "y", "galaxy_amplitude"]
    columns = [templates["delta_temperature_K"], templates["mu"], templates["y"], data.galaxy_kjy_sr]
    if injection_template_kjy_sr is not None:
        injection = np.asarray(injection_template_kjy_sr, dtype=float)
        if injection.shape != data.frequency_ghz.shape or not np.all(np.isfinite(injection)):
            raise ValueError("Invalid physical injection template")
        names.append("physical_injection_amplitude")
        columns.append(injection)
    design = np.column_stack(columns)
    precision = np.linalg.inv(covariance)
    fisher = design.T @ precision @ design
    parameter_covariance = np.linalg.inv(fisher)
    estimate = parameter_covariance @ design.T @ precision @ data.residual_kjy_sr
    residual = data.residual_kjy_sr - design @ estimate
    return {
        "parameter_names": names,
        "estimates": estimate,
        "parameter_covariance": parameter_covariance,
        "chi_square": float(residual @ precision @ residual),
        "degrees_of_freedom": int(len(residual) - len(estimate)),
    }
