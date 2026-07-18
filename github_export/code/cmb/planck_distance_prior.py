#!/usr/bin/env python3
"""Load and evaluate the registered Planck three-variable compression."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray


FloatArray = NDArray[np.float64]
C_KM_S = 299792.458


@dataclass(frozen=True)
class PlanckDistancePrior:
    labels: tuple[str, str, str]
    mean: FloatArray
    covariance: FloatArray
    precision: FloatArray
    claim_boundary: str

    @classmethod
    def load(cls, path: Path) -> "PlanckDistancePrior":
        payload = json.loads(path.read_text(encoding="utf-8"))
        labels = tuple(payload["vector_order"])
        if labels != ("R", "l_A", "omega_b_h2"):
            raise ValueError(f"Unexpected Planck vector order: {labels}")
        mean = np.array([payload["means"][label] for label in labels], dtype=float)
        covariance = np.asarray(payload["cov"], dtype=float)
        if covariance.shape != (3, 3) or not np.all(np.isfinite(covariance)):
            raise ValueError("Invalid Planck covariance")
        covariance = 0.5 * (covariance + covariance.T)
        if not np.all(np.linalg.eigvalsh(covariance) > 0):
            raise ValueError("Planck covariance is not positive definite")
        return cls(
            labels=labels,
            mean=mean,
            covariance=covariance,
            precision=np.linalg.inv(covariance),
            claim_boundary=payload["claim_boundary"],
        )

    def chi_square(self, prediction: FloatArray) -> float:
        vector = np.asarray(prediction, dtype=float)
        if vector.shape != (3,) or not np.all(np.isfinite(vector)):
            return float("inf")
        residual = vector - self.mean
        return float(residual @ self.precision @ residual)


def vector_from_camb_results(results: object) -> FloatArray:
    """Calculate (R, l_A, omega_b h^2) from one CAMB result."""
    derived = results.get_derived_params()
    params = results.Params
    theta_star = float(derived["thetastar"]) / 100.0
    r_star = float(derived["rstar"])
    dm_star_mpc = r_star / theta_star
    h = float(params.H0) / 100.0
    omega_m = float(params.ombh2 + params.omch2 + params.omnuh2) / (h * h)
    shift = np.sqrt(omega_m) * float(params.H0) * dm_star_mpc / C_KM_S
    return np.array([shift, np.pi / theta_star, float(params.ombh2)], dtype=float)
