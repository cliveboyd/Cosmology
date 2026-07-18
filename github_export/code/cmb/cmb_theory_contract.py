#!/usr/bin/env python3
"""Typed contracts shared by the registered CMB analysis branches."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Mapping

import numpy as np
from numpy.typing import NDArray


FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class TheoryCapabilities:
    background: bool = False
    thermodynamics: bool = False
    scalar_perturbations: bool = False
    tensor_perturbations: bool = False
    parity_odd_spectra: bool = False
    spectral_distortions: bool = False


@dataclass(frozen=True)
class BackgroundPrediction:
    z: FloatArray
    hubble_km_s_mpc: FloatArray
    comoving_distance_mpc: FloatArray
    z_star: float
    sound_horizon_star_mpc: float
    omega_b_h2: float
    omega_m: float

    def validate(self) -> None:
        if self.z.ndim != 1 or len(self.z) < 2:
            raise ValueError("Background redshift grid must be one-dimensional")
        if not np.all(np.diff(self.z) > 0):
            raise ValueError("Background redshift grid must be strictly increasing")
        for name, values in (
            ("hubble", self.hubble_km_s_mpc),
            ("comoving distance", self.comoving_distance_mpc),
        ):
            if values.shape != self.z.shape or not np.all(np.isfinite(values)):
                raise ValueError(f"Invalid {name} prediction")
        if np.any(self.hubble_km_s_mpc <= 0) or np.any(self.comoving_distance_mpc < 0):
            raise ValueError("Background predictions violate positivity")
        if not (
            self.z_star > 0
            and self.sound_horizon_star_mpc > 0
            and self.omega_b_h2 > 0
            and self.omega_m > 0
        ):
            raise ValueError("Recombination and baryon quantities must be positive")
        if not self.z[0] <= self.z_star <= self.z[-1]:
            raise ValueError("Background grid does not bracket recombination redshift")


@dataclass(frozen=True)
class CmbSpectra:
    ell: NDArray[np.int64]
    cl: Mapping[str, FloatArray]
    units: str = "C_ell_uK2"

    def validate(self, required: tuple[str, ...]) -> None:
        if self.ell.ndim != 1 or len(self.ell) < 3 or self.ell[0] != 0:
            raise ValueError("CMB multipoles must be a one-dimensional grid beginning at zero")
        if not np.array_equal(self.ell, np.arange(len(self.ell))):
            raise ValueError("CMB multipoles must be contiguous")
        for spectrum in required:
            if spectrum not in self.cl:
                raise ValueError(f"Missing required CMB spectrum: {spectrum}")
            values = np.asarray(self.cl[spectrum], dtype=float)
            if values.shape != self.ell.shape or not np.all(np.isfinite(values)):
                raise ValueError(f"Invalid CMB spectrum: {spectrum}")


@dataclass(frozen=True)
class SpectralDistortion:
    frequency_ghz: FloatArray
    delta_intensity_kjy_sr: FloatArray

    def validate(self) -> None:
        if self.frequency_ghz.ndim != 1 or len(self.frequency_ghz) < 2:
            raise ValueError("Distortion frequency grid must be one-dimensional")
        if not np.all(np.diff(self.frequency_ghz) > 0):
            raise ValueError("Distortion frequency grid must be strictly increasing")
        if self.delta_intensity_kjy_sr.shape != self.frequency_ghz.shape:
            raise ValueError("Distortion intensity and frequency grids differ")
        if not np.all(np.isfinite(self.delta_intensity_kjy_sr)):
            raise ValueError("Distortion prediction contains non-finite values")


class CmbTheoryAdapter(ABC):
    """Theory boundary; branch code must not infer missing physics silently."""

    capabilities: TheoryCapabilities

    @abstractmethod
    def background(self, parameters: Mapping[str, float], z: FloatArray) -> BackgroundPrediction:
        raise NotImplementedError

    def spectra(self, parameters: Mapping[str, float], lmax: int) -> CmbSpectra:
        raise NotImplementedError("This theory has no registered perturbation implementation")

    def distortion(self, parameters: Mapping[str, float], frequency_ghz: FloatArray) -> SpectralDistortion:
        raise NotImplementedError("This theory has no registered photon-thermalisation implementation")


def planck_distance_vector(prediction: BackgroundPrediction) -> FloatArray:
    """Return (R, l_A, omega_b h^2) for a standard early-physics branch."""
    prediction.validate()
    c_km_s = 299792.458
    dm_star = float(np.interp(prediction.z_star, prediction.z, prediction.comoving_distance_mpc))
    h0 = float(prediction.hubble_km_s_mpc[0])
    acoustic_scale = np.pi * dm_star / prediction.sound_horizon_star_mpc
    shift = np.sqrt(prediction.omega_m) * h0 * dm_star / c_km_s
    return np.array([shift, acoustic_scale, prediction.omega_b_h2], dtype=float)
