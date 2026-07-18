#!/usr/bin/env python3
"""Unit-explicit parity and CMB-lensing spectrum transformations."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


FloatArray = NDArray[np.float64]


def phi_phi_to_kappa_kappa(ell: FloatArray, cl_phi_phi: FloatArray) -> FloatArray:
    multipole = np.asarray(ell, dtype=float)
    potential = np.asarray(cl_phi_phi, dtype=float)
    if multipole.shape != potential.shape:
        raise ValueError("Multipole and lensing-potential spectra differ in shape")
    if np.any(multipole < 0) or not np.all(np.isfinite(potential)):
        raise ValueError("Invalid lensing-potential spectrum")
    factor = np.square(multipole * (multipole + 1.0)) / 4.0
    return factor * potential


def rotate_even_parity_spectra(
    alpha_rad: float,
    cl_te: FloatArray,
    cl_ee: FloatArray,
    cl_bb: FloatArray,
) -> dict[str, FloatArray]:
    """Rotate zero-intrinsic-TB/EB spectra by a uniform birefringence angle."""
    te = np.asarray(cl_te, dtype=float)
    ee = np.asarray(cl_ee, dtype=float)
    bb = np.asarray(cl_bb, dtype=float)
    if te.shape != ee.shape or te.shape != bb.shape:
        raise ValueError("TE, EE and BB spectra differ in shape")
    angle2 = 2.0 * float(alpha_rad)
    cosine = np.cos(angle2)
    sine = np.sin(angle2)
    return {
        "TE": te * cosine,
        "TB": te * sine,
        "EE": ee * cosine * cosine + bb * sine * sine,
        "BB": bb * cosine * cosine + ee * sine * sine,
        "EB": 0.5 * (ee - bb) * np.sin(2.0 * angle2),
    }


def tensor_chirality(power_right: FloatArray, power_left: FloatArray) -> FloatArray:
    right = np.asarray(power_right, dtype=float)
    left = np.asarray(power_left, dtype=float)
    if right.shape != left.shape or np.any(right < 0) or np.any(left < 0):
        raise ValueError("Invalid right/left tensor power")
    total = right + left
    return np.divide(right - left, total, out=np.zeros_like(total), where=total > 0)
