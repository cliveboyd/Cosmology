#!/usr/bin/env python3
"""Analytic PLAMB clock-law path distances with explicit named branches."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ClockLaw:
    key: str
    label: str
    gamma_c: float
    redshift_rate_power: float
    description: str


PETER_LINEAR_REDSHIFT = ClockLaw(
    key="peter_linear_redshift",
    label="Peter linear-redshift clock law (p=0)",
    gamma_c=1.0,
    redshift_rate_power=0.0,
    description=(
        "c(z)/c0=1+z and |dz/dT|=H0; gives "
        "D_path=(c0/H0) z(1+z/2)."
    ),
)

FRACTIONAL_CLOCK = ClockLaw(
    key="fractional_clock",
    label="Constant fractional clock law (p=1)",
    gamma_c=1.0,
    redshift_rate_power=1.0,
    description=(
        "c(z)/c0=1+z and |dz/dT|=H0(1+z); gives "
        "D_path=(c0/H0) z."
    ),
)

CONSTANT_C_LINEAR_REDSHIFT = ClockLaw(
    key="constant_c_linear_redshift",
    label="Constant-c linear-redshift control",
    gamma_c=0.0,
    redshift_rate_power=0.0,
    description="c(z)=c0 and |dz/dT|=H0; gives D_path=(c0/H0) z.",
)

CLOCK_LAWS = {
    law.key: law
    for law in (
        PETER_LINEAR_REDSHIFT,
        FRACTIONAL_CLOCK,
        CONSTANT_C_LINEAR_REDSHIFT,
    )
}


def _validated_redshift(z: np.ndarray | float) -> tuple[np.ndarray, bool]:
    scalar = np.ndim(z) == 0
    values = np.asarray(z, dtype=float)
    if np.any(~np.isfinite(values)):
        raise ValueError("Redshift must be finite")
    if np.any(values < 0.0):
        raise ValueError("Redshift must be non-negative")
    return values, scalar


def _power_integral(z: np.ndarray, exponent: float) -> np.ndarray:
    """Return [(1+z)^exponent-1]/exponent with its logarithmic limit."""

    log1pz = np.log1p(z)
    if abs(float(exponent)) <= 1e-10:
        return log1pz
    return np.expm1(float(exponent) * log1pz) / float(exponent)


def clock_path_integrand(
    z: np.ndarray | float,
    gamma_c: float,
    redshift_rate_power: float,
    minimum_light_speed_ratio: float | None = None,
    maximum_light_speed_ratio: float | None = None,
) -> np.ndarray | float:
    """Return [c(z)/c0] / [|dz/dT|/H0]."""

    values, scalar = _validated_redshift(z)
    gamma_c = float(gamma_c)
    redshift_rate_power = float(redshift_rate_power)
    if not np.isfinite(gamma_c) or not np.isfinite(redshift_rate_power):
        raise ValueError("Clock-law parameters must be finite")
    numerator = 1.0 + gamma_c * values
    minimum: float | None = None
    if minimum_light_speed_ratio is not None:
        minimum = float(minimum_light_speed_ratio)
        if not np.isfinite(minimum) or minimum <= 0.0:
            raise ValueError("The minimum light-speed ratio must be positive")
        numerator = np.maximum(numerator, minimum)
    if maximum_light_speed_ratio is not None:
        maximum = float(maximum_light_speed_ratio)
        if not np.isfinite(maximum) or maximum <= 0.0:
            raise ValueError("The maximum light-speed ratio must be positive")
        if minimum is not None and maximum < minimum:
            raise ValueError("The maximum light-speed ratio must not be below the minimum")
        numerator = np.minimum(numerator, maximum)
    if np.any(numerator <= 0.0):
        raise ValueError("The registered light-speed ratio must remain positive")
    result = numerator / np.power(1.0 + values, redshift_rate_power)
    return float(result) if scalar else result


def clock_path_integral(
    z: np.ndarray | float,
    gamma_c: float,
    redshift_rate_power: float,
) -> np.ndarray | float:
    """Return the dimensionless analytic integral of the clock-law kernel."""

    values, scalar = _validated_redshift(z)
    gamma_c = float(gamma_c)
    redshift_rate_power = float(redshift_rate_power)
    if not np.isfinite(gamma_c) or not np.isfinite(redshift_rate_power):
        raise ValueError("Clock-law parameters must be finite")
    if np.any(1.0 + gamma_c * values <= 0.0):
        raise ValueError("The registered light-speed ratio must remain positive")

    # 1 + gamma*z = gamma*(1+z) + (1-gamma).
    result = gamma_c * _power_integral(values, 2.0 - redshift_rate_power)
    result += (1.0 - gamma_c) * _power_integral(
        values, 1.0 - redshift_rate_power
    )
    return float(result) if scalar else result


def clock_path_distance(
    z: np.ndarray | float,
    h0: float,
    c0: float,
    law: ClockLaw,
) -> np.ndarray | float:
    """Return D_path=(c0/H0) times the registered dimensionless integral."""

    h0 = float(h0)
    c0 = float(c0)
    if not np.isfinite(h0) or h0 <= 0.0:
        raise ValueError("H0 must be finite and positive")
    if not np.isfinite(c0) or c0 <= 0.0:
        raise ValueError("c0 must be finite and positive")
    return (c0 / h0) * clock_path_integral(
        z, law.gamma_c, law.redshift_rate_power
    )
