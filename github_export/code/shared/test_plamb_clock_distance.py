#!/usr/bin/env python3
"""Contract tests for the shared PLAMB clock-law path distance."""

from __future__ import annotations

import unittest

import numpy as np
from scipy.integrate import quad

from plamb_clock_distance import (
    CLOCK_LAWS,
    CONSTANT_C_LINEAR_REDSHIFT,
    FRACTIONAL_CLOCK,
    PETER_LINEAR_REDSHIFT,
    clock_path_distance,
    clock_path_integral,
    clock_path_integrand,
)


class ClockPathTests(unittest.TestCase):
    def setUp(self) -> None:
        self.z = np.array([0.0, 0.01, 0.2, 1.0, 2.4])

    def test_peter_linear_redshift_identity(self) -> None:
        measured = clock_path_integral(self.z, gamma_c=1.0, redshift_rate_power=0.0)
        expected = self.z * (1.0 + 0.5 * self.z)
        np.testing.assert_allclose(measured, expected, rtol=2e-15, atol=2e-15)

    def test_fractional_clock_identity(self) -> None:
        measured = clock_path_integral(self.z, gamma_c=1.0, redshift_rate_power=1.0)
        np.testing.assert_allclose(measured, self.z, rtol=2e-15, atol=2e-15)

    def test_constant_c_linear_redshift_identity(self) -> None:
        measured = clock_path_integral(self.z, gamma_c=0.0, redshift_rate_power=0.0)
        np.testing.assert_allclose(measured, self.z, rtol=2e-15, atol=2e-15)

    def test_analytic_integral_matches_quadrature(self) -> None:
        maximum_difference = 0.0
        for gamma_c in (0.0, 0.25, 1.0):
            for power in (-0.5, 0.0, 0.7, 1.0, 1.3, 2.0, 2.4):
                for redshift in (0.01, 0.4, 1.0, 2.4):
                    expected = quad(
                        lambda value: float(
                            clock_path_integrand(value, gamma_c, power)
                        ),
                        0.0,
                        redshift,
                        epsabs=2e-13,
                        epsrel=2e-13,
                    )[0]
                    measured = float(clock_path_integral(redshift, gamma_c, power))
                    maximum_difference = max(maximum_difference, abs(measured - expected))
        self.assertLess(maximum_difference, 2e-12)

    def test_named_laws_are_distinct_and_stable(self) -> None:
        self.assertEqual(set(CLOCK_LAWS), {
            "peter_linear_redshift",
            "fractional_clock",
            "constant_c_linear_redshift",
        })
        self.assertEqual(PETER_LINEAR_REDSHIFT.redshift_rate_power, 0.0)
        self.assertEqual(FRACTIONAL_CLOCK.redshift_rate_power, 1.0)
        self.assertEqual(CONSTANT_C_LINEAR_REDSHIFT.gamma_c, 0.0)

    def test_distance_scale_and_scalar_contract(self) -> None:
        value = clock_path_distance(1.0, 67.5, 299792.458, PETER_LINEAR_REDSHIFT)
        self.assertIsInstance(value, float)
        self.assertAlmostEqual(value, (299792.458 / 67.5) * 1.5, places=11)

    def test_invalid_inputs_are_refused(self) -> None:
        with self.assertRaises(ValueError):
            clock_path_integral(-0.1, 1.0, 0.0)
        with self.assertRaises(ValueError):
            clock_path_integrand(2.0, -1.0, 0.0)
        with self.assertRaises(ValueError):
            clock_path_distance(1.0, 0.0, 299792.458, PETER_LINEAR_REDSHIFT)

    def test_historical_pilot_clipping_is_explicit(self) -> None:
        measured = clock_path_integrand(
            np.array([0.0, 2.0]),
            gamma_c=-1.0,
            redshift_rate_power=1.0,
            minimum_light_speed_ratio=0.05,
            maximum_light_speed_ratio=20.0,
        )
        np.testing.assert_allclose(measured, np.array([1.0, 0.05 / 3.0]))
        with self.assertRaises(ValueError):
            clock_path_integrand(
                1.0,
                gamma_c=0.0,
                redshift_rate_power=1.0,
                minimum_light_speed_ratio=2.0,
                maximum_light_speed_ratio=1.0,
            )


if __name__ == "__main__":
    unittest.main()
