#!/usr/bin/env python3

from __future__ import annotations

import unittest
from pathlib import Path

import numpy as np

from cmb_spectrum_transforms import phi_phi_to_kappa_kappa, rotate_even_parity_spectra, tensor_chirality
from cmb_theory_contract import BackgroundPrediction, planck_distance_vector
from firas_monopole_likelihood import FirasMonopoleData, load_covariance, spectral_templates_kjy_sr
from low_l_global_null import calibrate_global_null


ROOT = Path(__file__).resolve().parents[3]
FIRAS_PATH = (
    ROOT
    / "external_datasets"
    / "current_cosmology_sources"
    / "COBE_FIRAS"
    / "firas_monopole_spec_v1.txt"
)


class CmbContractTests(unittest.TestCase):
    def test_phi_to_kappa(self) -> None:
        ell = np.arange(5, dtype=float)
        phi = np.ones(5)
        expected = np.square(ell * (ell + 1.0)) / 4.0
        np.testing.assert_allclose(phi_phi_to_kappa_kappa(ell, phi), expected)

    def test_zero_birefringence(self) -> None:
        te = np.array([0.0, 1.0, 2.0])
        ee = np.array([0.0, 3.0, 4.0])
        bb = np.array([0.0, 0.2, 0.3])
        rotated = rotate_even_parity_spectra(0.0, te, ee, bb)
        np.testing.assert_allclose(rotated["TE"], te)
        np.testing.assert_allclose(rotated["EE"], ee)
        np.testing.assert_allclose(rotated["BB"], bb)
        np.testing.assert_allclose(rotated["TB"], 0.0)
        np.testing.assert_allclose(rotated["EB"], 0.0)

    def test_chirality_bounds(self) -> None:
        value = tensor_chirality(np.array([1.0, 0.0, 2.0]), np.array([0.0, 1.0, 2.0]))
        np.testing.assert_allclose(value, [1.0, -1.0, 0.0])

    def test_firas_templates_are_finite(self) -> None:
        templates = spectral_templates_kjy_sr(np.linspace(70.0, 640.0, 43))
        for values in templates.values():
            self.assertTrue(np.all(np.isfinite(values)))

    @unittest.skipUnless(FIRAS_PATH.exists(), "FIRAS table not downloaded")
    def test_firas_table_and_covariance_gate(self) -> None:
        data = FirasMonopoleData.load(FIRAS_PATH)
        self.assertEqual(data.frequency_ghz.shape, (43,))
        self.assertAlmostEqual(float(data.frequency_ghz[0]), 68.0529, places=3)
        self.assertAlmostEqual(float(data.frequency_ghz[-1]), 639.4573, places=3)
        with self.assertRaises(RuntimeError):
            load_covariance(data, None)
        covariance, label = load_covariance(data, None, allow_diagonal_diagnostic=True)
        self.assertEqual(covariance.shape, (43, 43))
        self.assertEqual(label, "diagonal_diagnostic_only")

    def test_distance_prior_requires_recombination_bracket(self) -> None:
        prediction = BackgroundPrediction(
            z=np.array([0.0, 1.0]),
            hubble_km_s_mpc=np.array([70.0, 120.0]),
            comoving_distance_mpc=np.array([0.0, 3300.0]),
            z_star=1090.0,
            sound_horizon_star_mpc=144.0,
            omega_b_h2=0.0224,
            omega_m=0.3,
        )
        with self.assertRaises(ValueError):
            planck_distance_vector(prediction)

    def test_low_l_global_null_is_rank_calibrated(self) -> None:
        simulations = np.column_stack(
            [np.linspace(-2.0, 2.0, 101), np.linspace(3.0, -3.0, 101)]
        )
        result = calibrate_global_null(
            np.array([0.0, 0.0]), simulations, ["two-sided", "greater"]
        )
        self.assertTrue(0.0 < result.global_p <= 1.0)
        self.assertEqual(result.observed_local_p.shape, (2,))
        self.assertEqual(result.simulation_max_scores.shape, (101,))


if __name__ == "__main__":
    unittest.main()
