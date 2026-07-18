#!/usr/bin/env python3
"""Audit the reconciled p=0 and p=1 PLAMB clock-law implementations."""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad


ROOT = Path(__file__).resolve().parents[3]
CODE_DIR = Path(__file__).resolve().parent
SHARED_DIR = ROOT / "github_export" / "code" / "shared"
for path in (CODE_DIR, SHARED_DIR, ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import run_rawmu_release_grounded_holdouts_2026_07_18 as rawmu  # noqa: E402
from plamb_clock_distance import (  # noqa: E402
    CONSTANT_C_LINEAR_REDSHIFT,
    FRACTIONAL_CLOCK,
    PETER_LINEAR_REDSHIFT,
    clock_path_distance,
    clock_path_integral,
    clock_path_integrand,
)


DATE = "2026-07-18"
OUT = ROOT / "github_export" / "results" / DATE / "rawmu_plamb_axis"
SUMMARY = OUT / f"plamb_axis_comparison_summary_{DATE}.json"
EXPECTED_DELTA_BIC = 94.34498313933864
TOLERANCE = 2e-12


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_pilot_module():
    path = SHARED_DIR / "PLamb_Test_10V6c_plus.py"
    spec = importlib.util.spec_from_file_location("plamb_clock_audit_pilot", path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def quadrature_maximum() -> float:
    maximum = 0.0
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
                maximum = max(maximum, abs(measured - expected))
    return maximum


def write_manifest(paths: list[Path]) -> Path:
    manifest = OUT / f"plamb_clock_law_consistency_manifest_{DATE}.csv"
    rows = [
        {
            "path": path.relative_to(ROOT).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in sorted(set(paths), key=lambda item: item.as_posix())
    ]
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["path", "bytes", "sha256"])
        writer.writeheader()
        writer.writerows(rows)
    return manifest


def main() -> None:
    grid = np.linspace(0.0, 2.4, 241)
    peter_expected = grid * (1.0 + 0.5 * grid)
    fractional_expected = grid
    constant_c_expected = grid
    peter_difference = float(
        np.max(
            np.abs(
                clock_path_integral(
                    grid,
                    PETER_LINEAR_REDSHIFT.gamma_c,
                    PETER_LINEAR_REDSHIFT.redshift_rate_power,
                )
                - peter_expected
            )
        )
    )
    fractional_difference = float(
        np.max(
            np.abs(
                clock_path_integral(
                    grid,
                    FRACTIONAL_CLOCK.gamma_c,
                    FRACTIONAL_CLOCK.redshift_rate_power,
                )
                - fractional_expected
            )
        )
    )
    constant_c_difference = float(
        np.max(
            np.abs(
                clock_path_integral(
                    grid,
                    CONSTANT_C_LINEAR_REDSHIFT.gamma_c,
                    CONSTANT_C_LINEAR_REDSHIFT.redshift_rate_power,
                )
                - constant_c_expected
            )
        )
    )
    quadrature_difference = quadrature_maximum()

    positive_grid = grid[1:]
    rawmu_measured = rawmu.model_mu(
        positive_grid, "FR_C1PZ_ALPHA0", None, 67.5
    )
    rawmu_expected = (
        5.0
        * np.log10(
            clock_path_distance(
                positive_grid, 67.5, rawmu.C_KMS, PETER_LINEAR_REDSHIFT
            )
        )
        + 25.0
    )
    rawmu_wiring_difference = float(np.max(np.abs(rawmu_measured - rawmu_expected)))

    pilot = load_pilot_module()
    _dm, _dh, _dv, pilot_distance = pilot.plamb_noexp_distances(
        positive_grid,
        67.5,
        A_clump=0.0,
        n_clump=1.0,
        c_drift=1.0,
        dimming_alpha=0.0,
    )
    pilot_expected = (pilot.C_KMS / 67.5) * positive_grid
    pilot_wiring_difference = float(
        np.max(np.abs(np.asarray(pilot_distance) - pilot_expected))
    )

    comparison = json.loads(SUMMARY.read_text(encoding="utf-8"))
    measured_delta_bic = float(
        comparison["primary_fit"]["delta_BIC_FR_minus_LCDM"]
    )
    delta_bic_regression_difference = abs(measured_delta_bic - EXPECTED_DELTA_BIC)

    checks = {
        "peter_p0_identity": peter_difference <= TOLERANCE,
        "fractional_p1_identity": fractional_difference <= TOLERANCE,
        "constant_c_identity": constant_c_difference <= TOLERANCE,
        "analytic_vs_quadrature": quadrature_difference <= TOLERANCE,
        "rawmu_p0_wiring": rawmu_wiring_difference <= TOLERANCE,
        "historical_pilot_p1_wiring": pilot_wiring_difference <= 2e-9,
        "fixed_likelihood_regression": delta_bic_regression_difference <= 1e-9,
    }
    result = {
        "date": DATE,
        "decision": "PASS_CLOCK_LAW_RECONCILIATION" if all(checks.values()) else "FAIL",
        "checks": checks,
        "maximum_absolute_differences": {
            "peter_p0_identity": peter_difference,
            "fractional_p1_identity": fractional_difference,
            "constant_c_identity": constant_c_difference,
            "analytic_vs_quadrature": quadrature_difference,
            "rawmu_p0_wiring_mag": rawmu_wiring_difference,
            "historical_pilot_p1_wiring_mpc": pilot_wiring_difference,
            "delta_bic_regression": delta_bic_regression_difference,
        },
        "registered_branches": {
            "rawmu_fixed_plamb": {
                "clock_law": PETER_LINEAR_REDSHIFT.key,
                "gamma_c": PETER_LINEAR_REDSHIFT.gamma_c,
                "p": PETER_LINEAR_REDSHIFT.redshift_rate_power,
                "alpha": 0.0,
            },
            "historical_noexp_pilot": {
                "clock_law": "fractional_clock_pilot",
                "gamma_c": "c_drift parameter",
                "p": float(pilot.PLAMB_PILOT_REDSHIFT_RATE_POWER),
                "alpha": "plamb_alpha parameter; historical default 0.5",
            },
        },
        "fixed_plamb_regression": {
            "N": int(comparison["primary_fit"]["N"]),
            "delta_BIC_PLAMB_minus_LCDM": measured_delta_bic,
            "expected_before_refactor": EXPECTED_DELTA_BIC,
        },
        "claim_boundary": (
            "The audit proves analytic and implementation consistency only. It does not "
            "select p=0 or p=1 physically; that requires the nested likelihood and external "
            "clock/flux constraints."
        ),
    }
    if not all(checks.values()):
        failed = [name for name, passed in checks.items() if not passed]
        raise RuntimeError(f"Clock-law audit failed: {failed}")

    json_path = OUT / f"plamb_clock_law_consistency_audit_{DATE}.json"
    json_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    report = f"""# PLAMB Clock-law Consistency Audit

Date: 18 July 2026

## Decision

**PASS CLOCK-LAW RECONCILIATION.**

The raw-MU and historical no-expansion branches now use one shared clock-law
kernel while retaining distinct, explicit physical assumptions.

## Registered branches

| branch | gamma_c | p | consequence at gamma_c=1 |
| --- | ---: | ---: | --- |
| fixed raw-MU PLAMB | 1 | 0 | `D_path=(c0/H0) z(1+z/2)` |
| historical no-expansion pilot | parameter | 1 | `D_path=(c0/H0) z` when unclumped |
| constant-c linear-redshift control | 0 | 0 | `D_path=(c0/H0) z` |

The branches are no longer mathematically conflated. The raw-MU model uses a
linear redshift rate, whereas the earlier pilot uses a fractional redshift-rate
factor.

## Numerical gates

| gate | maximum absolute difference | status |
| --- | ---: | --- |
| Peter p=0 identity | `{peter_difference:.3e}` | PASS |
| fractional-clock p=1 identity | `{fractional_difference:.3e}` | PASS |
| constant-c control | `{constant_c_difference:.3e}` | PASS |
| analytic expression versus quadrature | `{quadrature_difference:.3e}` | PASS |
| raw-MU p=0 wiring [mag] | `{rawmu_wiring_difference:.3e}` | PASS |
| historical pilot p=1 wiring [Mpc] | `{pilot_wiring_difference:.3e}` | PASS |

## Likelihood regression

Re-running the three-release fixed-law comparison gives
`Delta BIC(PLAMB-Lambda-CDM)={measured_delta_bic:.12f}` for `N={int(comparison['primary_fit']['N'])}`.
The difference from the pre-refactor value is
`{delta_bic_regression_difference:.3e}`. The refactor therefore changes no
scientific result.

## Claim boundary

This audit proves analytic and implementation consistency. It does not decide
whether `p=0` or `p=1` is physically correct. That question belongs to the
separately preregistered nested clock-law likelihood with external time-dilation
and flux constraints.

## Reproduction

```powershell
python github_export/code/shared/test_plamb_clock_distance.py
python github_export/code/rawmu/run_rawmu_release_grounded_holdouts_2026_07_18.py --self-test-only
python github_export/code/rawmu/run_plamb_axis_comparison_2026_07_18.py
python github_export/code/rawmu/audit_plamb_clock_law_consistency_2026_07_18.py
```
"""
    report_path = OUT / f"plamb_clock_law_consistency_report_{DATE}.md"
    report_path.write_text(report, encoding="utf-8")

    source_paths = [
        Path(__file__).resolve(),
        SHARED_DIR / "plamb_clock_distance.py",
        SHARED_DIR / "test_plamb_clock_distance.py",
        SHARED_DIR / "PLamb_Test_10V6c_plus.py",
        CODE_DIR / "run_rawmu_release_grounded_holdouts_2026_07_18.py",
        CODE_DIR / "run_plamb_axis_comparison_2026_07_18.py",
        CODE_DIR / "analyze_fr_c1pz_noloss_lcdm_calibration_ladder_2026-07-16.py",
        CODE_DIR / "analyze_fr_c1pz_noloss_dimming_sequence_2026-07-16.py",
        json_path,
        report_path,
        SUMMARY,
    ]
    manifest = write_manifest(source_paths)
    print(f"Decision: {result['decision']}")
    print(f"Saved: {report_path}")
    print(f"Saved: {manifest}")


if __name__ == "__main__":
    main()
