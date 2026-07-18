#!/usr/bin/env python3
"""Run the preregistered ACT DR6-lite collaboration regression test."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
DATE = "2026-07-18"
ACT_ROOT = (
    ROOT
    / "external_datasets"
    / "current_cosmology_sources"
    / "ACT_DR6_lite"
)
ACT_REPO = ACT_ROOT / "DR6-ACT-lite"
PACKAGES = ACT_ROOT / "cobaya_packages"
DATA = PACKAGES / "data" / "ACTDR6CMBonly" / "v1.0" / "dr6_data_cmbonly.fits"
RESULT_DIR = ROOT / "github_export" / "results" / DATE / "cmb"
PREREG = RESULT_DIR / f"cmb_programme_preregistration_{DATE}.json"
EXPECTED = -395.48
TOLERANCE = 0.01


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_regression() -> dict[str, object]:
    if not DATA.exists():
        raise FileNotFoundError(DATA)
    sys.path.insert(0, str(ACT_REPO))
    from cobaya.model import get_model

    info = {
        "packages_path": str(PACKAGES),
        "params": {
            "ombh2": 0.022,
            "omch2": 0.117,
            "ns": 0.96,
            "As": 2e-9,
            "tau": 0.065,
            "cosmomc_theta": 104.09e-4,
        },
        "theory": {
            "camb": {
                "extra_args": {
                    "lmax": 9000,
                    "lens_potential_accuracy": 8,
                    "min_l_logl_sampling": 6000,
                }
            }
        },
        "likelihood": {
            "act_dr6_cmbonly": {
                "stop_at_error": True,
                "input_file": "dr6_data_cmbonly.fits",
                "params": {"A_act": 1.0, "P_act": 1.0},
            }
        },
        "sampler": {"evaluate": None},
        "debug": False,
    }
    model = get_model(info)
    try:
        values, derived = model.loglikes()
        names = list(model.likelihood)
        total = float(np.sum(values))
    finally:
        model.close()
    difference = total - EXPECTED
    return {
        "analysis": "act_dr6_lite_regression",
        "date": DATE,
        "likelihood_names": names,
        "individual_loglikes": [float(value) for value in values],
        "derived_values": [float(value) for value in derived if np.isscalar(value)],
        "total_loglike": total,
        "expected_loglike": EXPECTED,
        "absolute_tolerance": TOLERANCE,
        "difference": difference,
        "pass": bool(abs(difference) <= TOLERANCE),
        "data": {
            "path": str(DATA),
            "bytes": DATA.stat().st_size,
            "sha256": sha256_file(DATA),
        },
        "preregistration_sha256": sha256_file(PREREG),
        "claim_boundary": "Software and data regression only; not a cosmological fit.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=RESULT_DIR)
    args = parser.parse_args()
    if not PREREG.exists():
        raise FileNotFoundError(PREREG)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    result = run_regression()
    json_path = args.out_dir / f"act_dr6_lite_regression_{DATE}.json"
    json_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    report = [
        "# ACT DR6-lite Regression",
        "",
        "Date: 18 July 2026",
        "",
        f"- measured log-likelihood: `{result['total_loglike']:.8f}`",
        f"- collaboration target: `{result['expected_loglike']:.8f}`",
        f"- difference: `{result['difference']:.8g}`",
        f"- tolerance: `{result['absolute_tolerance']}`",
        f"- gate: **{'PASS' if result['pass'] else 'FAIL'}**",
        "",
        "This is a software/data regression only. It verifies that the local CAMB, Cobaya, SACC data and ACT likelihood reproduce the collaboration test; it is not a cosmological parameter fit or evidence comparison.",
    ]
    report_path = args.out_dir / f"act_dr6_lite_regression_report_{DATE}.md"
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"ACT DR6-lite loglike: {result['total_loglike']:.8f}")
    print(f"Gate: {'PASS' if result['pass'] else 'FAIL'}")
    print(f"Saved: {report_path}")
    if not result["pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
