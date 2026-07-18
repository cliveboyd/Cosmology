#!/usr/bin/env python3
"""Audit which registered CMB branches are executable without silent proxies."""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
DATE = "2026-07-18"
RESULT_DIR = ROOT / "github_export" / "results" / DATE / "cmb"
EXTERNAL = ROOT / "external_datasets" / "current_cosmology_sources"


def json_pass(path: Path) -> bool:
    if not path.exists():
        return False
    return bool(json.loads(path.read_text(encoding="utf-8")).get("pass"))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=RESULT_DIR)
    args = parser.parse_args()
    planck_prior = RESULT_DIR / f"planck_pr3_distance_prior_{DATE}.json"
    planck_validation = RESULT_DIR / f"planck_pr3_distance_prior_camb_validation_{DATE}.json"
    act_lite = RESULT_DIR / f"act_dr6_lite_regression_{DATE}.json"
    act_lensing = RESULT_DIR / f"act_dr6_lensing_regression_{DATE}.json"
    firas = EXTERNAL / "COBE_FIRAS" / "firas_monopole_spec_v1.txt"
    firas_covariance = EXTERNAL / "COBE_FIRAS" / "firas_monopole_covariance.npy"
    pr4_maps = EXTERNAL / "Planck_PR4_NPIPE" / "maps"
    pr4_simulations = EXTERNAL / "Planck_PR4_NPIPE" / "simulations"
    angle_priors = EXTERNAL / "CMB_polarisation_angle_priors" / "angle_priors.json"
    act_lite_data = (
        EXTERNAL
        / "ACT_DR6_lite"
        / "cobaya_packages"
        / "data"
        / "ACTDR6CMBonly"
        / "v1.0"
        / "dr6_data_cmbonly.fits"
    )

    rows = [
        {
            "branch": "planck_distance_prior",
            "data_ready": planck_prior.exists(),
            "software_ready": json_pass(planck_validation),
            "theory_ready": True,
            "scientific_ready": planck_prior.exists() and json_pass(planck_validation),
            "boundary": "standard early physics and recombination only",
        },
        {
            "branch": "act_dr6_primary_cmb",
            "data_ready": act_lite_data.exists(),
            "software_ready": json_pass(act_lite),
            "theory_ready": importlib.util.find_spec("camb") is not None,
            "scientific_ready": act_lite_data.exists() and json_pass(act_lite),
            "boundary": "new-model use requires CAMB/CLASS-compatible spectra",
        },
        {
            "branch": "firas_spectral_distortion",
            "data_ready": firas.exists() and firas_covariance.exists(),
            "software_ready": firas.exists(),
            "theory_ready": False,
            "scientific_ready": False,
            "boundary": "blocked by full covariance and registered dQ/dz or photon-kinetic mapping",
        },
        {
            "branch": "planck_low_l",
            "data_ready": pr4_maps.exists() and pr4_simulations.exists(),
            "software_ready": importlib.util.find_spec("healpy") is not None,
            "theory_ready": True,
            "scientific_ready": False,
            "boundary": "blocked by PR4 maps/masks/simulations and frozen global statistic family",
        },
        {
            "branch": "tb_eb_parity",
            "data_ready": pr4_maps.exists() and angle_priors.exists(),
            "software_ready": importlib.util.find_spec("pymaster") is not None,
            "theory_ready": False,
            "scientific_ready": False,
            "boundary": "blocked by split maps, angle priors, foreground model and parity-odd transfer spectra",
        },
        {
            "branch": "act_dr6_lensing",
            "data_ready": (EXTERNAL / "ACT_DR6_lensing" / "data" / "v1.2").exists(),
            "software_ready": json_pass(act_lensing),
            "theory_ready": False,
            "scientific_ready": False,
            "boundary": "Lambda-CDM regression ready; physical SU2R inference is blocked by perturbation equations",
        },
    ]
    result = {"date": DATE, "branches": rows}
    args.out_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.out_dir / f"cmb_branch_readiness_{DATE}.json"
    json_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    report = [
        "# CMB Branch Readiness",
        "",
        "Date: 18 July 2026",
        "",
        "| branch | data | software | theory | scientific | boundary |",
        "| --- | --- | --- | --- | --- | --- |",
        *[
            f"| {row['branch']} | {row['data_ready']} | {row['software_ready']} | {row['theory_ready']} | {row['scientific_ready']} | {row['boundary']} |"
            for row in rows
        ],
        "",
        "A false scientific-readiness flag is a hard implementation block. It must not be bypassed with a phenomenological proxy carrying a physical model label.",
    ]
    report_path = args.out_dir / f"cmb_branch_readiness_report_{DATE}.md"
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"Saved: {report_path}")


if __name__ == "__main__":
    main()
