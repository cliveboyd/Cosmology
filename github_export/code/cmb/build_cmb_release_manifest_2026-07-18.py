#!/usr/bin/env python3
"""Build the dated CMB machine summary and SHA-256 manifest."""

from __future__ import annotations

import csv
import hashlib
import json
import platform
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
DATE = "2026-07-18"
CODE_DIR = ROOT / "github_export" / "code" / "cmb"
RESULT_DIR = ROOT / "github_export" / "results" / DATE / "cmb"
MANIFEST = RESULT_DIR / f"cmb_release_manifest_{DATE}.csv"
SUMMARY = RESULT_DIR / f"cmb_implementation_summary_{DATE}.json"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(name: str) -> dict[str, object]:
    return json.loads((RESULT_DIR / name).read_text(encoding="utf-8"))


def package_version(name: str) -> str | None:
    try:
        return version(name)
    except PackageNotFoundError:
        return None


def main() -> None:
    planck = load_json(f"planck_pr3_distance_prior_{DATE}.json")
    planck_validation = load_json(f"planck_pr3_distance_prior_camb_validation_{DATE}.json")
    act_lite = load_json(f"act_dr6_lite_regression_{DATE}.json")
    firas = load_json(f"firas_ingestion_regression_{DATE}.json")
    readiness_path = RESULT_DIR / f"cmb_branch_readiness_{DATE}.json"
    act_lensing_path = RESULT_DIR / f"act_dr6_lensing_regression_{DATE}.json"
    su2r_readiness_path = RESULT_DIR / f"su2r_physical_readiness_audit_{DATE}.json"
    summary = {
        "date": DATE,
        "software_environment": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "numpy": package_version("numpy"),
            "scipy": package_version("scipy"),
            "camb": package_version("camb"),
            "cobaya": package_version("cobaya"),
            "sacc": package_version("sacc"),
            "act_dr6_lenslike": package_version("act_dr6_lenslike"),
        },
        "planck_distance_prior": {
            "posterior_rows": planck["posterior_rows"],
            "means": planck["means"],
            "camb_validation_pass": planck_validation["pass"],
            "maximum_camb_discrepancy_sigma": planck_validation[
                "maximum_absolute_difference_in_prior_sigma"
            ],
        },
        "act_dr6_lite": {
            "loglike": act_lite["total_loglike"],
            "target": act_lite["expected_loglike"],
            "pass": act_lite["pass"],
        },
        "act_dr6_lensing": (
            {
                "pass": load_json(act_lensing_path.name)["pass"],
                "cells": load_json(act_lensing_path.name)["cells"],
            }
            if act_lensing_path.exists()
            else {"pass": False, "status": "official data regression pending"}
        ),
        "firas": {
            "ingestion_pass": firas["pass"],
            "scientific_fit_permitted": firas["scientific_fit_permitted"],
            "frequency_range_ghz": firas["frequency_range_ghz"],
        },
        "readiness": (
            load_json(readiness_path.name)["branches"] if readiness_path.exists() else []
        ),
        "su2r_physical_adapter": (
            {
                "physical_primary_cmb_ready": load_json(su2r_readiness_path.name)[
                    "physical_primary_cmb_ready"
                ],
                "physical_lensing_ready": load_json(su2r_readiness_path.name)[
                    "physical_lensing_ready"
                ],
                "decision": load_json(su2r_readiness_path.name)["decision"],
            }
            if su2r_readiness_path.exists()
            else {"status": "physical-readiness audit pending"}
        ),
        "claim_boundary": (
            "These are software, data and background-compression foundations. Physical SU(2), "
            "SU2R, magnetic, parity or FR interpretation requires the registered missing theory "
            "and data capabilities."
        ),
    }
    SUMMARY.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    paths = sorted(
        [path for path in CODE_DIR.rglob("*") if path.is_file() and "__pycache__" not in path.parts]
        + [
            path
            for path in RESULT_DIR.rglob("*")
            if path.is_file() and path != MANIFEST
        ],
        key=lambda path: path.relative_to(ROOT).as_posix(),
    )
    rows = [
        {
            "path": path.relative_to(ROOT).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in paths
    ]
    with MANIFEST.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["path", "bytes", "sha256"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved: {SUMMARY}")
    print(f"Saved: {MANIFEST}")
    print(f"Manifest files: {len(rows)}")


if __name__ == "__main__":
    main()
