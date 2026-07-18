#!/usr/bin/env python3
"""Run the registered ACT DR6 lensing v1.2 likelihood regressions."""

from __future__ import annotations

import argparse
import hashlib
import json
from importlib.metadata import version
from pathlib import Path

import act_dr6_lenslike as alike
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
DATE = "2026-07-18"
RESULT_DIR = ROOT / "github_export" / "results" / DATE / "cmb"
REGISTRATION = RESULT_DIR / f"act_dr6_lensing_regression_registration_{DATE}.json"
DATA_DIR = (
    ROOT
    / "external_datasets"
    / "current_cosmology_sources"
    / "ACT_DR6_lensing"
    / "data"
    / "v1.2"
)
ARCHIVE = DATA_DIR.parents[1] / "ACT_dr6_likelihood_v1.2.tgz"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fiducial_spectra(data_dir: Path) -> tuple[np.ndarray, ...]:
    ell, cl_tt, cl_ee, cl_bb, cl_te = np.loadtxt(
        data_dir / "like_corrs" / "cosmo2017_10K_acc3_lensedCls.dat", unpack=True
    )
    ell_p, cl_pp = np.loadtxt(
        data_dir / "like_corrs" / "cosmo2017_10K_acc3_lenspotentialCls.dat",
        unpack=True,
        usecols=[0, 5],
    )
    prefactor = 2.0 * np.pi / ell / (ell + 1.0)
    cl_kk = cl_pp * 2.0 * np.pi / 4.0
    return ell_p, cl_kk, ell, cl_tt * prefactor, cl_ee * prefactor, cl_te * prefactor, cl_bb * prefactor


def evaluate_cell(cell: dict[str, object], tolerance: float, spectra: tuple[np.ndarray, ...]) -> dict[str, object]:
    ell_kk, cl_kk, ell, cl_tt, cl_ee, cl_te, cl_bb = spectra
    lens_only = bool(cell["lens_only"])
    data = alike.load_data(
        str(cell["variant"]),
        ddir=str(DATA_DIR),
        lens_only=lens_only,
        like_corrections=not lens_only,
        version="v1.2",
    )
    loglike = alike.generic_lnlike(data, ell_kk, cl_kk, ell, cl_tt, cl_ee, cl_te, cl_bb)
    chi_square = float(-2.0 * loglike)
    expected = float(cell["expected_chi_square"])
    difference = chi_square - expected
    return {
        **cell,
        "measured_chi_square": chi_square,
        "difference": difference,
        "pass": bool(abs(difference) <= tolerance),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=RESULT_DIR)
    args = parser.parse_args()
    registration = json.loads(REGISTRATION.read_text(encoding="utf-8"))
    if not DATA_DIR.exists():
        raise FileNotFoundError(DATA_DIR)
    spectra = fiducial_spectra(DATA_DIR)
    tolerance = float(registration["absolute_chi_square_tolerance"])
    rows = [evaluate_cell(cell, tolerance, spectra) for cell in registration["cells"]]
    passed = all(bool(row["pass"]) for row in rows)
    result = {
        "analysis": "act_dr6_lensing_v1p2_regression",
        "date": DATE,
        "package_version": version("act_dr6_lenslike"),
        "data_archive_sha256": sha256_file(ARCHIVE),
        "registration_sha256": sha256_file(REGISTRATION),
        "cells": rows,
        "pass": passed,
        "claim_boundary": registration["claim_boundary"],
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.out_dir / f"act_dr6_lensing_regression_{DATE}.json"
    json_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    report = [
        "# ACT DR6 Lensing Regression",
        "",
        "Date: 18 July 2026",
        "",
        "| variant | lens only | expected chi-squared | measured chi-squared | difference | pass |",
        "| --- | --- | ---: | ---: | ---: | --- |",
        *[
            f"| {row['variant']} | {row['lens_only']} | {row['expected_chi_square']:.5g} | {row['measured_chi_square']:.8g} | {row['difference']:.8g} | {row['pass']} |"
            for row in rows
        ],
        "",
        f"Overall software/data gate: **{'PASS' if passed else 'FAIL'}**.",
        "",
        "This is a package regression, not a cosmological fit. SU2R inference remains blocked until registered perturbation equations predict convergence and lensed primary-CMB spectra.",
    ]
    report_path = args.out_dir / f"act_dr6_lensing_regression_report_{DATE}.md"
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"Overall gate: {'PASS' if passed else 'FAIL'}")
    for row in rows:
        print(f"{row['variant']} lens_only={row['lens_only']}: chi2={row['measured_chi_square']:.8g}")
    print(f"Saved: {report_path}")
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
