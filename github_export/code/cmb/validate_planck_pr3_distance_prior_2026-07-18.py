#!/usr/bin/env python3
"""Recompute a registered Planck chain sample with CAMB and audit units."""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
from pathlib import Path

import camb
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
DATE = "2026-07-18"
RESULT_DIR = ROOT / "github_export" / "results" / DATE / "cmb"
PRIOR_PATH = RESULT_DIR / f"planck_pr3_distance_prior_{DATE}.json"
MODULE_PATH = Path(__file__).with_name("planck_distance_prior.py")
MAX_ABS_DIFFERENCE_SIGMA = 0.25


def load_module():
    spec = importlib.util.spec_from_file_location("registered_planck_distance_prior", MODULE_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prior", type=Path, default=PRIOR_PATH)
    parser.add_argument("--out-dir", type=Path, default=RESULT_DIR)
    args = parser.parse_args()
    payload = json.loads(args.prior.read_text(encoding="utf-8"))
    reference = payload["source"]["maximum_posterior_reference_sample"]

    parameters = camb.CAMBparams()
    parameters.set_cosmology(
        H0=reference["H0"],
        ombh2=reference["ombh2"],
        omch2=reference["omch2"],
        tau=reference["tau"],
        mnu=0.06,
        nnu=3.046,
        num_massive_neutrinos=1,
    )
    parameters.InitPower.set_params(As=math.exp(reference["logA"]) / 1.0e10, ns=reference["ns"])
    results = camb.get_background(parameters)
    module = load_module()
    camb_vector = module.vector_from_camb_results(results)
    chain_vector = np.asarray(reference["chain_vector"], dtype=float)
    sigma = np.sqrt(np.diag(np.asarray(payload["cov"], dtype=float)))
    difference_sigma = (camb_vector - chain_vector) / sigma
    passed = bool(np.max(np.abs(difference_sigma)) <= MAX_ABS_DIFFERENCE_SIGMA)

    prior = module.PlanckDistancePrior.load(args.prior)
    audit = {
        "analysis": "planck_pr3_distance_prior_camb_validation",
        "date": DATE,
        "reference_sample": reference,
        "chain_vector": chain_vector.tolist(),
        "camb_vector": camb_vector.tolist(),
        "difference": (camb_vector - chain_vector).tolist(),
        "difference_in_prior_sigma": difference_sigma.tolist(),
        "maximum_absolute_difference_in_prior_sigma": float(np.max(np.abs(difference_sigma))),
        "tolerance_in_prior_sigma": MAX_ABS_DIFFERENCE_SIGMA,
        "pass": passed,
        "reference_sample_chi_square_under_compression": prior.chi_square(camb_vector),
        "claim_boundary": "Software unit/recombination regression only; not an independent cosmological constraint.",
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.out_dir / f"planck_pr3_distance_prior_camb_validation_{DATE}.json"
    json_path.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    rows = [
        "# Planck Distance-prior CAMB Validation",
        "",
        "Date: 18 July 2026",
        "",
        "| parameter | chain | CAMB | difference / prior sigma |",
        "| --- | ---: | ---: | ---: |",
        *[
            f"| {label} | {chain:.10g} | {theory:.10g} | {pull:.8g} |"
            for label, chain, theory, pull in zip(
                payload["vector_order"], chain_vector, camb_vector, difference_sigma
            )
        ],
        "",
        f"Maximum absolute discrepancy: `{audit['maximum_absolute_difference_in_prior_sigma']:.8g}` prior sigma.",
        f"Registered software tolerance: `{MAX_ABS_DIFFERENCE_SIGMA}` prior sigma.",
        f"Gate: **{'PASS' if passed else 'FAIL'}**.",
        "",
        "This verifies the distance-vector units and current CAMB recombination calculation against one public Planck chain sample. It is not independent evidence.",
    ]
    report_path = args.out_dir / f"planck_pr3_distance_prior_camb_validation_report_{DATE}.md"
    report_path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(f"Maximum discrepancy: {audit['maximum_absolute_difference_in_prior_sigma']:.8g} sigma")
    print(f"Gate: {'PASS' if passed else 'FAIL'}")
    print(f"Saved: {report_path}")
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
