#!/usr/bin/env python3
"""Validate FIRAS ingestion and the primary-fit covariance safety gate."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from firas_monopole_likelihood import FirasMonopoleData, load_covariance, spectral_templates_kjy_sr


ROOT = Path(__file__).resolve().parents[3]
DATE = "2026-07-18"
DATA = (
    ROOT
    / "external_datasets"
    / "current_cosmology_sources"
    / "COBE_FIRAS"
    / "firas_monopole_spec_v1.txt"
)
OUT = ROOT / "github_export" / "results" / DATE / "cmb"


def main() -> None:
    data = FirasMonopoleData.load(DATA)
    templates = spectral_templates_kjy_sr(data.frequency_ghz)
    primary_gate_refused = False
    try:
        load_covariance(data, None)
    except RuntimeError:
        primary_gate_refused = True
    covariance, covariance_label = load_covariance(data, None, allow_diagonal_diagnostic=True)
    checks = {
        "channels": len(data.frequency_ghz) == 43,
        "strictly_increasing_frequency": bool(np.all(np.diff(data.frequency_ghz) > 0)),
        "finite_templates": all(bool(np.all(np.isfinite(values))) for values in templates.values()),
        "positive_diagonal_diagnostic": bool(np.all(np.linalg.eigvalsh(covariance) > 0)),
        "primary_fit_without_covariance_refused": primary_gate_refused,
    }
    result = {
        "analysis": "cobe_firas_ingestion_regression",
        "date": DATE,
        "data_path": str(DATA),
        "frequency_range_ghz": [float(data.frequency_ghz.min()), float(data.frequency_ghz.max())],
        "covariance_mode_exercised": covariance_label,
        "checks": checks,
        "pass": all(checks.values()),
        "scientific_fit_permitted": False,
        "claim_boundary": (
            "Ingestion and template software only. A primary spectral-distortion fit remains blocked "
            "until a full registered covariance and physical distortion prediction are supplied."
        ),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    json_path = OUT / f"firas_ingestion_regression_{DATE}.json"
    json_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    report = [
        "# COBE/FIRAS Ingestion Regression",
        "",
        "Date: 18 July 2026",
        "",
        f"- Channels: {len(data.frequency_ghz)}",
        f"- Frequency range: {data.frequency_ghz.min():.3f}-{data.frequency_ghz.max():.3f} GHz",
        f"- Software regression: {'PASS' if result['pass'] else 'FAIL'}",
        "- Scientific fit permitted: no",
        "",
        "The diagonal channel-error matrix was exercised only as a non-gating software diagnostic. "
        "The likelihood refused an unlabelled primary fit without a full covariance matrix.",
    ]
    report_path = OUT / f"firas_ingestion_regression_report_{DATE}.md"
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"Saved: {report_path}")


if __name__ == "__main__":
    main()
