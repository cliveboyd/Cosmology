#!/usr/bin/env python3
"""Audit the existing SU2R closures against the physical CMB adapter contract."""

from __future__ import annotations

import json
from pathlib import Path

from su2r_physical_perturbation_adapter import PhysicalTheoryIncomplete, SU2RRegistry, sha256_file


ROOT = Path(__file__).resolve().parents[3]
DATE = "2026-07-18"
RESULT_DIR = ROOT / "github_export" / "results" / DATE / "cmb"
REGISTRY = RESULT_DIR / f"su2r_physical_equation_registry_{DATE}.json"
GROWTH_CLOSURE = ROOT / "diagnose_su2r_perturbations.py"
PPF_PROXY = ROOT / "run_act_dr6_lite_su2r_ppf_proxy.py"
GROWTH_REPORT = (
    ROOT
    / "plamb_runs"
    / "diagnostics"
    / "su2r_growth_likelihood"
    / "nested_closure"
    / "su2r_growth_nested_closure_summary.md"
)
PPF_REPORT = (
    ROOT
    / "plamb_runs"
    / "diagnostics"
    / "act_dr6_lite_su2r_ppf_proxy"
    / "act_dr6_lite_su2r_ppf_proxy.md"
)


def file_record(path: Path) -> dict[str, object]:
    return {
        "path": str(path),
        "exists": path.exists(),
        "bytes": path.stat().st_size if path.exists() else None,
        "sha256": sha256_file(path) if path.exists() else None,
    }


def main() -> None:
    registry = SU2RRegistry.load(REGISTRY)
    missing_lensing = registry.missing("lensing")
    missing_primary = registry.missing("primary_cmb")
    refusal_verified = False
    try:
        registry.require("lensing")
    except PhysicalTheoryIncomplete:
        refusal_verified = True
    result = {
        "date": DATE,
        "physical_lensing_ready": not missing_lensing,
        "physical_primary_cmb_ready": not missing_primary,
        "incomplete_registry_refusal_verified": refusal_verified,
        "missing_lensing_requirements": missing_lensing,
        "missing_primary_cmb_requirements": missing_primary,
        "existing_assets": {
            "growth_closure_code": file_record(GROWTH_CLOSURE),
            "growth_closure_report": file_record(GROWTH_REPORT),
            "ppf_proxy_code": file_record(PPF_PROXY),
            "ppf_proxy_report": file_record(PPF_REPORT),
        },
        "audit_findings": [
            "The growth closure inserts freely chosen clustering, slip, sound-speed, gauge-range and helicity responses into a standard quasi-static matter-growth equation.",
            "The closure does not derive those responses from an SU2R action, perturbed stress-energy tensor or gauge-field constraints.",
            "The CAMB PPF run is an effective dark-energy-fluid proxy and does not implement non-Abelian SU2R perturbations.",
            "The historical PPF proxy fixed an ACT acoustic/nuisance anchor rather than fitting a complete source-background-matched SU2R parameter vector.",
            "Neither diagnostic can predict physical SU2R CMB lensing or parity spectra.",
        ],
        "decision": "BLOCK_PHYSICAL_SU2R_CMB",
    }
    json_path = RESULT_DIR / f"su2r_physical_readiness_audit_{DATE}.json"
    json_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    lines = [
        "# SU2R Physical Perturbation Readiness Audit",
        "",
        "Date: 18 July 2026",
        "",
        "## Decision",
        "",
        "**BLOCK PHYSICAL SU2R CMB/LENSING INTERPRETATION.**",
        "",
        f"The adapter refusal gate passed: `{refusal_verified}`. The existing growth closure and CAMB PPF proxy remain useful diagnostics, but neither supplies the field-level equations required for a physical SU2R spectrum.",
        "",
        "## Existing evidence",
        "",
        "- The 66,150-row growth scan found only a `Delta chi2=-0.0971` raw improvement over its smooth closure; parameter penalties retained the smooth branch.",
        "- The historical ACT PPF proxy strongly penalised its effective-w row, but it was an unmatched fluid proxy and is not a physical SU2R exclusion.",
        "",
        "## Missing lensing requirements",
        "",
        *[f"- `{item}`" for item in missing_lensing],
        "",
        "## Required next input from the model",
        "",
        "Provide the action or complete field equations; define the background SU2R field; derive scalar dynamical and constraint equations in a stated gauge; give delta-rho, delta-p, velocity and anisotropic-stress terms; specify adiabatic/isocurvature initial conditions; and derive the Weyl-potential relation used for lensing. The implementation module must then be hash-locked in the registry before ACT or Planck evaluation.",
    ]
    report_path = RESULT_DIR / f"su2r_physical_readiness_audit_report_{DATE}.md"
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Saved: {report_path}")


if __name__ == "__main__":
    main()
