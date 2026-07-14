r"""
Export Pantheon+SH0ES model residual series for HHT resonance diagnostics.

This helper creates a real CSV input for diagnose_hht_resonance.py --mode csv.
It evaluates selected best-fit cosmology branches on the SN redshift grid and
writes residual columns ordered by log(1+z).

Outputs:
    plamb_runs/diagnostics/hht_resonance/sn_residuals_for_hht.csv
    plamb_runs/diagnostics/hht_resonance/sn_residuals_for_hht_summary.md

Example:
    python export_sn_residuals_for_hht.py
    python diagnose_hht_resonance.py --mode csv --input-csv plamb_runs/diagnostics/hht_resonance/sn_residuals_for_hht.csv --x-col log1pz --y-col residual_SU2_baocov_4D --label SU2_residual
"""

from __future__ import annotations

import argparse
import ast
import csv
import importlib.util
import math
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parent
PROGRAM = ROOT / "PLamb_Test_10V6c_plus.py"
DATA = ROOT / "Pantheon+SH0ES.dat"
OUTDIR = ROOT / "plamb_runs" / "diagnostics" / "hht_resonance"


BASE_MODELS = [
    {
        "label": "LCDM_baocov",
        "model": "LCDM",
        "path": ROOT / "plamb_runs" / "next_su2_tests" / "runs" / "baocov_LCDM_SN_BAO_Planck" / "LCDM_emcee_bestfit.txt",
        "note": "Matched LCDM SN+BAO covariance+Planck baseline.",
    },
    {
        "label": "SU2_baocov_4D",
        "model": "SU2",
        "path": ROOT / "plamb_runs" / "next_su2_tests" / "runs" / "baocov_SU2_SN_BAO_Planck" / "SU2_emcee_bestfit.txt",
        "note": "Matched SU2 4D SN+BAO covariance+Planck branch.",
    },
    {
        "label": "SU2R_baocov_4D",
        "model": "SU2R",
        "path": ROOT / "plamb_runs" / "next_su2_tests" / "runs" / "baocov_SU2R_SN_BAO_Planck" / "SU2R_emcee_bestfit.txt",
        "note": "Matched reparameterized SU2R 4D branch.",
    },
    {
        "label": "SU2R_long_baocov",
        "model": "SU2R",
        "path": ROOT / "plamb_runs" / "next_su2_tests" / "runs" / "baocov_SU2R_long_SN_BAO_Planck" / "SU2R_emcee_bestfit.txt",
        "note": "Long SU2R chain, included when available.",
    },
    {
        "label": "SU2_wa_free_baocov",
        "model": "SU2",
        "path": ROOT / "plamb_runs" / "next_su2_tests" / "runs" / "baocov_SU2_wa_free_SN_BAO_Planck" / "SU2_emcee_bestfit.txt",
        "note": "Matched SU2 run with wa_chi free.",
    },
]


PLAMB_MODELS = [
    {
        "label": "PLAMB_H0free_alpha_free",
        "model": "PLAMB",
        "path": ROOT / "plamb_runs" / "plamb_noexp_pilot" / "runs" / "PLAMB_SN_BAO_H0free_alpha_free" / "PLAMB_emcee_bestfit.txt",
        "note": "PLAMB SN+BAO with H0 and alpha free.",
    },
    {
        "label": "PLAMB_SN_only_H0free",
        "model": "PLAMB",
        "path": ROOT / "plamb_runs" / "plamb_noexp_pilot" / "runs" / "PLAMB_SN_only_H0free" / "PLAMB_emcee_bestfit.txt",
        "note": "PLAMB SN-only with H0 free.",
    },
]


def import_plamb_module():
    spec = importlib.util.spec_from_file_location("plamb_fit", PROGRAM)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {PROGRAM}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_bestfit(path: Path) -> dict[str, Any]:
    values: dict[str, Any] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "=" not in line:
            continue
        key, raw = line.split("=", 1)
        key = key.strip()
        raw = raw.strip()
        try:
            values[key] = ast.literal_eval(raw)
        except Exception:
            try:
                values[key] = float(raw)
            except Exception:
                values[key] = raw
    return values


def fval(values: dict[str, Any], key: str, default: float = 0.0) -> float:
    try:
        raw = values.get(key, default)
        if raw in ("", None):
            return default
        return float(raw)
    except Exception:
        return default


def safe_col(label: str) -> str:
    return "".join(ch if ch.isalnum() or ch == "_" else "_" for ch in label)


def weighted_rms(resid: np.ndarray, sigma: np.ndarray) -> float:
    weights = 1.0 / np.maximum(sigma, 1e-30) ** 2
    return float(math.sqrt(np.sum(weights * resid * resid) / np.sum(weights)))


def evaluate_model(module, data: dict[str, np.ndarray], spec: dict[str, Any], nint: int) -> tuple[np.ndarray, dict[str, Any]]:
    vals = parse_bestfit(spec["path"])
    z = data["z"]
    mu_th = module.model_mu(
        z,
        spec["model"],
        fval(vals, "H0"),
        fval(vals, "Om"),
        fval(vals, "Ol"),
        fval(vals, "A_acc"),
        fval(vals, "n_acc"),
        fval(vals, "gamma_c"),
        fval(vals, "epsilon_M"),
        fval(vals, "Omega_chi0"),
        fval(vals, "w0_chi", -1.0),
        fval(vals, "wa_chi"),
        fval(vals, "plamb_alpha", 0.5),
        "linear_z",
        nint,
    )
    resid = np.asarray(mu_th, dtype=float) - np.asarray(data["mu"], dtype=float)
    return resid, vals


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_summary(
    path: Path,
    output_csv: Path,
    used: list[dict[str, Any]],
    skipped: list[dict[str, Any]],
    metrics: list[dict[str, Any]],
) -> None:
    lines = [
        "# SN residual export for HHT",
        "",
        f"CSV: `{output_csv}`",
        "",
        "## Included models",
        "",
        "| label | model | residual column | RMS [mag] | weighted RMS [mag] | note |",
        "|---|---|---|---:|---:|---|",
    ]
    for metric in metrics:
        lines.append(
            f"| {metric['label']} | {metric['model']} | `{metric['column']}` | "
            f"{metric['rms']:.5g} | {metric['weighted_rms']:.5g} | {metric['note']} |"
        )
    if skipped:
        lines.extend(["", "## Skipped models", "", "| label | reason |", "|---|---|"])
        for item in skipped:
            lines.append(f"| {item['label']} | {item['reason']} |")
    lines.extend(
        [
            "",
            "## Example HHT commands",
            "",
            "```powershell",
            f"python C:\\Users\\clive\\Documents\\Cosmology\\diagnose_hht_resonance.py --mode csv --input-csv {output_csv} --x-col log1pz --y-col residual_SU2_baocov_4D --label SU2_residual",
            f"python C:\\Users\\clive\\Documents\\Cosmology\\diagnose_hht_resonance.py --mode csv --input-csv {output_csv} --x-col log1pz --y-col delta_residual_SU2_baocov_4D_minus_LCDM_baocov --label SU2_minus_LCDM_residual",
            "```",
            "",
            "Use these as exploratory diagnostics only. Any apparent residual oscillation should be checked against shuffled-redshift controls and multiple model baselines.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Export SN residual columns for HHT diagnostics.")
    parser.add_argument("--include-plamb", action="store_true", help="Include PLAMB residual columns as well as LCDM/SU2.")
    parser.add_argument("--nint", type=int, default=200, help="Integration panel base passed to PLamb model_mu.")
    args = parser.parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)
    module = import_plamb_module()
    data = module.load_pantheon(str(DATA))
    order = np.argsort(data["z"])
    z = np.asarray(data["z"], dtype=float)[order]
    mu = np.asarray(data["mu"], dtype=float)[order]
    sigma = np.asarray(data["sigma"], dtype=float)[order]
    sorted_data = {"z": z, "mu": mu, "sigma": sigma}

    specs = list(BASE_MODELS)
    if args.include_plamb:
        specs.extend(PLAMB_MODELS)

    residual_columns: dict[str, np.ndarray] = {}
    used: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []
    metrics: list[dict[str, Any]] = []

    for spec in specs:
        if not spec["path"].exists():
            skipped.append({"label": spec["label"], "reason": f"missing {spec['path']}"})
            continue
        resid, vals = evaluate_model(module, sorted_data, spec, args.nint)
        col = "residual_" + safe_col(spec["label"])
        residual_columns[col] = resid
        used.append(spec)
        metrics.append(
            {
                "label": spec["label"],
                "model": spec["model"],
                "column": col,
                "rms": float(np.sqrt(np.mean(resid * resid))),
                "weighted_rms": weighted_rms(resid, sigma),
                "note": spec["note"],
                "H0": fval(vals, "H0"),
                "Om": fval(vals, "Om"),
                "Omega_chi0": fval(vals, "Omega_chi0"),
                "w0_chi": fval(vals, "w0_chi", -1.0),
            }
        )

    if not residual_columns:
        raise RuntimeError("No residual columns were created; check best-fit paths.")

    if "residual_LCDM_baocov" in residual_columns:
        lcdm = residual_columns["residual_LCDM_baocov"]
        for col, values in list(residual_columns.items()):
            if col == "residual_LCDM_baocov":
                continue
            delta_col = "delta_" + col + "_minus_LCDM_baocov"
            residual_columns[delta_col] = values - lcdm

    rows: list[dict[str, Any]] = []
    for i in range(z.size):
        row: dict[str, Any] = {
            "index_sorted": i,
            "z": float(z[i]),
            "log1pz": float(np.log1p(z[i])),
            "mu_obs": float(mu[i]),
            "sigma_mu": float(sigma[i]),
        }
        for col, values in residual_columns.items():
            row[col] = float(values[i])
        rows.append(row)

    output_csv = OUTDIR / "sn_residuals_for_hht.csv"
    summary_md = OUTDIR / "sn_residuals_for_hht_summary.md"
    write_csv(output_csv, rows)
    write_summary(summary_md, output_csv, used, skipped, metrics)

    print(f"Saved residual CSV: {output_csv}")
    print(f"Saved summary: {summary_md}")
    print("Available HHT y-columns:")
    for col in residual_columns:
        print(f"  {col}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
