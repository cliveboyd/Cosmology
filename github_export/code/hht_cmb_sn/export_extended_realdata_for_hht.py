r"""
Export extended real-data pull series for HHT resonance diagnostics.

This script expands the HHT real-data test beyond SN-only residuals. It builds
standardized pull series from:

    * Pantheon+SH0ES SN residuals, binned in log(1+z)
    * BAO DM/DH residual pulls from bao_MASTER_LONG.csv and bao_MASTER_cov.txt
    * Cosmic chronometer H(z) residual pulls from chronometers/data_CC.dat

The output is deliberately a diagnostic screening series, not a formal joint
likelihood. The probes have different observables and covariances, so the
script converts each row to a normalized pull before combining.

Outputs:
    plamb_runs/diagnostics/hht_resonance/extended_realdata_pulls_for_hht.csv
    plamb_runs/diagnostics/hht_resonance/extended_realdata_pulls_for_hht_summary.md

Example:
    python export_extended_realdata_for_hht.py
    python diagnose_hht_resonance.py --mode csv --input-csv plamb_runs/diagnostics/hht_resonance/extended_realdata_pulls_for_hht.csv --x-col x_hht --y-col joint_delta_pull_SU2_baocov_4D_minus_LCDM_baocov --label extended_SU2_minus_LCDM
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
SN_DATA = ROOT / "Pantheon+SH0ES.dat"
BAO_DATA = ROOT / "bao_MASTER_LONG.csv"
BAO_COV = ROOT / "bao_MASTER_cov.txt"
CC_DATA = ROOT / "chronometers" / "data_CC.dat"
OUTDIR = ROOT / "plamb_runs" / "diagnostics" / "hht_resonance"

RD = 147.09
NINT = 200


MODELS = [
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
        "label": "SU2_wa_free_baocov",
        "model": "SU2",
        "path": ROOT / "plamb_runs" / "next_su2_tests" / "runs" / "baocov_SU2_wa_free_SN_BAO_Planck" / "SU2_emcee_bestfit.txt",
        "note": "Matched SU2 run with wa_chi free.",
    },
    {
        "label": "SU2R_long_baocov",
        "model": "SU2R",
        "path": ROOT / "plamb_runs" / "next_su2_tests" / "runs" / "baocov_SU2R_long_SN_BAO_Planck" / "SU2R_emcee_bestfit.txt",
        "note": "Long SU2R chain, included when available.",
    },
]


PROBE_OFFSETS = {
    "SN_bin": 0.00000,
    "CC": 0.00008,
    "BAO_DH_over_rd": 0.00016,
    "BAO_DM_over_rd": 0.00024,
    "BAO_DV_over_rd": 0.00032,
}


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


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def load_cc_table(path: Path) -> dict[str, np.ndarray]:
    rows: list[tuple[float, float, float, str]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        text = line.strip()
        if not text or text.startswith("#"):
            continue
        parts = [p.strip() for p in text.split(",")]
        if len(parts) < 3:
            continue
        try:
            z = float(parts[0])
            h = float(parts[1])
            sigma = float(parts[2])
        except Exception:
            continue
        ref = parts[3] if len(parts) > 3 else ""
        rows.append((z, h, sigma, ref))
    if not rows:
        raise ValueError(f"No CC rows parsed from {path}")
    rows.sort(key=lambda r: r[0])
    return {
        "z": np.asarray([r[0] for r in rows], dtype=float),
        "H": np.asarray([r[1] for r in rows], dtype=float),
        "sigma": np.asarray([r[2] for r in rows], dtype=float),
        "ref": np.asarray([r[3] for r in rows], dtype=object),
    }


def model_mu(module, z: np.ndarray, spec: dict[str, Any], vals: dict[str, Any]) -> np.ndarray:
    return module.model_mu(
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
        NINT,
    )


def model_hz(module, z: np.ndarray, spec: dict[str, Any], vals: dict[str, Any]) -> np.ndarray:
    if spec["model"].upper() == "PLAMB":
        _dm, dh, _dv, _dl = module.plamb_noexp_distances(
            z,
            fval(vals, "H0"),
            fval(vals, "A_acc"),
            fval(vals, "n_acc"),
            fval(vals, "gamma_c"),
            nint_base=NINT,
        )
        return module.C_KMS / np.maximum(dh, 1e-12)
    ez = module.E_of_z(
        z,
        fval(vals, "Om"),
        fval(vals, "Ol"),
        fval(vals, "A_acc"),
        fval(vals, "n_acc"),
        fval(vals, "Omega_chi0"),
        fval(vals, "w0_chi", -1.0),
        fval(vals, "wa_chi"),
    )
    return fval(vals, "H0") * ez


def model_bao(module, bao: dict[str, Any], spec: dict[str, Any], vals: dict[str, Any]) -> np.ndarray:
    return module.predict_bao_vector(
        bao,
        fval(vals, "H0"),
        fval(vals, "Om"),
        fval(vals, "Ol"),
        fval(vals, "A_acc"),
        fval(vals, "n_acc"),
        fval(vals, "gamma_c"),
        NINT,
        RD,
        fval(vals, "Omega_chi0"),
        fval(vals, "w0_chi", -1.0),
        fval(vals, "wa_chi"),
        spec["model"],
    )


def available_models() -> tuple[list[dict[str, Any]], dict[str, dict[str, Any]], list[dict[str, str]]]:
    used: list[dict[str, Any]] = []
    vals_by_label: dict[str, dict[str, Any]] = {}
    skipped: list[dict[str, str]] = []
    for spec in MODELS:
        if not spec["path"].exists():
            skipped.append({"label": spec["label"], "reason": f"missing {spec['path']}"})
            continue
        used.append(spec)
        vals_by_label[spec["label"]] = parse_bestfit(spec["path"])
    return used, vals_by_label, skipped


def build_sn_bin_rows(
    module,
    used: list[dict[str, Any]],
    vals_by_label: dict[str, dict[str, Any]],
    n_bins: int,
    min_count: int,
) -> list[dict[str, Any]]:
    data = module.load_pantheon(str(SN_DATA))
    z = np.asarray(data["z"], dtype=float)
    mu = np.asarray(data["mu"], dtype=float)
    sigma = np.asarray(data["sigma"], dtype=float)
    x = np.log1p(z)

    model_pull: dict[str, np.ndarray] = {}
    for spec in used:
        label = spec["label"]
        pred = model_mu(module, z, spec, vals_by_label[label])
        model_pull[label] = (pred - mu) / np.maximum(sigma, 1e-30)

    edges = np.linspace(float(np.min(x)), float(np.max(x)), n_bins + 1)
    rows: list[dict[str, Any]] = []
    for i in range(n_bins):
        mask = (x >= edges[i]) & (x < edges[i + 1] if i < n_bins - 1 else x <= edges[i + 1])
        if int(np.sum(mask)) < min_count:
            continue
        weights = 1.0 / np.maximum(sigma[mask], 1e-30) ** 2
        x_mid = float(np.average(x[mask], weights=weights))
        z_mid = float(np.expm1(x_mid))
        row: dict[str, Any] = {
            "probe": "SN_bin",
            "kind": "mu_binned_pull",
            "z": z_mid,
            "log1pz": x_mid,
            "x_hht": x_mid + PROBE_OFFSETS["SN_bin"],
            "n_in_bin": int(np.sum(mask)),
            "observed": math.nan,
            "sigma": 1.0,
        }
        for spec in used:
            label = spec["label"]
            col = "joint_pull_" + safe_col(label)
            row[col] = float(np.average(model_pull[label][mask], weights=weights))
        rows.append(row)
    return rows


def build_bao_rows(
    module,
    used: list[dict[str, Any]],
    vals_by_label: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    bao = module.load_bao(str(BAO_DATA), str(BAO_COV))
    cov = np.asarray(bao["C"], dtype=float)
    sigma = np.sqrt(np.diag(cov))
    obs = np.asarray(bao["y"], dtype=float)
    z = np.asarray(bao["z"], dtype=float)
    kinds = [str(k) for k in bao["kind"]]

    preds_by_label: dict[str, np.ndarray] = {}
    for spec in used:
        label = spec["label"]
        preds_by_label[label] = model_bao(module, bao, spec, vals_by_label[label])

    rows: list[dict[str, Any]] = []
    for i, (zi, kind) in enumerate(zip(z, kinds)):
        probe = "BAO_" + kind
        x = float(np.log1p(zi) + PROBE_OFFSETS.get(probe, 0.0004))
        row: dict[str, Any] = {
            "probe": "BAO",
            "kind": kind,
            "z": float(zi),
            "log1pz": float(np.log1p(zi)),
            "x_hht": x,
            "n_in_bin": 1,
            "observed": float(obs[i]),
            "sigma": float(sigma[i]),
        }
        for spec in used:
            label = spec["label"]
            col = "joint_pull_" + safe_col(label)
            row[col] = float((preds_by_label[label][i] - obs[i]) / max(sigma[i], 1e-30))
        rows.append(row)
    return rows


def build_cc_rows(
    module,
    used: list[dict[str, Any]],
    vals_by_label: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    cc = load_cc_table(CC_DATA)
    z = np.asarray(cc["z"], dtype=float)
    obs = np.asarray(cc["H"], dtype=float)
    sigma = np.asarray(cc["sigma"], dtype=float)

    preds_by_label: dict[str, np.ndarray] = {}
    for spec in used:
        label = spec["label"]
        preds_by_label[label] = model_hz(module, z, spec, vals_by_label[label])

    rows: list[dict[str, Any]] = []
    for i, zi in enumerate(z):
        x0 = float(np.log1p(zi))
        row: dict[str, Any] = {
            "probe": "CC",
            "kind": "H_z",
            "z": float(zi),
            "log1pz": x0,
            "x_hht": x0 + PROBE_OFFSETS["CC"],
            "n_in_bin": 1,
            "observed": float(obs[i]),
            "sigma": float(sigma[i]),
        }
        for spec in used:
            label = spec["label"]
            col = "joint_pull_" + safe_col(label)
            row[col] = float((preds_by_label[label][i] - obs[i]) / max(sigma[i], 1e-30))
        rows.append(row)
    return rows


def add_delta_columns(rows: list[dict[str, Any]], used: list[dict[str, Any]]) -> None:
    lcdm_col = "joint_pull_LCDM_baocov"
    if not rows or lcdm_col not in rows[0]:
        return
    for spec in used:
        label = spec["label"]
        col = "joint_pull_" + safe_col(label)
        if col == lcdm_col or col not in rows[0]:
            continue
        delta_col = "joint_delta_pull_" + safe_col(label) + "_minus_LCDM_baocov"
        for row in rows:
            row[delta_col] = float(row[col]) - float(row[lcdm_col])


def write_summary(
    path: Path,
    output_csv: Path,
    rows: list[dict[str, Any]],
    used: list[dict[str, Any]],
    skipped: list[dict[str, str]],
) -> None:
    probe_counts: dict[str, int] = {}
    for row in rows:
        probe_counts[str(row["probe"])] = probe_counts.get(str(row["probe"]), 0) + 1

    z_values = [float(row["z"]) for row in rows]
    lines = [
        "# Extended real-data HHT export",
        "",
        f"CSV: `{output_csv}`",
        "",
        "## Scope",
        "",
        "This file combines standardized pull series from SN bins, BAO rows, and cosmic chronometer H(z) rows. It is intended for HHT screening only, not as a replacement for the formal likelihood.",
        "",
        f"- z range: {min(z_values):.5g} to {max(z_values):.5g}",
        f"- HHT rows: {len(rows)}",
    ]
    for probe, count in sorted(probe_counts.items()):
        lines.append(f"- {probe}: {count} rows")

    lines.extend(["", "## Included models", "", "| label | model | note |", "|---|---|---|"])
    for spec in used:
        lines.append(f"| {spec['label']} | {spec['model']} | {spec['note']} |")

    if skipped:
        lines.extend(["", "## Skipped models", "", "| label | reason |", "|---|---|"])
        for item in skipped:
            lines.append(f"| {item['label']} | {item['reason']} |")

    lines.extend(
        [
            "",
            "## Recommended HHT commands",
            "",
            "```powershell",
            f"python C:\\Users\\clive\\Documents\\Cosmology\\diagnose_hht_resonance.py --mode csv --input-csv {output_csv} --x-col x_hht --y-col joint_delta_pull_SU2_baocov_4D_minus_LCDM_baocov --label extended_SU2_minus_LCDM",
            f"python C:\\Users\\clive\\Documents\\Cosmology\\diagnose_hht_resonance.py --mode csv --input-csv {output_csv} --x-col x_hht --y-col joint_delta_pull_SU2R_baocov_4D_minus_LCDM_baocov --label extended_SU2R_minus_LCDM",
            f"python C:\\Users\\clive\\Documents\\Cosmology\\diagnose_hht_resonance.py --mode csv --input-csv {output_csv} --x-col x_hht --y-col joint_pull_SU2_baocov_4D --label extended_SU2_pull",
            "```",
            "",
            "## Caution",
            "",
            "The combined series mixes different probes after pull-standardization. A detected IMF would be a feature candidate only; it must be checked probe-by-probe and against shuffled-order controls.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Export extended real-data pull series for HHT diagnostics.")
    parser.add_argument("--sn-bins", type=int, default=60, help="Number of log(1+z) SN bins before min-count filtering.")
    parser.add_argument("--sn-min-count", type=int, default=8, help="Minimum SNe per retained SN bin.")
    args = parser.parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)
    module = import_plamb_module()
    used, vals_by_label, skipped = available_models()
    if not used:
        raise RuntimeError("No model best-fit files available.")

    rows: list[dict[str, Any]] = []
    rows.extend(build_sn_bin_rows(module, used, vals_by_label, args.sn_bins, args.sn_min_count))
    rows.extend(build_bao_rows(module, used, vals_by_label))
    rows.extend(build_cc_rows(module, used, vals_by_label))
    add_delta_columns(rows, used)
    rows.sort(key=lambda r: (float(r["x_hht"]), str(r["probe"]), str(r["kind"])))

    output_csv = OUTDIR / "extended_realdata_pulls_for_hht.csv"
    summary_md = OUTDIR / "extended_realdata_pulls_for_hht_summary.md"
    write_csv(output_csv, rows)
    write_summary(summary_md, output_csv, rows, used, skipped)

    print(f"Saved extended real-data pulls: {output_csv}")
    print(f"Saved summary: {summary_md}")
    print(f"Rows: {len(rows)}  z=[{min(float(r['z']) for r in rows):.5g}, {max(float(r['z']) for r in rows):.5g}]")
    print("Available HHT y-columns:")
    for key in rows[0]:
        if key.startswith("joint_pull_") or key.startswith("joint_delta_pull_"):
            print(f"  {key}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
