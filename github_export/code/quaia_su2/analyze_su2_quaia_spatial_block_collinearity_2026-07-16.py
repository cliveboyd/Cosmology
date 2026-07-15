from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import shutil
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np


ROOT = Path(r"C:\Users\clive\Documents\Cosmology")
CODE_DIR = ROOT / "github_export" / "code" / "quaia_su2"
OUTPUTS = Path(r"C:\Users\clive\Documents\Codex\2026-07-14\ok-continu\outputs")
DEFAULT_OUT = ROOT / "plamb_runs" / "diagnostics" / "su2_quaia_spatial_block_collinearity_20260716"
INJECTION_SCRIPT = "analyze_su2_quaia_injection_recovery_2026-07-15.py"
CMB_L_DEG = 264.021
CMB_B_DEG = 48.253


def load_module(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


INJ = load_module("su2_quaia_injection_recovery", CODE_DIR / INJECTION_SCRIPT)


def parse_csv_floats(text: str) -> list[float]:
    return [float(part.strip()) for part in text.split(",") if part.strip()]


def fmt(value: Any) -> str:
    if isinstance(value, float):
        return f"{value:.6g}" if math.isfinite(value) else ""
    return str(value)


def markdown_table(rows: list[dict[str, Any]], columns: list[str]) -> str:
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(fmt(row.get(column, "")) for column in columns) + " |")
    return "\n".join(lines)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def p_ge(values: np.ndarray, observed: float) -> float:
    return float((1 + np.sum(values >= observed)) / (len(values) + 1))


def unit_from_lb(l_deg: float, b_deg: float) -> np.ndarray:
    l_rad = math.radians(l_deg)
    b_rad = math.radians(b_deg)
    cos_b = math.cos(b_rad)
    return np.array([cos_b * math.cos(l_rad), cos_b * math.sin(l_rad), math.sin(b_rad)], dtype=float)


def standardize_columns(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    center = np.nanmedian(values, axis=0)
    scale = np.nanstd(values, axis=0)
    scale = np.where(np.isfinite(scale) & (scale > 0.0), scale, 1.0)
    return (values - center) / scale


def direction_collinearity(projection: np.ndarray, controls: np.ndarray) -> tuple[float, float, float]:
    y = np.asarray(projection, dtype=float)
    y = y - np.mean(y)
    x = np.column_stack([np.ones(len(y)), standardize_columns(controls)])
    beta = np.linalg.lstsq(x, y, rcond=None)[0]
    residual = y - x @ beta
    total = float(np.dot(y, y))
    r2 = 1.0 - float(np.dot(residual, residual)) / max(total, 1e-300)
    vif = 1.0 / max(1.0 - r2, 1e-12)
    y_sd = float(np.std(y))
    if y_sd <= 0.0:
        max_abs_corr = float("nan")
    else:
        y_std = y / y_sd
        corr = np.mean(standardize_columns(controls) * y_std[:, None], axis=0)
        max_abs_corr = float(np.max(np.abs(corr)))
    return float(r2), float(vif), max_abs_corr


def canonical_correlation(dipole: np.ndarray, controls: np.ndarray) -> float:
    d = standardize_columns(dipole)
    c = standardize_columns(controls)
    qd = np.linalg.qr(d, mode="reduced")[0]
    qc = np.linalg.qr(c, mode="reduced")[0]
    singular = np.linalg.svd(qd.T @ qc, compute_uv=False)
    return float(np.clip(singular[0], 0.0, 1.0))


def equal_area_cells(l_deg: np.ndarray, b_deg: np.ndarray, cell_deg: float) -> tuple[np.ndarray, int, int, int]:
    n_lon = max(int(round(360.0 / cell_deg)), 4)
    n_sin_b = max(int(round(180.0 / cell_deg)), 2)
    lon_bin = np.floor((np.mod(l_deg, 360.0) / 360.0) * n_lon).astype(int)
    sin_b = np.sin(np.radians(b_deg))
    sin_bin = np.floor(((sin_b + 1.0) / 2.0) * n_sin_b).astype(int)
    lon_bin = np.clip(lon_bin, 0, n_lon - 1)
    sin_bin = np.clip(sin_bin, 0, n_sin_b - 1)
    raw = sin_bin * n_lon + lon_bin
    unique, inverse = np.unique(raw, return_inverse=True)
    return inverse.astype(np.int64), len(unique), n_lon, n_sin_b


def prepare_cut(
    q: dict[str, np.ndarray],
    control_values: list[np.ndarray],
    bcut: float,
    z_min: float,
    z_max: float,
) -> dict[str, Any]:
    base = (q["z"] >= z_min) & (q["z"] < z_max) & (np.abs(q["b"]) >= bcut)
    finite = base & np.isfinite(q["z"])
    finite &= np.all(np.isfinite(q["vec"]), axis=1)
    for values in control_values:
        finite &= np.isfinite(values)
    idx = np.flatnonzero(finite)
    z = q["z"][idx].astype(float)
    vec = q["vec"][idx].astype(float)
    controls = np.column_stack([values[idx] for values in control_values]).astype(float)
    x_control = np.column_stack([np.ones(len(idx)), standardize_columns(controls)])
    control_beta = np.linalg.lstsq(x_control, z, rcond=None)[0]
    residual = z - x_control @ control_beta
    control_total = float(np.dot(z - np.mean(z), z - np.mean(z)))
    control_r2 = 1.0 - float(np.dot(residual, residual)) / max(control_total, 1e-300)

    x_dipole = np.column_stack([np.ones(len(idx)), vec])
    xtx_inv = np.linalg.pinv(x_dipole.T @ x_dipole)
    beta = xtx_inv @ (x_dipole.T @ residual)
    observed_fit_residual = residual - x_dipole @ beta
    dof = max(len(idx) - x_dipole.shape[1], 1)
    sigma2 = float(np.dot(observed_fit_residual, observed_fit_residual) / dof)
    dipole = beta[1:4]
    dipole_cov = xtx_inv[1:4, 1:4]
    info = np.linalg.pinv(dipole_cov)
    observed_snr = math.sqrt(max(float(dipole @ info @ dipole) / max(sigma2, 1e-300), 0.0))
    observed_amp = float(np.linalg.norm(dipole))
    l_deg, b_deg = INJ.unit_to_lb(dipole)

    full_design = np.column_stack([np.ones(len(idx)), standardize_columns(vec), standardize_columns(controls)])
    condition_number = float(np.linalg.cond(full_design))
    return {
        "bcut_deg": bcut,
        "idx": idx,
        "z": z,
        "vec": vec,
        "l": q["l"][idx].astype(float),
        "b": q["b"][idx].astype(float),
        "controls": controls,
        "residual": residual,
        "x_dipole": x_dipole,
        "xtx_inv": xtx_inv,
        "dipole_info": info,
        "dof": dof,
        "yty": float(np.dot(residual, residual)),
        "observed_snr": observed_snr,
        "observed_amp": observed_amp,
        "observed_l_deg": l_deg,
        "observed_b_deg": b_deg,
        "control_r2": control_r2,
        "canonical_corr": canonical_correlation(vec, controls),
        "condition_number": condition_number,
    }


def add_cluster_summaries(prepared: list[dict[str, Any]], q: dict[str, np.ndarray], cell_deg: float) -> tuple[int, int, int]:
    target = (q["z"] >= min(item["z"].min() for item in prepared)) & (q["z"] <= max(item["z"].max() for item in prepared))
    global_idx = np.flatnonzero(target & np.isfinite(q["l"]) & np.isfinite(q["b"]))
    global_cells, n_cells, n_lon, n_sin_b = equal_area_cells(q["l"][global_idx], q["b"][global_idx], cell_deg)
    cell_map = {int(source_idx): int(cell) for source_idx, cell in zip(global_idx, global_cells)}
    for item in prepared:
        cell = np.array([cell_map[int(source_idx)] for source_idx in item["idx"]], dtype=np.int64)
        aggregate = np.zeros((item["x_dipole"].shape[1], n_cells), dtype=float)
        xr = item["x_dipole"] * item["residual"][:, None]
        for column in range(xr.shape[1]):
            aggregate[column] = np.bincount(cell, weights=xr[:, column], minlength=n_cells)
        item["cluster_xt_residual"] = aggregate
        item["n_occupied_cells"] = int(len(np.unique(cell)))
    return n_cells, n_lon, n_sin_b


def fit_cluster_batches(item: dict[str, Any], signs: np.ndarray) -> np.ndarray:
    xty = item["cluster_xt_residual"] @ signs
    beta = item["xtx_inv"] @ xty
    rss = np.maximum(item["yty"] - np.sum(beta * xty, axis=0), 0.0)
    sigma2 = rss / item["dof"]
    dipole = beta[1:4]
    snr2 = np.einsum("ib,ij,jb->b", dipole, item["dipole_info"], dipole) / np.maximum(sigma2, 1e-300)
    return np.sqrt(np.maximum(snr2, 0.0))


def run_spatial_null(
    prepared: list[dict[str, Any]],
    q: dict[str, np.ndarray],
    cell_deg: float,
    n_mocks: int,
    batch_size: int,
    seed: int,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    n_cells, n_lon, n_sin_b = add_cluster_summaries(prepared, q, cell_deg)
    rng = np.random.default_rng(seed)
    maxima = np.empty(n_mocks, dtype=float)
    cursor = 0
    while cursor < n_mocks:
        batch = min(batch_size, n_mocks - cursor)
        signs = rng.choice(np.array([-1.0, 1.0]), size=(n_cells, batch), replace=True)
        cut_snrs = np.vstack([fit_cluster_batches(item, signs) for item in prepared])
        maxima[cursor : cursor + batch] = np.max(cut_snrs, axis=0)
        cursor += batch

    observed_max = max(item["observed_snr"] for item in prepared)
    row = {
        "cell_deg": cell_deg,
        "n_cells": n_cells,
        "n_lon_bins": n_lon,
        "n_sin_b_bins": n_sin_b,
        "n_mocks": n_mocks,
        "observed_family_max_snr": observed_max,
        "null_max_snr_mean": float(np.mean(maxima)),
        "null_max_snr_p95": float(np.quantile(maxima, 0.95)),
        "null_max_snr_p99": float(np.quantile(maxima, 0.99)),
        "family_p": p_ge(maxima, observed_max),
        "monte_carlo_floor": 1.0 / (n_mocks + 1),
    }
    for item in prepared:
        row[f"occupied_cells_b{int(item['bcut_deg'])}"] = item["n_occupied_cells"]
    return [row], {"maxima": maxima}


def build_collinearity_rows(
    prepared: list[dict[str, Any]],
    direction_map: dict[str, np.ndarray],
    control_names: list[str],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    summary_rows: list[dict[str, Any]] = []
    control_rows: list[dict[str, Any]] = []
    for item in prepared:
        controls_z = standardize_columns(item["controls"])
        for label, unit in direction_map.items():
            projection = item["vec"] @ unit
            r2, vif, max_abs_corr = direction_collinearity(projection, item["controls"])
            correlations = np.mean(controls_z * ((projection - np.mean(projection)) / np.std(projection))[:, None], axis=0)
            max_index = int(np.argmax(np.abs(correlations)))
            summary_rows.append(
                {
                    "direction_label": label,
                    "bcut_deg": item["bcut_deg"],
                    "N": len(item["idx"]),
                    "projection_control_r2": r2,
                    "projection_vif": vif,
                    "max_abs_control_corr": max_abs_corr,
                    "max_corr_control": control_names[max_index],
                    "max_corr_signed": float(correlations[max_index]),
                    "dipole_control_canonical_corr": item["canonical_corr"],
                    "full_design_condition_number": item["condition_number"],
                    "control_only_r2_z": item["control_r2"],
                }
            )
            for name, corr in zip(control_names, correlations):
                control_rows.append(
                    {
                        "direction_label": label,
                        "bcut_deg": item["bcut_deg"],
                        "control": name,
                        "projection_correlation": float(corr),
                    }
                )
    return summary_rows, control_rows


def write_reports(
    out_dir: Path,
    prepared: list[dict[str, Any]],
    collinearity: list[dict[str, Any]],
    spatial_rows: list[dict[str, Any]],
    config: dict[str, Any],
) -> tuple[Path, Path]:
    target = min(collinearity, key=lambda row: (0 if row["direction_label"] == "locked_b15" else 1, abs(float(row["bcut_deg"]) - 15.0)))
    all_pass = all(float(row["family_p"]) <= 0.01 for row in spatial_rows)
    if all_pass:
        bottom = (
            "The z=1.0-1.5 residual family remains unusual under every tested spatial-cluster scale. "
            "This strengthens its status as a follow-up candidate, but it is still not a discovery because this is a targeted-family test and the target was selected from earlier scans."
        )
    else:
        bottom = (
            "The z=1.0-1.5 residual family does not remain below the 1% gate at every spatial-cluster scale. "
            "The earlier object-shuffle p-value is therefore not robust enough for promotion."
        )
    methods = (
        "For each latitude cut, redshift was regressed on the seven Gaia scan-law and catalogue-colour controls. "
        "Control residuals were multiplied by shared Rademacher signs on an equal-area Galactic longitude x sin(latitude) grid. "
        "The same cell signs were used for all nested latitude cuts in each mock, preserving their covariance and within-cell residual structure."
    )
    observed_rows = [
        {
            "bcut_deg": item["bcut_deg"],
            "N": len(item["idx"]),
            "residual_snr": item["observed_snr"],
            "residual_amp": item["observed_amp"],
            "residual_l_deg": item["observed_l_deg"],
            "residual_b_deg": item["observed_b_deg"],
            "control_r2_z": item["control_r2"],
            "canonical_corr": item["canonical_corr"],
            "condition_number": item["condition_number"],
        }
        for item in prepared
    ]
    report = [
        "# SU2 / Quaia Spatial-Block and Template-Collinearity Audit",
        "",
        "Date: 2026-07-16",
        "",
        "## Bottom Line",
        "",
        bottom,
        "",
        "## Why This Replaces the Simple Target-Family Shuffle Readout",
        "",
        "The earlier residual scan used 1,000 object-level shuffles and independently permuted overlapping windows. Its p=0.000999 value was the Monte Carlo floor, not a resolved probability. The present audit addresses the targeted z=1.0-1.5 family with shared spatial blocks and a substantially larger null ensemble.",
        "",
        "## Spatial Null",
        "",
        methods,
        "",
        markdown_table(
            spatial_rows,
            ["cell_deg", "n_cells", "n_mocks", "observed_family_max_snr", "null_max_snr_p95", "null_max_snr_p99", "family_p", "monte_carlo_floor"],
        ),
        "",
        "## Observed Controlled Dipoles",
        "",
        markdown_table(
            observed_rows,
            ["bcut_deg", "N", "residual_snr", "residual_amp", "residual_l_deg", "residual_b_deg", "control_r2_z", "canonical_corr", "condition_number"],
        ),
        "",
        "## Direction-Specific Template Collinearity",
        "",
        markdown_table(
            collinearity,
            ["direction_label", "bcut_deg", "projection_control_r2", "projection_vif", "max_corr_control", "max_corr_signed", "dipole_control_canonical_corr", "full_design_condition_number"],
        ),
        "",
        "## Interpretation",
        "",
        f"For the locked b=15-degree direction at the matching cut, the control-template R2 is {float(target['projection_control_r2']):.4g} (VIF {float(target['projection_vif']):.4g}); its largest single-template correlation is {float(target['max_corr_signed']):.4g} with {target['max_corr_control']}.",
        "",
        "The collinearity table is diagnostic rather than a correction. A high canonical correlation or direction-specific VIF explains why unconstrained controls can rotate or inflate recovered injections; it does not prove that the angular mode is physical.",
        "",
        "## Claim Boundary",
        "",
        "Treat the z=1.0-1.5 mode as a targeted follow-up candidate only. Promotion requires an independently specified catalogue or release, a pre-registered angular statistic, and a survey-selection likelihood whose injection power and false-positive rate are calibrated before examining the new data.",
        "",
        "## Configuration",
        "",
        "```json",
        json.dumps(config, indent=2),
        "```",
    ]
    readout = [
        "# SU2 / Quaia Spatial-Block and Collinearity Readout",
        "",
        "Date: 2026-07-16",
        "",
        bottom,
        "",
        markdown_table(
            spatial_rows,
            ["cell_deg", "n_mocks", "observed_family_max_snr", "null_max_snr_p99", "family_p"],
        ),
        "",
        f"Locked b=15 direction at |b|>15: template R2 `{float(target['projection_control_r2']):.4g}`, VIF `{float(target['projection_vif']):.4g}`, strongest control `{target['max_corr_control']}` with signed correlation `{float(target['max_corr_signed']):.4g}`.",
        "",
        "Claim boundary: targeted follow-up candidate only; no SU2 detection or headline anisotropy claim.",
    ]
    report_path = out_dir / "su2_quaia_spatial_block_collinearity_report.md"
    readout_path = out_dir / "su2_quaia_spatial_block_collinearity_readout.md"
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")
    readout_path.write_text("\n".join(readout) + "\n", encoding="utf-8")
    return report_path, readout_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Spatial-block null and template-collinearity audit for the SU2/Quaia z=1.0-1.5 family.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--quaia", type=Path, default=INJ.EXT.DEFAULT_QUAIA)
    parser.add_argument("--randoms", type=Path, default=INJ.EXT.DEFAULT_RANDOMS)
    parser.add_argument("--selection", type=Path, default=INJ.EXT.DEFAULT_SELECTION)
    parser.add_argument("--sfd-dir", type=Path, default=INJ.EXT.DEFAULT_SFD_DIR)
    parser.add_argument("--gaiascanlaw-data-dir", type=Path, default=None)
    parser.add_argument("--z-min", type=float, default=1.0)
    parser.add_argument("--z-max", type=float, default=1.5)
    parser.add_argument("--b-cuts", type=float, nargs="+", default=[10.0, 15.0, 20.0, 25.0, 30.0, 35.0])
    parser.add_argument("--cell-degrees", default="8,12,16")
    parser.add_argument("--n-mocks", type=int, default=50000)
    parser.add_argument("--batch-size", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=200716)
    parser.add_argument("--copy-to-outputs", action=argparse.BooleanOptionalAction, default=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    load_args = SimpleNamespace(
        quaia=args.quaia,
        randoms=args.randoms,
        selection=args.selection,
        sfd_dir=args.sfd_dir,
        gaiascanlaw_data_dir=args.gaiascanlaw_data_dir,
        z_min=args.z_min,
        z_max=args.z_max,
        b_cuts=args.b_cuts,
    )
    q, load_meta = INJ.load_data(load_args)
    control_pairs = INJ.scanlaw_colour_controls(q)
    control_names = [name for name, _values in control_pairs]
    control_values = [values for _name, values in control_pairs]
    prepared = [prepare_cut(q, control_values, bcut, args.z_min, args.z_max) for bcut in args.b_cuts]

    raw_15 = INJ.fit_single_dipole(prepared[args.b_cuts.index(15.0)]["z"], prepared[args.b_cuts.index(15.0)]["vec"])
    raw_25 = INJ.fit_single_dipole(prepared[args.b_cuts.index(25.0)]["z"], prepared[args.b_cuts.index(25.0)]["vec"])
    direction_map = {
        "locked_b15": raw_15["dvec"] / raw_15["amp"],
        "locked_b25": raw_25["dvec"] / raw_25["amp"],
        "cmb": unit_from_lb(CMB_L_DEG, CMB_B_DEG),
        "anti_cmb": -unit_from_lb(CMB_L_DEG, CMB_B_DEG),
    }
    collinearity_rows, control_rows = build_collinearity_rows(prepared, direction_map, control_names)

    spatial_rows: list[dict[str, Any]] = []
    for index, cell_deg in enumerate(parse_csv_floats(args.cell_degrees)):
        rows, _details = run_spatial_null(
            prepared,
            q,
            cell_deg,
            args.n_mocks,
            args.batch_size,
            args.seed + index * 100003,
        )
        spatial_rows.extend(rows)
        print(f"cell_deg={cell_deg:g} family_p={rows[0]['family_p']:.6g}", flush=True)

    config = {
        "date": "2026-07-16",
        "analysis": "su2_quaia_spatial_block_collinearity",
        "z_min": args.z_min,
        "z_max": args.z_max,
        "b_cuts": args.b_cuts,
        "cell_degrees": parse_csv_floats(args.cell_degrees),
        "n_mocks_per_cell_scale": args.n_mocks,
        "batch_size": args.batch_size,
        "seed": args.seed,
        "controls": control_names,
        "null": "Shared Rademacher wild bootstrap on equal-area Galactic longitude x sin(latitude) cells; common signs across nested latitude cuts.",
        "load_meta": load_meta,
    }
    observed_rows = [
        {
            "bcut_deg": item["bcut_deg"],
            "N": len(item["idx"]),
            "residual_snr": item["observed_snr"],
            "residual_amp": item["observed_amp"],
            "residual_l_deg": item["observed_l_deg"],
            "residual_b_deg": item["observed_b_deg"],
            "control_r2_z": item["control_r2"],
            "canonical_corr": item["canonical_corr"],
            "condition_number": item["condition_number"],
        }
        for item in prepared
    ]
    write_csv(args.out_dir / "su2_quaia_spatial_block_summary.csv", spatial_rows)
    write_csv(args.out_dir / "su2_quaia_template_collinearity_summary.csv", collinearity_rows)
    write_csv(args.out_dir / "su2_quaia_template_projection_correlations.csv", control_rows)
    write_csv(args.out_dir / "su2_quaia_spatial_block_observed_dipoles.csv", observed_rows)
    (args.out_dir / "su2_quaia_spatial_block_collinearity_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")
    report_path, readout_path = write_reports(args.out_dir, prepared, collinearity_rows, spatial_rows, config)

    output_files = [
        args.out_dir / "su2_quaia_spatial_block_summary.csv",
        args.out_dir / "su2_quaia_template_collinearity_summary.csv",
        args.out_dir / "su2_quaia_template_projection_correlations.csv",
        args.out_dir / "su2_quaia_spatial_block_observed_dipoles.csv",
        args.out_dir / "su2_quaia_spatial_block_collinearity_config.json",
        report_path,
        readout_path,
    ]
    manifest_rows = [
        {"file": str(path), "bytes": path.stat().st_size, "sha256": sha256(path), "role": "result"}
        for path in output_files
    ]
    manifest_path = args.out_dir / "su2_quaia_spatial_block_collinearity_manifest.csv"
    write_csv(manifest_path, manifest_rows)
    output_files.append(manifest_path)
    if args.copy_to_outputs:
        OUTPUTS.mkdir(parents=True, exist_ok=True)
        for path in output_files:
            shutil.copy2(path, OUTPUTS / f"{path.stem}_2026-07-16{path.suffix}")
    print(f"Saved report: {report_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
