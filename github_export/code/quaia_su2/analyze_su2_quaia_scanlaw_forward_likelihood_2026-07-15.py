from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits
from scipy.stats import chi2

ROOT = Path(r"C:\Users\clive\Documents\Cosmology")
CODE_DIR = ROOT / "github_export" / "code" / "quaia_su2"
OUTPUTS = Path(r"C:\Users\clive\Documents\Codex\2026-07-14\ok-continu\outputs")
DEFAULT_OUT = ROOT / "plamb_runs" / "diagnostics" / "su2_quaia_scanlaw_forward_likelihood_20260715"
EXT_SCRIPT = "analyze_su2_quaia_external_dust_gate_2026-07-15.py"


def load_module(name: str, candidates: list[Path]) -> Any:
    for path in candidates:
        if path.exists():
            spec = importlib.util.spec_from_file_location(name, path)
            if spec is None or spec.loader is None:
                continue
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            return module
    raise FileNotFoundError(f"Could not find module source for {name}")


EXT = load_module("su2_external_dust_gate", [OUTPUTS / EXT_SCRIPT, CODE_DIR / EXT_SCRIPT])


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


def fmt(value: Any) -> str:
    if isinstance(value, float):
        return f"{value:.6g}" if math.isfinite(value) else ""
    return str(value)


def markdown_table(rows: list[dict[str, Any]], columns: list[str]) -> str:
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(fmt(row.get(col, "")) for col in columns) + " |")
    return "\n".join(lines)


def finite_mask(*arrays: np.ndarray) -> np.ndarray:
    mask = np.ones(len(arrays[0]), dtype=bool)
    for arr in arrays:
        mask &= np.isfinite(arr)
    return mask


def standardize(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    med = float(np.nanmedian(values))
    scale = float(np.nanstd(values))
    if not math.isfinite(scale) or scale <= 0.0:
        scale = 1.0
    return (values - med) / scale


def unit_from_lb(l_deg: float, b_deg: float) -> np.ndarray:
    l = math.radians(l_deg)
    b = math.radians(b_deg)
    cb = math.cos(b)
    return np.array([cb * math.cos(l), cb * math.sin(l), math.sin(b)], dtype=float)


def angular_sep_deg(a: np.ndarray, b: np.ndarray) -> float:
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denom <= 0.0:
        return float("nan")
    dot = float(np.clip(np.dot(a, b) / denom, -1.0, 1.0))
    return math.degrees(math.acos(dot))


def load_data(args: argparse.Namespace) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    gaiascanlaw_data_dir = args.gaiascanlaw_data_dir or EXT.resolve_gaiascanlaw_data_dir()
    q = EXT.load_data(args.quaia, args.selection, args.randoms, args.sfd_dir, gaiascanlaw_data_dir)
    with fits.open(args.quaia, memmap=True) as hdul:
        data = hdul[1].data
        q["bp_rp"] = np.asarray(data["phot_bp_mean_mag"], dtype=float) - np.asarray(data["phot_rp_mean_mag"], dtype=float)
    q["w1_w2"] = q["w1"] - q["w2"]
    base_mask = (
        (q["z"] >= args.z_min)
        & (q["z"] < args.z_max)
        & (np.abs(q["b"]) >= min(args.b_cuts))
        & finite_mask(q["bp_rp"], q["w1_w2"])
    )
    bp_med = float(np.nanmedian(q["bp_rp"][base_mask]))
    bp_sd = float(np.nanstd(q["bp_rp"][base_mask]))
    if not math.isfinite(bp_sd) or bp_sd <= 0.0:
        bp_sd = 1.0
    w_med = float(np.nanmedian(q["w1_w2"][base_mask]))
    w_sd = float(np.nanstd(q["w1_w2"][base_mask]))
    if not math.isfinite(w_sd) or w_sd <= 0.0:
        w_sd = 1.0
    q["bp_rp_z"] = (q["bp_rp"] - bp_med) / bp_sd
    q["w1_w2_z"] = (q["w1_w2"] - w_med) / w_sd
    cross = q["bp_rp_z"] * q["w1_w2_z"]
    c_med = float(np.nanmedian(cross[base_mask]))
    c_sd = float(np.nanstd(cross[base_mask]))
    if not math.isfinite(c_sd) or c_sd <= 0.0:
        c_sd = 1.0
    q["bp_rp_x_w1_w2_z"] = (cross - c_med) / c_sd
    meta = {
        "gaiascanlaw_data_dir": str(gaiascanlaw_data_dir) if gaiascanlaw_data_dir else None,
        "scanlaw_templates_present": bool(np.all(np.isfinite(q["gaia_scan_count_log1p_dr3"][base_mask]))),
    }
    return q, meta


def named_controls(q: dict[str, np.ndarray]) -> dict[str, list[tuple[str, np.ndarray]]]:
    scan = [
        ("gaia_scan_count_log1p_dr3", q["gaia_scan_count_log1p_dr3"]),
        ("gaia_scan_angle_cos2_mean_dr3", q["gaia_scan_angle_cos2_mean_dr3"]),
        ("gaia_scan_angle_sin2_mean_dr3", q["gaia_scan_angle_sin2_mean_dr3"]),
        ("gaia_scan_angle_resultant_dr3", q["gaia_scan_angle_resultant_dr3"]),
    ]
    colour = [
        ("bp_rp_z", q["bp_rp_z"]),
        ("w1_w2_z", q["w1_w2_z"]),
        ("bp_rp_x_w1_w2_z", q["bp_rp_x_w1_w2_z"]),
    ]
    quality = [
        ("zerr", q["zerr"]),
        ("g", q["g"]),
        ("w1", q["w1"]),
        ("w2", q["w2"]),
        ("pmra_error", q["pmra_error"]),
        ("pmdec_error", q["pmdec_error"]),
    ]
    dust_selection = [
        ("ebv_sfd", q["ebv_sfd"]),
        ("sfd_log1p", q["sfd_log1p"]),
        ("sfd_sq", q["sfd_sq"]),
        ("selection_T", q["selection_T"]),
        ("random_density_log1p", q["random_density_log1p"]),
    ]
    return {
        "null_intercept": [],
        "dipole_only": [],
        "scanlaw_colour": scan + colour,
        "scanlaw_colour_plus_dipole": scan + colour,
        "scanlaw_colour_quality": scan + colour + quality,
        "scanlaw_colour_quality_plus_dipole": scan + colour + quality,
        "all_external_colour": scan + colour + dust_selection,
        "all_external_colour_plus_dipole": scan + colour + dust_selection,
    }


def fit_model(z: np.ndarray, vec: np.ndarray, controls: list[tuple[str, np.ndarray]], include_dipole: bool) -> dict[str, Any]:
    columns = [np.ones(len(z))]
    names = ["intercept"]
    dipole_slice: slice | None = None
    if include_dipole:
        start = len(columns)
        columns.extend([vec[:, 0], vec[:, 1], vec[:, 2]])
        names.extend(["n_x", "n_y", "n_z"])
        dipole_slice = slice(start, start + 3)
    for name, values in controls:
        columns.append(standardize(values))
        names.append(name)
    x = np.column_stack(columns)
    beta, *_ = np.linalg.lstsq(x, z, rcond=None)
    fitted = x @ beta
    resid = z - fitted
    rss = float(np.dot(resid, resid))
    tss = float(np.dot(z - np.mean(z), z - np.mean(z)))
    n = int(len(z))
    k = int(x.shape[1])
    sigma2_mle = max(rss / max(n, 1), 1e-300)
    loglike = -0.5 * n * (math.log(2.0 * math.pi * sigma2_mle) + 1.0)
    row: dict[str, Any] = {
        "N": n,
        "k": k,
        "rss": rss,
        "r2": 1.0 - rss / tss if tss > 0.0 else float("nan"),
        "aic": 2.0 * k - 2.0 * loglike,
        "bic": math.log(max(n, 2)) * k - 2.0 * loglike,
        "loglike": loglike,
        "rms_resid": math.sqrt(rss / max(n - k, 1)),
        "columns": ";".join(names),
    }
    if dipole_slice is not None:
        dvec = beta[dipole_slice]
        amp = float(np.linalg.norm(dvec))
        l_deg, b_deg = EXT.unit_to_lb(dvec)
        dof = max(n - k, 1)
        sigma2 = rss / dof
        try:
            cov = sigma2 * np.linalg.pinv(x.T @ x)
            cov_d = cov[dipole_slice, dipole_slice]
            snr = float(math.sqrt(max(float(dvec @ np.linalg.pinv(cov_d) @ dvec), 0.0)))
        except Exception:
            snr = float("nan")
        row.update({"dipole_amp": amp, "dipole_l_deg": l_deg, "dipole_b_deg": b_deg, "dipole_snr": snr})
    else:
        row.update({"dipole_amp": float("nan"), "dipole_l_deg": float("nan"), "dipole_b_deg": float("nan"), "dipole_snr": float("nan")})
    return row


def compare_nested(restricted: dict[str, Any], full: dict[str, Any]) -> dict[str, float]:
    df = int(full["k"] - restricted["k"])
    lrt = float(max(2.0 * (full["loglike"] - restricted["loglike"]), 0.0))
    partial_r2 = 1.0 - float(full["rss"]) / float(restricted["rss"]) if restricted["rss"] > 0.0 else float("nan")
    return {
        "delta_k": df,
        "lrt_stat": lrt,
        "lrt_p": float(chi2.sf(lrt, df)) if df > 0 else float("nan"),
        "partial_r2": partial_r2,
        "delta_aic_full_minus_restricted": float(full["aic"] - restricted["aic"]),
        "delta_bic_full_minus_restricted": float(full["bic"] - restricted["bic"]),
    }


def run(args: argparse.Namespace) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    q, meta = load_data(args)
    control_groups = named_controls(q)
    fit_rows: list[dict[str, Any]] = []
    comparison_rows: list[dict[str, Any]] = []
    for bcut in args.b_cuts:
        base_mask = (
            (q["z"] >= args.z_min)
            & (q["z"] < args.z_max)
            & (np.abs(q["b"]) >= bcut)
            & finite_mask(q["z"], q["vec"][:, 0], q["vec"][:, 1], q["vec"][:, 2], q["bp_rp_z"], q["w1_w2_z"], q["bp_rp_x_w1_w2_z"])
        )
        idx0 = np.flatnonzero(base_mask)
        controls_all = []
        for controls in control_groups.values():
            controls_all.extend(values for _, values in controls)
        finite = np.ones(len(idx0), dtype=bool)
        for values in controls_all:
            finite &= np.isfinite(values[idx0])
        idx = idx0[finite]
        z = q["z"][idx]
        vec = q["vec"][idx]
        model_results: dict[str, dict[str, Any]] = {}
        for model_name, controls_full in control_groups.items():
            include_dipole = model_name in {
                "dipole_only",
                "scanlaw_colour_plus_dipole",
                "scanlaw_colour_quality_plus_dipole",
                "all_external_colour_plus_dipole",
            }
            controls = [(name, values[idx]) for name, values in controls_full]
            result = fit_model(z, vec, controls, include_dipole=include_dipole)
            result.update({"model": model_name, "bcut_deg": bcut})
            fit_rows.append(result)
            model_results[model_name] = result

        comparisons = [
            ("dipole_only_vs_null", "null_intercept", "dipole_only"),
            ("dipole_added_to_scanlaw_colour", "scanlaw_colour", "scanlaw_colour_plus_dipole"),
            ("dipole_added_to_scanlaw_colour_quality", "scanlaw_colour_quality", "scanlaw_colour_quality_plus_dipole"),
            ("dipole_added_to_all_external_colour", "all_external_colour", "all_external_colour_plus_dipole"),
        ]
        raw_axis = unit_from_lb(model_results["dipole_only"]["dipole_l_deg"], model_results["dipole_only"]["dipole_b_deg"])
        for label, restricted_name, full_name in comparisons:
            restricted = model_results[restricted_name]
            full = model_results[full_name]
            comp = compare_nested(restricted, full)
            full_axis = unit_from_lb(full["dipole_l_deg"], full["dipole_b_deg"]) if math.isfinite(full["dipole_l_deg"]) else np.array([np.nan, np.nan, np.nan])
            comp.update(
                {
                    "comparison": label,
                    "bcut_deg": bcut,
                    "N": full["N"],
                    "restricted_model": restricted_name,
                    "full_model": full_name,
                    "restricted_bic": restricted["bic"],
                    "full_bic": full["bic"],
                    "full_dipole_amp": full["dipole_amp"],
                    "full_dipole_snr": full["dipole_snr"],
                    "raw_dipole_amp": model_results["dipole_only"]["dipole_amp"],
                    "raw_dipole_snr": model_results["dipole_only"]["dipole_snr"],
                    "amp_ratio_vs_raw_dipole": full["dipole_amp"] / model_results["dipole_only"]["dipole_amp"] if model_results["dipole_only"]["dipole_amp"] > 0.0 else float("nan"),
                    "direction_sep_vs_raw_dipole_deg": angular_sep_deg(full_axis, raw_axis),
                }
            )
            comparison_rows.append(comp)
        print(f"completed bcut {bcut:g}: N={len(idx)}", flush=True)

    summary_rows: list[dict[str, Any]] = []
    for label in sorted(set(row["comparison"] for row in comparison_rows)):
        rows = [row for row in comparison_rows if row["comparison"] == label]
        summary_rows.append(
            {
                "comparison": label,
                "n_bcuts": len(rows),
                "N_min": int(min(row["N"] for row in rows)),
                "N_max": int(max(row["N"] for row in rows)),
                "mean_delta_bic_full_minus_restricted": float(np.mean([row["delta_bic_full_minus_restricted"] for row in rows])),
                "max_delta_bic_full_minus_restricted": float(np.max([row["delta_bic_full_minus_restricted"] for row in rows])),
                "min_delta_bic_full_minus_restricted": float(np.min([row["delta_bic_full_minus_restricted"] for row in rows])),
                "mean_partial_r2": float(np.mean([row["partial_r2"] for row in rows])),
                "max_full_dipole_snr": float(np.max([row["full_dipole_snr"] for row in rows])),
                "mean_amp_ratio_vs_raw_dipole": float(np.mean([row["amp_ratio_vs_raw_dipole"] for row in rows])),
                "max_direction_sep_vs_raw_dipole_deg": float(np.max([row["direction_sep_vs_raw_dipole_deg"] for row in rows])),
                "bic_readout": "full model favoured" if float(np.mean([row["delta_bic_full_minus_restricted"] for row in rows])) < 0.0 else "restricted model favoured",
            }
        )
    config = {
        "date": "2026-07-15",
        "analysis": "scanlaw_forward_template_likelihood",
        "z_min": args.z_min,
        "z_max": args.z_max,
        "b_cuts": args.b_cuts,
        "external_inputs": meta,
        "models": {
            "dipole_only": "z = alpha + d.n + epsilon",
            "scanlaw_colour": "z = alpha + Gaia_scanlaw beta_s + colour beta_c + epsilon",
            "scanlaw_colour_plus_dipole": "z = alpha + d.n + Gaia_scanlaw beta_s + colour beta_c + epsilon",
        },
        "gate": "The dipole term should add material BIC-supported explanatory power after scan-law and colour controls before the angular mode can be promoted.",
    }
    return fit_rows, comparison_rows, summary_rows, config


def write_reports(out_dir: Path, comparison_rows: list[dict[str, Any]], summary_rows: list[dict[str, Any]], config: dict[str, Any]) -> None:
    summary_cols = [
        "comparison",
        "mean_delta_bic_full_minus_restricted",
        "mean_partial_r2",
        "max_full_dipole_snr",
        "mean_amp_ratio_vs_raw_dipole",
        "max_direction_sep_vs_raw_dipole_deg",
        "bic_readout",
    ]
    detail_cols = [
        "comparison",
        "bcut_deg",
        "N",
        "delta_bic_full_minus_restricted",
        "partial_r2",
        "lrt_p",
        "full_dipole_snr",
        "amp_ratio_vs_raw_dipole",
        "direction_sep_vs_raw_dipole_deg",
    ]
    primary = [row for row in summary_rows if row["comparison"] == "dipole_added_to_scanlaw_colour"]
    if primary:
        p = primary[0]
        if p["mean_delta_bic_full_minus_restricted"] < -10.0:
            bottom = "The dipole term is BIC-supported after scan-law plus colour controls."
        elif p["mean_delta_bic_full_minus_restricted"] < 0.0:
            bottom = "The dipole term gives weak-to-moderate BIC improvement after scan-law plus colour controls."
        else:
            bottom = "Scan-law plus colour controls are favoured over adding a dipole; the angular mode does not pass this mechanism gate."
    else:
        bottom = "Primary scan-law plus colour comparison did not run."
    report = [
        "# SU2 / Quaia Gaia Scan-Law Forward-Template Likelihood",
        "",
        "Date: July 15, 2026",
        "",
        "## Purpose",
        "",
        "This test treats Gaia scan-law and catalogue-colour fields as an explicit forward-template model, then asks whether a redshift dipole adds explanatory power after those fields are already in the likelihood.",
        "",
        "## Model Set",
        "",
        "```text",
        "M0: z_i = alpha + epsilon_i",
        "Md: z_i = alpha + d . n_i + epsilon_i",
        "Mt: z_i = alpha + S_i beta_s + C_i beta_c + epsilon_i",
        "Mtd: z_i = alpha + d . n_i + S_i beta_s + C_i beta_c + epsilon_i",
        "```",
        "",
        "`S_i` contains Gaia scan count and scan-angle geometry. `C_i` contains Gaia BP-RP, WISE W1-W2, and their colour cross-term.",
        "",
        "## Bottom Line",
        "",
        bottom,
        "",
        "## Summary",
        "",
        markdown_table(summary_rows, summary_cols),
        "",
        "## B-Cut Details",
        "",
        markdown_table(comparison_rows, detail_cols),
        "",
        "## Configuration",
        "",
        "```json",
        json.dumps(config, indent=2),
        "```",
    ]
    (out_dir / "su2_quaia_scanlaw_forward_likelihood_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    readout = [
        "# SU2 / Quaia Gaia Scan-Law Forward-Template Likelihood Readout",
        "",
        "Date: July 15, 2026",
        "",
        "## Bottom Line",
        "",
        bottom,
        "",
        "## Summary",
        "",
        markdown_table(summary_rows, summary_cols),
    ]
    (out_dir / "su2_quaia_scanlaw_forward_likelihood_readout.md").write_text("\n".join(readout) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Gaia scan-law forward-template likelihood for the locked SU2/Quaia mode.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--quaia", type=Path, default=EXT.DEFAULT_QUAIA)
    parser.add_argument("--randoms", type=Path, default=EXT.DEFAULT_RANDOMS)
    parser.add_argument("--selection", type=Path, default=EXT.DEFAULT_SELECTION)
    parser.add_argument("--sfd-dir", type=Path, default=EXT.DEFAULT_SFD_DIR)
    parser.add_argument("--gaiascanlaw-data-dir", type=Path, default=None)
    parser.add_argument("--z-min", type=float, default=0.95)
    parser.add_argument("--z-max", type=float, default=1.45)
    parser.add_argument("--b-cuts", type=float, nargs="+", default=[15.0, 20.0, 25.0, 30.0, 35.0])
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    fit_rows, comparison_rows, summary_rows, config = run(args)
    write_csv(args.out_dir / "su2_quaia_scanlaw_forward_likelihood_fit_rows.csv", fit_rows)
    write_csv(args.out_dir / "su2_quaia_scanlaw_forward_likelihood_comparisons.csv", comparison_rows)
    write_csv(args.out_dir / "su2_quaia_scanlaw_forward_likelihood_summary.csv", summary_rows)
    (args.out_dir / "su2_quaia_scanlaw_forward_likelihood_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")
    write_reports(args.out_dir, comparison_rows, summary_rows, config)
    print(f"Saved report: {args.out_dir / 'su2_quaia_scanlaw_forward_likelihood_report.md'}", flush=True)
    print(f"Saved readout: {args.out_dir / 'su2_quaia_scanlaw_forward_likelihood_readout.md'}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
