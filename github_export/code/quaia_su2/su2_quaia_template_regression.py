from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

try:
    from astropy.io import fits
except Exception as exc:  # pragma: no cover
    fits = None
    FITS_IMPORT_ERROR = exc
else:
    FITS_IMPORT_ERROR = None


ROOT = Path(r"C:\Users\clive\Documents\Cosmology")
DEFAULT_QUAIA_FITS = ROOT / "quaia_G20.0.fits"
DEFAULT_OUT = ROOT / "plamb_runs" / "diagnostics" / "su2_quaia_template_regression_20260714"


def gal_unit_array(l_deg: np.ndarray, b_deg: np.ndarray) -> np.ndarray:
    l = np.deg2rad(l_deg)
    b = np.deg2rad(b_deg)
    cb = np.cos(b)
    return np.column_stack([cb * np.cos(l), cb * np.sin(l), np.sin(b)])


def unit_to_lb(vec: np.ndarray) -> tuple[float, float]:
    norm = float(np.linalg.norm(vec))
    if norm <= 0.0 or not math.isfinite(norm):
        return float("nan"), float("nan")
    unit = vec / norm
    return math.degrees(math.atan2(float(unit[1]), float(unit[0]))) % 360.0, math.degrees(math.asin(float(np.clip(unit[2], -1.0, 1.0))))


def standardize(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    med = np.nanmedian(values)
    scale = np.nanstd(values)
    if not math.isfinite(scale) or scale <= 0.0:
        scale = 1.0
    return (values - med) / scale


def fit_model(z: np.ndarray, vec: np.ndarray, extras: list[np.ndarray]) -> dict[str, float]:
    parts = [np.ones(len(z)), vec[:, 0], vec[:, 1], vec[:, 2]]
    parts.extend(extras)
    x = np.column_stack(parts)
    beta, *_ = np.linalg.lstsq(x, z, rcond=None)
    fitted = x @ beta
    resid = z - fitted
    dof = max(len(z) - x.shape[1], 1)
    sigma2 = float(np.dot(resid, resid) / dof)
    dvec = beta[1:4]
    amp = float(np.linalg.norm(dvec))
    l_deg, b_deg = unit_to_lb(dvec)
    xtx_inv = np.linalg.pinv(x.T @ x)
    cov_d = sigma2 * xtx_inv[1:4, 1:4]
    try:
        snr = float(math.sqrt(max(float(dvec @ np.linalg.pinv(cov_d) @ dvec), 0.0)))
    except Exception:
        snr = float("nan")
    return {
        "N": int(len(z)),
        "n_params": int(x.shape[1]),
        "amp": amp,
        "l_deg": l_deg,
        "b_deg": b_deg,
        "formal_snr": snr,
        "rms_resid": float(math.sqrt(sigma2)),
    }


def build_quadrupole(vec: np.ndarray) -> list[np.ndarray]:
    nx, ny, nz = vec[:, 0], vec[:, 1], vec[:, 2]
    return [
        nx * nx - 1.0 / 3.0,
        ny * ny - 1.0 / 3.0,
        nz * nz - 1.0 / 3.0,
        nx * ny,
        nx * nz,
        ny * nz,
    ]


def finite_mask(*arrays: np.ndarray) -> np.ndarray:
    mask = np.ones(len(arrays[0]), dtype=bool)
    for arr in arrays:
        mask &= np.isfinite(arr)
    return mask


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


def markdown_table(rows: list[dict[str, Any]], columns: list[str]) -> str:
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in rows:
        vals = []
        for col in columns:
            value = row.get(col, "")
            if isinstance(value, float):
                vals.append(f"{value:.6g}" if math.isfinite(value) else "")
            else:
                vals.append(str(value))
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)


def load_quaia(path: Path) -> dict[str, np.ndarray]:
    if fits is None:
        raise RuntimeError(f"astropy.io.fits unavailable: {FITS_IMPORT_ERROR}")
    with fits.open(path, memmap=True) as hdul:
        data = hdul[1].data
        out = {
            "z": np.asarray(data["redshift_quaia"], dtype=float),
            "zerr": np.asarray(data["redshift_quaia_err"], dtype=float),
            "l": np.asarray(data["l"], dtype=float),
            "b": np.asarray(data["b"], dtype=float),
            "g": np.asarray(data["phot_g_mean_mag"], dtype=float),
            "bp": np.asarray(data["phot_bp_mean_mag"], dtype=float),
            "rp": np.asarray(data["phot_rp_mean_mag"], dtype=float),
            "w1": np.asarray(data["mag_w1_vg"], dtype=float),
            "w2": np.asarray(data["mag_w2_vg"], dtype=float),
            "pm": np.asarray(data["pm"], dtype=float),
            "pmra": np.asarray(data["pmra"], dtype=float),
            "pmdec": np.asarray(data["pmdec"], dtype=float),
        }
    return out


def run(args: argparse.Namespace) -> list[dict[str, Any]]:
    q = load_quaia(args.quaia_fits)
    z = q["z"]
    vec_all = gal_unit_array(q["l"], q["b"])
    rows: list[dict[str, Any]] = []
    windows = [(1.0, 1.5, "z1p00_1p50")]
    bcuts = args.bcuts

    photometry = [
        standardize(q["zerr"]),
        standardize(q["g"]),
        standardize(q["bp"] - q["rp"]),
        standardize(q["w1"]),
        standardize(q["w2"]),
        standardize(q["w1"] - q["w2"]),
        standardize(q["pm"]),
        standardize(q["pmra"]),
        standardize(q["pmdec"]),
    ]

    for z_lo, z_hi, tag in windows:
        for bcut in bcuts:
            base_mask = (z >= z_lo) & (z < z_hi) & (np.abs(q["b"]) >= bcut) & np.isfinite(q["l"]) & np.isfinite(q["b"]) & np.isfinite(z)
            idx_all = np.flatnonzero(base_mask)
            if len(idx_all) < args.min_count:
                continue

            z_all = z[idx_all]
            vec = vec_all[idx_all]
            quad = [standardize(v[idx_all]) for v in build_quadrupole(vec_all)]
            phot = [p[idx_all] for p in photometry]
            common = finite_mask(z_all, vec[:, 0], vec[:, 1], vec[:, 2], *phot, *quad)

            model_specs = [
                ("baseline_all", np.ones(len(z_all), dtype=bool), []),
                ("baseline_common", common, []),
                ("photometry_common", common, phot),
                ("quadrupole_common", common, quad),
                ("photometry_plus_quadrupole_common", common, phot + quad),
            ]
            baseline_common_amp = float("nan")
            baseline_common_snr = float("nan")
            for model_name, mask, extras in model_specs:
                if np.sum(mask) < args.min_count:
                    continue
                fit = fit_model(z_all[mask], vec[mask], [e[mask] for e in extras])
                if model_name == "baseline_common":
                    baseline_common_amp = fit["amp"]
                    baseline_common_snr = fit["formal_snr"]
                amp_ratio = fit["amp"] / baseline_common_amp if math.isfinite(baseline_common_amp) and baseline_common_amp > 0 else float("nan")
                snr_delta = fit["formal_snr"] - baseline_common_snr if math.isfinite(baseline_common_snr) else float("nan")
                rows.append(
                    {
                        "tag": tag,
                        "z_lo": z_lo,
                        "z_hi": z_hi,
                        "bcut_deg": bcut,
                        "model": model_name,
                        "N": fit["N"],
                        "n_params": fit["n_params"],
                        "amp": fit["amp"],
                        "amp_ratio_vs_baseline_common": amp_ratio,
                        "formal_snr": fit["formal_snr"],
                        "snr_delta_vs_baseline_common": snr_delta,
                        "l_deg": fit["l_deg"],
                        "b_deg": fit["b_deg"],
                        "rms_resid": fit["rms_resid"],
                    }
                )
            print(f"{tag} |b|>{bcut:g}: N_all={len(idx_all)} N_common={int(np.sum(common))}", flush=True)
    return rows


def write_report(path: Path, rows: list[dict[str, Any]], args: argparse.Namespace) -> None:
    lines = [
        "# SU2 / Quaia Template Regression Control",
        "",
        "Date: July 14, 2026",
        "",
        "## Purpose",
        "",
        "This is a first-pass selection/systematics stress test for the surviving `z1p00_1p50` Quaia amplitude hook. It refits the redshift dipole after adding internal catalog-quality proxies and smooth Galactic quadrupole terms.",
        "",
        "This is not yet a full external dust, scanning-law, or survey-depth likelihood. It is an immediate internal-control regression using object-level columns available in `quaia_G20.0.fits`.",
        "",
        "## Controls",
        "",
        "- photometry/common controls: redshift error, Gaia G, BP-RP, W1, W2, W1-W2, total proper motion, pmra, pmdec.",
        "- quadrupole/common controls: six smooth Galactic quadrupole basis terms.",
        "- common-sample rows compare against `baseline_common`, so amplitude changes are not just missing-data artifacts.",
        "",
        "## Summary",
        "",
        markdown_table(
            rows,
            [
                "tag",
                "bcut_deg",
                "model",
                "N",
                "amp",
                "amp_ratio_vs_baseline_common",
                "formal_snr",
                "snr_delta_vs_baseline_common",
                "l_deg",
                "b_deg",
            ],
        ),
        "",
        "## Interpretation Gate",
        "",
        "- A large dipole-amplitude drop under photometry or quadrupole controls points toward selection/depth leakage.",
        "- A stable amplitude under these internal controls still does not promote the result physically; it only justifies moving to external dust/scanning-law/depth templates.",
        "- The global look-elsewhere mock audit remains the first promotion gate.",
        "",
        "## Outputs",
        "",
        f"- summary CSV: `{args.out / 'su2_quaia_template_regression_summary.csv'}`",
        f"- report: `{args.out / 'su2_quaia_template_regression_report.md'}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Internal template-regression control for SU2/Quaia z1p00_1p50 targets.")
    parser.add_argument("--quaia-fits", type=Path, default=DEFAULT_QUAIA_FITS)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--bcuts", type=float, nargs="+", default=[10, 15, 20, 25, 30, 35])
    parser.add_argument("--min-count", type=int, default=1000)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    rows = run(args)
    write_csv(args.out / "su2_quaia_template_regression_summary.csv", rows)
    write_report(args.out / "su2_quaia_template_regression_report.md", rows, args)
    print(f"Saved summary: {args.out / 'su2_quaia_template_regression_summary.csv'}", flush=True)
    print(f"Saved report: {args.out / 'su2_quaia_template_regression_report.md'}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
