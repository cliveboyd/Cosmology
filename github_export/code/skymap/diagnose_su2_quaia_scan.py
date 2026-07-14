#!/usr/bin/env python3
r"""
Quaia angular/redshift diagnostic scan for the SU2/SU2R programme.

Quaia is not a scalar distance ladder like SN/BAO. This script treats it as an
independent angular quasar-redshift probe:

    * Load the local compact Quaia RA/DEC/Z table.
    * Convert ICRS RA/DEC to Galactic coordinates.
    * Scan Galactic latitude cuts and redshift windows.
    * Fit a simple dipole: z = a0 + D . n.
    * Compare the redshift windows with LCDM/SU2/SU2R scalar distance residuals.

The result is a diagnostic/targeting scan, not a formal SU2 likelihood. A
publication-grade use would need the Quaia selection function, mask, redshift
systematics, and mock catalogs in the likelihood.

Outputs:
    plamb_runs/diagnostics/su2_quaia_scan/su2_quaia_scan.csv
    plamb_runs/diagnostics/su2_quaia_scan/su2_quaia_scan_top.csv
    plamb_runs/diagnostics/su2_quaia_scan/su2_quaia_hht_series.csv
    plamb_runs/diagnostics/su2_quaia_scan/su2_quaia_scan_report.md
    plamb_runs/diagnostics/su2_quaia_scan/su2_quaia_scan_overview.png

Run:
    python diagnose_su2_quaia_scan.py
"""

from __future__ import annotations

import argparse
import ast
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as exc:  # pragma: no cover - optional plotting
    plt = None
    MATPLOTLIB_IMPORT_ERROR = exc
else:
    MATPLOTLIB_IMPORT_ERROR = None

from PLamb_Test_10V6c_plus import mu_model


ROOT = Path(__file__).resolve().parent
DEFAULT_QUAIA = ROOT / "quaia_G20p0_basic.csv"
OUTDIR = ROOT / "plamb_runs" / "diagnostics" / "su2_quaia_scan"

# Standard J2000/ICRS equatorial to Galactic rotation matrix.
# Rows are Galactic x/y/z axes in equatorial Cartesian coordinates.
EQ_TO_GAL = np.array(
    [
        [-0.0548755604, -0.8734370902, -0.4838350155],
        [0.4941094279, -0.4448296300, 0.7469822445],
        [-0.8676661490, -0.1980763734, 0.4559837762],
    ],
    dtype=float,
)

CMB_L_DEG = 264.021
CMB_B_DEG = 48.253


@dataclass
class QuaiaData:
    ra: np.ndarray
    dec: np.ndarray
    z: np.ndarray
    gal_l: np.ndarray
    gal_b: np.ndarray
    gal_vec: np.ndarray


@dataclass
class BestFit:
    label: str
    model: str
    path: Path
    values: dict[str, Any]


def gal_unit_from_lb(l_deg: float, b_deg: float) -> np.ndarray:
    l = math.radians(l_deg)
    b = math.radians(b_deg)
    cb = math.cos(b)
    return np.array([cb * math.cos(l), cb * math.sin(l), math.sin(b)], dtype=float)


def unit_to_lb(vec: np.ndarray) -> tuple[float, float]:
    norm = float(np.linalg.norm(vec))
    if norm <= 0.0 or not math.isfinite(norm):
        return float("nan"), float("nan")
    unit = vec / norm
    l = math.degrees(math.atan2(unit[1], unit[0])) % 360.0
    b = math.degrees(math.asin(float(np.clip(unit[2], -1.0, 1.0))))
    return l, b


def angular_sep_deg(a: np.ndarray, b: np.ndarray) -> float:
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denom <= 0.0:
        return float("nan")
    dot = float(np.clip(np.dot(a, b) / denom, -1.0, 1.0))
    return math.degrees(math.acos(dot))


def load_quaia(path: Path) -> QuaiaData:
    rows = np.genfromtxt(path, delimiter=",", names=True, dtype=float, encoding="utf-8")
    ra = np.asarray(rows["RA"], dtype=float)
    dec = np.asarray(rows["DEC"], dtype=float)
    z = np.asarray(rows["Z"], dtype=float)
    ok = np.isfinite(ra) & np.isfinite(dec) & np.isfinite(z) & (z > 0.0)
    ra = ra[ok]
    dec = dec[ok]
    z = z[ok]

    ra_rad = np.deg2rad(ra)
    dec_rad = np.deg2rad(dec)
    eq = np.column_stack(
        [
            np.cos(dec_rad) * np.cos(ra_rad),
            np.cos(dec_rad) * np.sin(ra_rad),
            np.sin(dec_rad),
        ]
    )
    gal = eq @ EQ_TO_GAL.T
    gal /= np.linalg.norm(gal, axis=1)[:, None]
    gal_l = (np.rad2deg(np.arctan2(gal[:, 1], gal[:, 0])) + 360.0) % 360.0
    gal_b = np.rad2deg(np.arcsin(np.clip(gal[:, 2], -1.0, 1.0)))
    return QuaiaData(ra=ra, dec=dec, z=z, gal_l=gal_l, gal_b=gal_b, gal_vec=gal)


def fit_redshift_dipole(z: np.ndarray, vec: np.ndarray) -> dict[str, float]:
    x = np.column_stack([np.ones(len(z)), vec])
    beta, *_ = np.linalg.lstsq(x, z, rcond=None)
    yfit = x @ beta
    resid = z - yfit
    dof = max(len(z) - 4, 1)
    sigma2 = float(np.sum(resid * resid) / dof)
    xtx_inv = np.linalg.pinv(x.T @ x)
    cov_d = sigma2 * xtx_inv[1:, 1:]
    dvec = np.asarray(beta[1:4], dtype=float)
    amp = float(np.linalg.norm(dvec))
    l_deg, b_deg = unit_to_lb(dvec)
    try:
        snr = float(math.sqrt(max(dvec @ np.linalg.pinv(cov_d) @ dvec, 0.0)))
    except Exception:
        snr = float("nan")
    cmb = gal_unit_from_lb(CMB_L_DEG, CMB_B_DEG)
    sep = angular_sep_deg(dvec, cmb)
    amp_par = float(np.dot(dvec, cmb))
    frac_par = amp_par / amp if amp > 0.0 else float("nan")
    return {
        "a0": float(beta[0]),
        "amp": amp,
        "l_deg": l_deg,
        "b_deg": b_deg,
        "sep_cmb_deg": sep,
        "amp_par_cmb": amp_par,
        "frac_par_cmb": frac_par,
        "formal_snr": snr,
        "rms_resid": float(math.sqrt(sigma2)),
    }


def read_bestfit(path: Path) -> dict[str, Any]:
    values: dict[str, Any] = {}
    if not path.exists():
        return values
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if "=" not in line:
                continue
            key, value = line.split("=", 1)
            key = key.strip()
            value = value.strip()
            try:
                values[key] = ast.literal_eval(value)
            except Exception:
                try:
                    values[key] = float(value)
                except ValueError:
                    values[key] = value
    return values


def load_bestfits() -> list[BestFit]:
    candidates = [
        (
            "LCDM_DESI_DR2",
            "LCDM",
            ROOT / "plamb_runs" / "updated_datasets_20260710" / "runs" / "desi_LCDM_Pantheon_SN_BAO_Planck" / "LCDM_emcee_bestfit.txt",
        ),
        (
            "SU2_DESI_DR2_fixed_posterior",
            "SU2",
            ROOT
            / "plamb_runs"
            / "updated_datasets_20260710"
            / "runs"
            / "desi_SU2_fixed_posterior_Pantheon_SN_BAO_Planck"
            / "SU2_emcee_bestfit.txt",
        ),
        (
            "SU2_DESI_DR2_legacy",
            "SU2",
            ROOT / "plamb_runs" / "updated_datasets_20260710" / "runs" / "desi_SU2_Pantheon_SN_BAO_Planck" / "SU2_emcee_bestfit.txt",
        ),
        (
            "SU2R_DESI_DR2",
            "SU2R",
            ROOT / "plamb_runs" / "updated_datasets_20260710" / "runs" / "desi_SU2R_Pantheon_SN_BAO_Planck" / "SU2R_emcee_bestfit.txt",
        ),
    ]
    bestfits: list[BestFit] = []
    for label, model, path in candidates:
        values = read_bestfit(path)
        if values:
            bestfits.append(BestFit(label=label, model=model, path=path, values=values))
    return bestfits


def model_mu_at(z: float, bestfit: BestFit, nint: int) -> float:
    v = bestfit.values
    return float(
        mu_model(
            np.array([z], dtype=float),
            float(v.get("H0", 70.0)),
            float(v.get("Om", 0.3)),
            float(v.get("Ol", 0.7)),
            0.0,
            0.0,
            0.0,
            0.0,
            nint_base=nint,
            Omega_chi0=float(v.get("Omega_chi0", 0.0)),
            w0_chi=float(v.get("w0_chi", -1.0)),
            wa_chi=float(v.get("wa_chi", 0.0)),
        )[0]
    )


def scalar_residual_context(z: float, bestfits: list[BestFit], nint: int) -> dict[str, float]:
    lcdm = next((bf for bf in bestfits if bf.label == "LCDM_DESI_DR2"), None)
    if lcdm is None:
        return {}
    mu_lcdm = model_mu_at(z, lcdm, nint)
    values: dict[str, float] = {"mu_LCDM_DESI_DR2": mu_lcdm}
    for bf in bestfits:
        if bf is lcdm:
            continue
        key = f"delta_mu_{bf.label}_minus_LCDM"
        values[key] = model_mu_at(z, bf, nint) - mu_lcdm
    return values


def make_bins(mode: str, z_min: float, z_max: float, dz: float) -> list[tuple[str, float | None, float | None]]:
    bins: list[tuple[str, float | None, float | None]] = []
    if mode in {"both", "broad"}:
        bins.extend(
            [
                ("full", None, None),
                ("z0p10_0p30", 0.10, 0.30),
                ("z0p30_0p50", 0.30, 0.50),
                ("z0p50_0p75", 0.50, 0.75),
                ("z0p75_1p00", 0.75, 1.00),
                ("z1p00_1p50", 1.00, 1.50),
                ("z1p50_2p50", 1.50, 2.50),
            ]
        )
    if mode in {"both", "narrow"}:
        lo = z_min
        while lo < z_max - 1e-12:
            hi = min(lo + dz, z_max)
            tag = f"z{lo:.2f}_{hi:.2f}".replace(".", "p")
            bins.append((tag, lo, hi))
            lo = hi
    return bins


def scan_quaia(data: QuaiaData, args: argparse.Namespace, bestfits: list[BestFit]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    bins = make_bins(args.bin_mode, args.z_min, args.z_max, args.dz)
    for bcut in args.bcuts:
        base_mask = np.abs(data.gal_b) >= bcut
        for tag, zlo, zhi in bins:
            mask = base_mask.copy()
            if zlo is not None:
                mask &= data.z >= zlo
            if zhi is not None:
                mask &= data.z < zhi
            n = int(np.sum(mask))
            if n < args.min_count:
                continue
            z = data.z[mask]
            vec = data.gal_vec[mask]
            fit = fit_redshift_dipole(z, vec)
            z_center = float(np.median(z)) if tag == "full" else float(0.5 * ((zlo or 0.0) + (zhi or 0.0)))
            row: dict[str, Any] = {
                "bcut_deg": bcut,
                "tag": tag,
                "z_lo": "" if zlo is None else zlo,
                "z_hi": "" if zhi is None else zhi,
                "z_center": z_center,
                "N": n,
                "z_mean": float(np.mean(z)),
                "z_median": float(np.median(z)),
            }
            row.update(fit)
            row.update(scalar_residual_context(z_center, bestfits, args.nint))
            deltas = [abs(value) for key, value in row.items() if key.startswith("delta_mu_SU2") and isinstance(value, float)]
            max_su2_delta = max(deltas) if deltas else 0.0
            row["max_abs_su2_delta_mu"] = max_su2_delta
            row["su2_quaia_priority_score"] = float(fit["formal_snr"] * max_su2_delta)
            rows.append(row)
    rows.sort(key=lambda row: (float(row["bcut_deg"]), str(row["tag"])))
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_hht_series(path: Path, rows: list[dict[str, Any]], bcut: float) -> None:
    selected = [
        row
        for row in rows
        if float(row["bcut_deg"]) == float(bcut)
        and str(row["tag"]).startswith("z")
        and float(row["N"]) >= 50
    ]
    selected.sort(key=lambda row: float(row["z_center"]))
    out_rows = [
        {
            "x": row["z_center"],
            "amp_par_cmb": row["amp_par_cmb"],
            "frac_par_cmb": row["frac_par_cmb"],
            "dipole_amp": row["amp"],
            "formal_snr": row["formal_snr"],
            "N": row["N"],
            "tag": row["tag"],
        }
        for row in selected
    ]
    write_csv(path, out_rows)


def write_plot(path: Path, rows: list[dict[str, Any]], bcut: float) -> str:
    if plt is None:
        return f"Plot skipped because matplotlib is unavailable: {MATPLOTLIB_IMPORT_ERROR}"
    selected = [
        row
        for row in rows
        if float(row["bcut_deg"]) == float(bcut)
        and str(row["tag"]).startswith("z")
        and float(row["N"]) >= 50
    ]
    selected.sort(key=lambda row: float(row["z_center"]))
    if not selected:
        return "Plot skipped: no selected rows."
    z = np.array([float(row["z_center"]) for row in selected])
    amp = np.array([float(row["amp"]) for row in selected])
    apar = np.array([float(row["amp_par_cmb"]) for row in selected])
    snr = np.array([float(row["formal_snr"]) for row in selected])
    delta_key = next((key for key in selected[0].keys() if key.startswith("delta_mu_SU2R")), None)
    delta = np.array([float(row.get(delta_key, 0.0) or 0.0) for row in selected]) if delta_key else np.zeros_like(z)
    scale = np.nanmax(np.abs(apar)) / max(np.nanmax(np.abs(delta)), 1e-12)

    fig, axes = plt.subplots(2, 1, figsize=(9.0, 7.0), dpi=150, sharex=True)
    axes[0].plot(z, amp, "o-", lw=1.5, label="Quaia redshift dipole amplitude")
    axes[0].plot(z, np.abs(apar), "s-", lw=1.3, label="|CMB-projected component|")
    axes[0].set_ylabel("Dipole amplitude in z")
    axes[0].grid(True, alpha=0.25)
    axes[0].legend(fontsize=8)
    axes[1].plot(z, apar, "o-", lw=1.5, label="CMB-projected dipole component")
    if delta_key:
        axes[1].plot(z, delta * scale, "--", lw=1.4, label=f"{delta_key} scaled")
    axes[1].set_xlabel("Redshift window center")
    axes[1].set_ylabel("Projected component / scaled SU2 residual")
    axes[1].grid(True, alpha=0.25)
    axes[1].legend(fontsize=8)
    ax2 = axes[1].twinx()
    ax2.plot(z, snr, ":", color="0.35", lw=1.2, label="formal dipole SNR")
    ax2.set_ylabel("Formal SNR")
    fig.suptitle(f"SU2-Quaia diagnostic scan, |b| >= {bcut:g} deg")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)
    return f"Saved plot: {path}"


def markdown_table(rows: list[dict[str, Any]], columns: list[str]) -> list[str]:
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in rows:
        vals: list[str] = []
        for col in columns:
            value = row.get(col, "")
            if isinstance(value, float):
                vals.append(f"{value:.6g}" if math.isfinite(value) else "")
            else:
                vals.append(str(value))
        lines.append("| " + " | ".join(vals) + " |")
    return lines


def write_report(
    path: Path,
    data: QuaiaData,
    rows: list[dict[str, Any]],
    top_rows: list[dict[str, Any]],
    bestfits: list[BestFit],
    args: argparse.Namespace,
    plot_status: str,
) -> None:
    full_rows = [row for row in rows if row["tag"] == "full"]
    lines = [
        "# SU2 / Quaia Angular Diagnostic Scan",
        "",
        "## Headline",
        "",
        f"- Loaded local Quaia compact table: `{args.quaia_csv}`.",
        f"- Quasar rows used: `{len(data.z)}`.",
        f"- Scan rows written: `{len(rows)}`.",
        "- Fit per window: `z = a0 + D.n` using object-level Galactic unit vectors.",
        "- This is a diagnostic scan, not a formal SU2 likelihood.",
        "",
        "## Scalar SU2 Context Used",
        "",
    ]
    if bestfits:
        lines.extend(markdown_table([{"label": bf.label, "model": bf.model, "path": str(bf.path)} for bf in bestfits], ["label", "model", "path"]))
    else:
        lines.append("- No scalar best-fit files were available for residual context.")
    lines.extend(
        [
            "",
            "## Full-Sample Latitude-Cut Dipoles",
            "",
        ]
    )
    lines.extend(markdown_table(full_rows, ["bcut_deg", "N", "amp", "l_deg", "b_deg", "sep_cmb_deg", "amp_par_cmb", "frac_par_cmb", "formal_snr", "max_abs_su2_delta_mu"]))
    lines.extend(
        [
            "",
            "## Highest SU2-Quaia Priority Windows",
            "",
            "The priority score is `formal_dipole_snr * max_abs_SU2_delta_mu`. It is a targeting score only: it flags redshift windows where angular Quaia structure overlaps the scalar redshift range where SU2/SU2R differs most from LCDM.",
            "",
        ]
    )
    lines.extend(markdown_table(top_rows, ["bcut_deg", "tag", "z_center", "N", "amp", "l_deg", "b_deg", "sep_cmb_deg", "formal_snr", "max_abs_su2_delta_mu", "su2_quaia_priority_score"]))
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- Quaia can test angular anisotropy, dipole/quadrupole structure, and possible parity/chirality signatures missing from scalar SN/BAO/Planck fits.",
            "- SU2/SU2R scalar distance fits have no sky direction, so this script does not claim a direct likelihood preference for SU2.",
            "- This current scan is object-level. Older Quaia summaries in `.spyder-py3/skymap2` used HEALPix/pixel-quality products, so `N` and dipole amplitudes are not expected to match one-for-one.",
            "- Rows with high formal SNR should be checked against Quaia selection, Galactic mask, survey footprint, and mock catalogs before being treated as physics.",
            "- The most useful next step is to rerun this scan with Quaia random mocks and export the CMB-projected series into `diagnose_hht_resonance.py --mode csv`.",
            "",
            "## Local Outputs",
            "",
            f"- Full scan: `{OUTDIR / 'su2_quaia_scan.csv'}`",
            f"- Top windows: `{OUTDIR / 'su2_quaia_scan_top.csv'}`",
            f"- HHT series: `{OUTDIR / 'su2_quaia_hht_series.csv'}`",
            f"- Plot: `{OUTDIR / 'su2_quaia_scan_overview.png'}`",
            f"- Plot status: {plot_status}",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quaia-csv", type=Path, default=DEFAULT_QUAIA)
    parser.add_argument("--out", type=Path, default=OUTDIR)
    parser.add_argument("--bcuts", type=float, nargs="+", default=[10, 15, 20, 25, 30, 35])
    parser.add_argument("--bin-mode", choices=["broad", "narrow", "both"], default="both")
    parser.add_argument("--z-min", type=float, default=0.30)
    parser.add_argument("--z-max", type=float, default=2.50)
    parser.add_argument("--dz", type=float, default=0.10)
    parser.add_argument("--min-count", type=int, default=50)
    parser.add_argument("--plot-bcut", type=float, default=25.0)
    parser.add_argument("--nint", type=int, default=180)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    if not args.quaia_csv.exists():
        print(f"Missing Quaia CSV: {args.quaia_csv}", file=sys.stderr)
        return 2

    data = load_quaia(args.quaia_csv)
    bestfits = load_bestfits()
    rows = scan_quaia(data, args, bestfits)
    if not rows:
        print("No scan rows passed the filters.", file=sys.stderr)
        return 2

    top_rows = sorted(
        [row for row in rows if row.get("tag") != "full"],
        key=lambda row: float(row.get("su2_quaia_priority_score", 0.0)),
        reverse=True,
    )[:25]

    scan_csv = args.out / "su2_quaia_scan.csv"
    top_csv = args.out / "su2_quaia_scan_top.csv"
    hht_csv = args.out / "su2_quaia_hht_series.csv"
    report_md = args.out / "su2_quaia_scan_report.md"
    plot_png = args.out / "su2_quaia_scan_overview.png"

    write_csv(scan_csv, rows)
    write_csv(top_csv, top_rows)
    write_hht_series(hht_csv, rows, args.plot_bcut)
    plot_status = write_plot(plot_png, rows, args.plot_bcut)
    write_report(report_md, data, rows, top_rows, bestfits, args, plot_status)

    best = top_rows[0]
    print(f"Loaded Quaia rows: {len(data.z)}")
    print(f"Saved scan: {scan_csv}")
    print(f"Saved top windows: {top_csv}")
    print(f"Saved HHT series: {hht_csv}")
    print(f"Saved report: {report_md}")
    print(plot_status)
    print(
        "Top SU2-Quaia window: "
        f"bcut={best['bcut_deg']} tag={best['tag']} z={best['z_center']:.3f} "
        f"N={best['N']} amp={best['amp']:.4g} formal_snr={best['formal_snr']:.3g} "
        f"score={best['su2_quaia_priority_score']:.4g}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
