from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import astropy.units as u
import requests
from astropy.io import fits
from astropy.wcs import WCS
from astropy_healpix import HEALPix
from scipy.spatial import cKDTree

import analyze_su2_quaia_external_dust_gate as ext


ROOT = Path(r"C:\Users\clive\Documents\Cosmology")
DEFAULT_WISE = ROOT / "external_datasets" / "allwise_atlas_inventory" / "allwise_p3am_cdd_w1w2_coverage.csv"
DEFAULT_GAIA_DENSITY = ROOT / "external_datasets" / "gaia_dr3_density" / "gaia_dr3_source_density_hpx4_glt20p5.csv"
DEFAULT_GAIA_HIPS_CACHE = ROOT / "external_datasets" / "gaia_edr3_density_hips"
DEFAULT_GAIA_HIPS_BASE = "https://alasky.cds.unistra.fr/ancillary/GaiaEDR3/density-map"
DEFAULT_OUT = ROOT / "plamb_runs" / "diagnostics" / "su2_quaia_wise_stellar_gate_20260715"


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


def eq_unit_array(ra_deg: np.ndarray, dec_deg: np.ndarray) -> np.ndarray:
    ra = np.deg2rad(ra_deg)
    dec = np.deg2rad(dec_deg)
    cd = np.cos(dec)
    return np.column_stack([cd * np.cos(ra), cd * np.sin(ra), np.sin(dec)])


def standardize(values: np.ndarray) -> np.ndarray:
    med = np.nanmedian(values)
    scale = np.nanstd(values)
    if not math.isfinite(scale) or scale <= 0.0:
        scale = 1.0
    return (values - med) / scale


def load_wise_depth(path: Path, q: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    df = pd.read_csv(path)
    out: dict[str, np.ndarray] = {}
    q_xyz = eq_unit_array(q["ra"], q["dec"])
    for band in (1, 2):
        sub = df[df["band"] == band].copy()
        xyz = eq_unit_array(sub["crval1"].to_numpy(float), sub["crval2"].to_numpy(float))
        tree = cKDTree(xyz)
        dist, idx = tree.query(q_xyz, k=1)
        nearest = sub.iloc[idx].reset_index(drop=True)
        prefix = f"wise_w{band}"
        out[f"{prefix}_tile_sep_deg"] = np.rad2deg(2.0 * np.arcsin(np.clip(dist / 2.0, 0.0, 1.0)))
        for col in ["medcov", "mincov", "maxcov", "lowcovpc1", "lowcovpc2", "nomcovpc", "mincovpc"]:
            out[f"{prefix}_{col}"] = nearest[col].to_numpy(float)
        out[f"{prefix}_medcov_log1p"] = np.log1p(np.clip(out[f"{prefix}_medcov"], 0.0, None))
        out[f"{prefix}_lowcov_log1p"] = np.log1p(np.clip(out[f"{prefix}_lowcovpc1"] + out[f"{prefix}_lowcovpc2"], 0.0, None))
    out["wise_w1_w2_medcov_diff"] = out["wise_w1_medcov"] - out["wise_w2_medcov"]
    out["wise_w1_w2_medcov_ratio"] = out["wise_w1_medcov"] / np.clip(out["wise_w2_medcov"], 1e-6, None)
    return out


def fetch_hips_tile(base_url: str, cache_dir: Path, order: int, tile: int) -> Path:
    dir_id = (tile // 10000) * 10000
    out_dir = cache_dir / f"Norder{order}" / f"Dir{dir_id}"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"Npix{tile}.fits"
    if out_path.exists() and out_path.stat().st_size > 1000:
        return out_path
    url = f"{base_url.rstrip('/')}/Norder{order}/Dir{dir_id}/Npix{tile}.fits"
    response = requests.get(url, timeout=90)
    response.raise_for_status()
    out_path.write_bytes(response.content)
    return out_path


def sample_gaia_density_hips(
    q: dict[str, np.ndarray],
    base_url: str,
    cache_dir: Path,
    order: int,
) -> np.ndarray:
    hp = HEALPix(nside=2**order, order="nested")
    hpx = hp.lonlat_to_healpix(q["ra"] * u.deg, q["dec"] * u.deg)
    density = np.full(len(q["z"]), np.nan, dtype=float)
    unique_tiles = np.unique(hpx)
    for i, tile in enumerate(unique_tiles, start=1):
        tile = int(tile)
        tile_path = fetch_hips_tile(base_url, cache_dir, order, tile)
        with fits.open(tile_path, memmap=True) as hdul:
            data = np.asarray(hdul[0].data, dtype=float)
            wcs = WCS(hdul[0].header)
            sub = hpx == tile
            x, y = wcs.wcs_world2pix(q["ra"][sub], q["dec"][sub], 0)
            ix = np.clip(np.rint(x).astype(int), 0, data.shape[1] - 1)
            iy = np.clip(np.rint(y).astype(int), 0, data.shape[0] - 1)
            density[sub] = data[iy, ix]
        if i == 1 or i % 25 == 0 or i == len(unique_tiles):
            print(f"Sampled Gaia EDR3 density HiPS tile {i}/{len(unique_tiles)}", flush=True)
    return density


def load_gaia_density(
    path: Path,
    q: dict[str, np.ndarray],
    hips_base_url: str,
    hips_cache_dir: Path,
    hips_order: int,
) -> tuple[dict[str, np.ndarray], str]:
    if not path.exists():
        density = sample_gaia_density_hips(q, hips_base_url, hips_cache_dir, hips_order)
        return (
            {
                "gaia_dr3_density_hpx4": density,
                "gaia_dr3_density_hpx4_log1p": np.log1p(np.clip(density, 0.0, None)),
            },
            f"CDS Gaia EDR3 density HiPS Norder{hips_order}",
        )
    df = pd.read_csv(path)
    hp = HEALPix(nside=16, order="nested")
    hpx = hp.lonlat_to_healpix(q["ra"] * u.deg, q["dec"] * u.deg)
    values = np.full(12 * 16 * 16, np.nan, dtype=float)
    values[df["hpx"].to_numpy(int)] = df["n_source"].to_numpy(float)
    density = values[hpx]
    return (
        {
            "gaia_dr3_density_hpx4": density,
            "gaia_dr3_density_hpx4_log1p": np.log1p(np.clip(density, 0.0, None)),
        },
        "Gaia DR3 TAP source-density HEALPix level 4",
    )


def summarize(control_group: str, fits: list[dict[str, Any]], baseline: dict[str, Any] | None) -> dict[str, Any]:
    row = ext.summarize(control_group, fits, baseline)
    return row


def corr_row(q: dict[str, np.ndarray], mask: np.ndarray, name: str, values: np.ndarray) -> dict[str, Any]:
    return ext.corr_row(q, mask, name, values)


def run(args: argparse.Namespace) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    args.out_dir.mkdir(parents=True, exist_ok=True)
    gaiascanlaw_data_dir = args.gaiascanlaw_data_dir or ext.resolve_gaiascanlaw_data_dir()
    q = ext.load_data(ext.DEFAULT_QUAIA, ext.DEFAULT_SELECTION, ext.DEFAULT_RANDOMS, ext.DEFAULT_SFD_DIR, gaiascanlaw_data_dir)
    q.update(load_wise_depth(args.wise_inventory, q))
    gaia_density_templates, gaia_density_source = load_gaia_density(
        args.gaia_density,
        q,
        args.gaia_density_hips_base_url,
        args.gaia_density_hips_cache_dir,
        args.gaia_density_hips_order,
    )
    q.update(gaia_density_templates)

    wise_controls = [
        ("wise_w1_medcov_log1p", q["wise_w1_medcov_log1p"]),
        ("wise_w2_medcov_log1p", q["wise_w2_medcov_log1p"]),
        ("wise_w1_mincov", q["wise_w1_mincov"]),
        ("wise_w2_mincov", q["wise_w2_mincov"]),
        ("wise_w1_lowcov_log1p", q["wise_w1_lowcov_log1p"]),
        ("wise_w2_lowcov_log1p", q["wise_w2_lowcov_log1p"]),
        ("wise_w1_w2_medcov_diff", q["wise_w1_w2_medcov_diff"]),
        ("wise_w1_w2_medcov_ratio", q["wise_w1_w2_medcov_ratio"]),
    ]
    stellar_controls = [
        ("gaia_dr3_density_hpx4_log1p", q["gaia_dr3_density_hpx4_log1p"]),
    ]
    base_external = [
        ("ebv_sfd", q["ebv_sfd"]),
        ("sfd_log1p", q["sfd_log1p"]),
        ("sfd_sq", q["sfd_sq"]),
        ("selection_T", q["selection_T"]),
        ("random_density_log1p", q["random_density_log1p"]),
        ("gaia_scan_count_log1p_dr3", q["gaia_scan_count_log1p_dr3"]),
        ("gaia_scan_angle_cos2_mean_dr3", q["gaia_scan_angle_cos2_mean_dr3"]),
        ("gaia_scan_angle_sin2_mean_dr3", q["gaia_scan_angle_sin2_mean_dr3"]),
        ("gaia_scan_angle_resultant_dr3", q["gaia_scan_angle_resultant_dr3"]),
    ]
    local_quality = [
        ("zerr", q["zerr"]),
        ("g", q["g"]),
        ("w1", q["w1"]),
        ("w2", q["w2"]),
        ("w1_w2", q["w1_w2"]),
        ("pmra_error", q["pmra_error"]),
        ("pmdec_error", q["pmdec_error"]),
    ]

    groups: dict[str, list[tuple[str, np.ndarray]]] = {
        "baseline_no_template": [],
        "wise_depth_external": wise_controls,
        "gaia_stellar_density_external": stellar_controls,
        "wise_plus_stellar_external": wise_controls + stellar_controls,
        "sfd_selection_scanlaw_external": base_external,
        "all_external_wise_stellar": base_external + wise_controls + stellar_controls,
        "local_catalog_quality_proxy": local_quality,
        "all_external_plus_local_quality": base_external + wise_controls + stellar_controls + local_quality,
    }
    if not np.all(np.isfinite(q["gaia_dr3_density_hpx4_log1p"])):
        groups.pop("gaia_stellar_density_external")
        groups.pop("wise_plus_stellar_external")
        groups.pop("all_external_wise_stellar")
        groups.pop("all_external_plus_local_quality")

    fit_rows: list[dict[str, Any]] = []
    for group_name, group_controls in groups.items():
        for b_cut in args.b_cuts:
            mask = (
                (q["z"] >= args.z_min)
                & (q["z"] < args.z_max)
                & (np.abs(q["b"]) >= b_cut)
                & ext.finite_mask(q["z"], q["vec"][:, 0], q["vec"][:, 1], q["vec"][:, 2])
            )
            controls = [values[mask] for _, values in group_controls]
            if controls and not all(np.all(np.isfinite(control)) for control in controls):
                raise RuntimeError(f"Non-finite controls in {group_name}")
            result = ext.fit_dipole(q["z"][mask], q["vec"][mask], controls=controls)
            fit_rows.append(
                {
                    "control_group": group_name,
                    "b_cut_deg": b_cut,
                    "N": int(np.sum(mask)),
                    "n_templates": len(group_controls),
                    "template_names": ";".join(name for name, _ in group_controls),
                    **result,
                }
            )

    summary_rows: list[dict[str, Any]] = []
    baseline: dict[str, Any] | None = None
    for group_name in groups:
        rows = [row for row in fit_rows if row["control_group"] == group_name]
        row = summarize(group_name, rows, baseline)
        summary_rows.append(row)
        if group_name == "baseline_no_template":
            baseline = row
            summary_rows[-1] = summarize(group_name, rows, None)

    primary_mask = (q["z"] >= args.z_min) & (q["z"] < args.z_max) & (np.abs(q["b"]) >= 15.0)
    corr_names = [
        "wise_w1_medcov_log1p",
        "wise_w2_medcov_log1p",
        "wise_w1_mincov",
        "wise_w2_mincov",
        "wise_w1_lowcov_log1p",
        "wise_w2_lowcov_log1p",
        "wise_w1_w2_medcov_diff",
        "gaia_dr3_density_hpx4_log1p",
        "ebv_sfd",
        "selection_T",
        "random_density_log1p",
        "gaia_scan_count_log1p_dr3",
        "gaia_scan_angle_resultant_dr3",
        "w1_w2",
    ]
    corr_rows = [corr_row(q, primary_mask, name, q[name]) for name in corr_names if name in q]

    config = {
        "date": "2026-07-15",
        "z_min": args.z_min,
        "z_max": args.z_max,
        "b_cuts": args.b_cuts,
        "wise_inventory": str(args.wise_inventory),
        "gaia_density": str(args.gaia_density),
        "gaia_density_present": bool(np.all(np.isfinite(q["gaia_dr3_density_hpx4_log1p"]))),
        "gaia_density_source": gaia_density_source,
        "gaia_density_hips_base_url": args.gaia_density_hips_base_url,
        "gaia_density_hips_order": args.gaia_density_hips_order,
        "method": "AllWISE Atlas W1/W2 tile coverage matched by nearest tile center; Gaia DR3 density sampled from TAP HEALPix level 4 when available, otherwise from CDS Gaia EDR3 density HiPS tiles.",
    }
    return fit_rows, summary_rows, corr_rows, config


def write_reports(out_dir: Path, summary_rows: list[dict[str, Any]], corr_rows: list[dict[str, Any]], config: dict[str, Any]) -> None:
    cols = [
        "control_group",
        "amp_ratio_vs_baseline",
        "snr_max",
        "snr_delta_vs_baseline",
        "max_pair_direction_sep_deg",
        "mean_resultant_length",
        "N_min",
        "N_max",
    ]
    corr_cols = ["template", "N", "mean", "std", "corr_z"]
    by = {row["control_group"]: row for row in summary_rows}
    wise = by.get("wise_depth_external")
    stellar = by.get("gaia_stellar_density_external")
    all_external = by.get("all_external_wise_stellar")
    if wise and wise["amp_ratio_vs_baseline"] >= 0.9 and wise["max_pair_direction_sep_deg"] <= 25.0:
        wise_line = "AllWISE W1/W2 tile-depth controls do not absorb the locked angular mode."
    else:
        wise_line = "AllWISE W1/W2 tile-depth controls materially suppress or rotate the locked angular mode."
    if stellar is None:
        stellar_line = "Gaia stellar-density control was not included because the external density map was not available at run time."
    elif stellar["amp_ratio_vs_baseline"] >= 0.9 and stellar["max_pair_direction_sep_deg"] <= 25.0:
        stellar_line = "Gaia DR3 stellar-density control does not absorb the locked angular mode."
    else:
        stellar_line = "Gaia DR3 stellar-density control materially suppresses or rotates the locked angular mode."
    if all_external is None:
        combined_line = "The full WISE+stellar external stack could not be evaluated because Gaia density was unavailable."
    elif all_external["amp_ratio_vs_baseline"] >= 0.75:
        combined_line = "The full external stack including WISE depth and stellar density retains a substantial residual mode."
    else:
        combined_line = "The full external stack including WISE depth and stellar density removes most of the baseline mode."

    report = [
        "# SU2 / Quaia WISE-Depth and Stellar-Density Gate",
        "",
        "Date: July 15, 2026",
        "",
        "## Bottom Line",
        "",
        f"- {wise_line}",
        f"- {stellar_line}",
        f"- {combined_line}",
        "",
        "## Control Summary",
        "",
        markdown_table(summary_rows, cols),
        "",
        "## Template Correlations In Primary Mask",
        "",
        markdown_table(corr_rows, corr_cols),
        "",
        "## External Inputs",
        "",
        "- AllWISE Atlas Inventory W1/W2 tile-level depth-of-coverage metadata (`medcov`, `mincov`, `maxcov`, low-coverage fractions), matched to Quaia objects by nearest Atlas tile center.",
        "- Gaia DR3 source-density map from `gaiadr3.gaia_source_lite`, `phot_g_mean_mag < 20.5`, HEALPix level 4, when the TAP aggregate is available; otherwise CDS Gaia EDR3 density HiPS tiles are sampled.",
        "- Existing SFD, Quaia selection/random-density, and Gaia scan-law templates are reused from the prior gate.",
        "",
        "## Configuration",
        "",
        "```json",
        json.dumps(config, indent=2),
        "```",
        "",
    ]
    (out_dir / "su2_quaia_wise_stellar_gate_report.md").write_text("\n".join(report), encoding="utf-8")

    readout = [
        "# SU2 / Quaia WISE-Depth and Stellar-Density Gate Readout",
        "",
        "Date: July 15, 2026",
        "",
        "## Bottom Line",
        "",
        wise_line,
        "",
        stellar_line,
        "",
        combined_line,
        "",
        "## Key Numbers",
        "",
        markdown_table(summary_rows, cols),
        "",
    ]
    (out_dir / "su2_quaia_wise_stellar_gate_readout.md").write_text("\n".join(readout), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run SU2/Quaia WISE-depth and Gaia stellar-density external gate.")
    parser.add_argument("--wise-inventory", type=Path, default=DEFAULT_WISE)
    parser.add_argument("--gaia-density", type=Path, default=DEFAULT_GAIA_DENSITY)
    parser.add_argument("--gaia-density-hips-base-url", default=DEFAULT_GAIA_HIPS_BASE)
    parser.add_argument("--gaia-density-hips-cache-dir", type=Path, default=DEFAULT_GAIA_HIPS_CACHE)
    parser.add_argument("--gaia-density-hips-order", type=int, default=2)
    parser.add_argument("--gaiascanlaw-data-dir", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--z-min", type=float, default=0.95)
    parser.add_argument("--z-max", type=float, default=1.45)
    parser.add_argument("--b-cuts", type=float, nargs="+", default=[15.0, 20.0, 25.0, 30.0, 35.0])
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    fit_rows, summary_rows, corr_rows, config = run(args)
    write_csv(args.out_dir / "su2_quaia_wise_stellar_gate_fit_rows.csv", fit_rows)
    write_csv(args.out_dir / "su2_quaia_wise_stellar_gate_summary.csv", summary_rows)
    write_csv(args.out_dir / "su2_quaia_wise_stellar_gate_template_correlations.csv", corr_rows)
    (args.out_dir / "su2_quaia_wise_stellar_gate_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")
    write_reports(args.out_dir, summary_rows, corr_rows, config)
    print(f"Saved readout: {args.out_dir / 'su2_quaia_wise_stellar_gate_readout.md'}")
    print(f"Saved report: {args.out_dir / 'su2_quaia_wise_stellar_gate_report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
