from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:
    from astropy.io import fits
except Exception as exc:  # pragma: no cover
    fits = None
    FITS_IMPORT_ERROR = exc
else:
    FITS_IMPORT_ERROR = None

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as exc:  # pragma: no cover
    plt = None
    MATPLOTLIB_IMPORT_ERROR = exc
else:
    MATPLOTLIB_IMPORT_ERROR = None


ROOT = Path(r"C:\Users\clive\Documents\Cosmology")
DEFAULT_QUAIA = ROOT / "quaia_G20p0_basic_gal.csv"
DEFAULT_RANDOMS = ROOT / "quaia_G20p0_randoms_simple.fits"
DEFAULT_SCAN = ROOT / "plamb_runs" / "diagnostics" / "su2_quaia_scan" / "su2_quaia_scan.csv"
DEFAULT_TOP = ROOT / "plamb_runs" / "diagnostics" / "su2_quaia_scan" / "su2_quaia_scan_top.csv"
DEFAULT_OUT = ROOT / "plamb_runs" / "diagnostics" / "quaia_redshift_window_full_mocks_20260714"

CMB_L_DEG = 264.021
CMB_B_DEG = 48.253

EQ_TO_GAL = np.array(
    [
        [-0.0548755604, -0.8734370902, -0.4838350155],
        [0.4941094279, -0.4448296300, 0.7469822445],
        [-0.8676661490, -0.1980763734, 0.4559837762],
    ],
    dtype=float,
)


@dataclass(frozen=True)
class Window:
    bcut_deg: float
    tag: str
    z_lo: float
    z_hi: float
    z_center: float
    observed_N: int
    observed_amp: float
    observed_l_deg: float
    observed_b_deg: float
    observed_sep_cmb_deg: float
    observed_amp_par_cmb: float
    observed_frac_par_cmb: float
    observed_formal_snr: float
    priority_score: float


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
    l_deg = math.degrees(math.atan2(float(unit[1]), float(unit[0]))) % 360.0
    b_deg = math.degrees(math.asin(float(np.clip(unit[2], -1.0, 1.0))))
    return l_deg, b_deg


def angular_sep_deg(a: np.ndarray, b: np.ndarray) -> float:
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denom <= 0.0:
        return float("nan")
    dot = float(np.clip(np.dot(a, b) / denom, -1.0, 1.0))
    return math.degrees(math.acos(dot))


def equatorial_to_galactic(ra_deg: np.ndarray, dec_deg: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    ra = np.deg2rad(ra_deg)
    dec = np.deg2rad(dec_deg)
    eq = np.column_stack(
        [
            np.cos(dec) * np.cos(ra),
            np.cos(dec) * np.sin(ra),
            np.sin(dec),
        ]
    )
    gal = eq @ EQ_TO_GAL.T
    gal /= np.linalg.norm(gal, axis=1)[:, None]
    l_gal = (np.rad2deg(np.arctan2(gal[:, 1], gal[:, 0])) + 360.0) % 360.0
    b_gal = np.rad2deg(np.arcsin(np.clip(gal[:, 2], -1.0, 1.0)))
    return l_gal, b_gal, gal


def load_quaia(path: Path) -> dict[str, np.ndarray]:
    df = pd.read_csv(path)
    cols = {c.upper(): c for c in df.columns}
    ra = df[cols["RA"]].to_numpy(float)
    dec = df[cols["DEC"]].to_numpy(float)
    z = df[cols["Z"]].to_numpy(float)
    if "L_GAL" in cols and "B_GAL" in cols:
        l_gal = df[cols["L_GAL"]].to_numpy(float)
        b_gal = df[cols["B_GAL"]].to_numpy(float)
        vec = gal_unit_array(l_gal, b_gal)
    else:
        l_gal, b_gal, vec = equatorial_to_galactic(ra, dec)
    ok = np.isfinite(ra) & np.isfinite(dec) & np.isfinite(z) & np.isfinite(l_gal) & np.isfinite(b_gal) & (z > 0)
    return {
        "ra": ra[ok],
        "dec": dec[ok],
        "z": z[ok],
        "l_gal": l_gal[ok],
        "b_gal": b_gal[ok],
        "vec": vec[ok],
    }


def gal_unit_array(l_deg: np.ndarray, b_deg: np.ndarray) -> np.ndarray:
    l = np.deg2rad(l_deg)
    b = np.deg2rad(b_deg)
    cb = np.cos(b)
    return np.column_stack([cb * np.cos(l), cb * np.sin(l), np.sin(b)])


def load_randoms(path: Path, cache_path: Path | None = None) -> dict[str, np.ndarray]:
    if cache_path is not None and cache_path.exists():
        z = np.load(cache_path)
        return {key: np.asarray(z[key]) for key in ["ra", "dec", "l_gal", "b_gal", "vec"]}
    if fits is None:
        raise RuntimeError(f"astropy.io.fits unavailable: {FITS_IMPORT_ERROR}")
    with fits.open(path, memmap=True) as hdul:
        data = hdul[1].data
        ra = np.asarray(data["RA"], dtype=float)
        dec = np.asarray(data["DEC"], dtype=float)
    l_gal, b_gal, vec = equatorial_to_galactic(ra, dec)
    out = {"ra": ra, "dec": dec, "l_gal": l_gal, "b_gal": b_gal, "vec": vec}
    if cache_path is not None:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(cache_path, **out)
    return out


def parse_float(value: Any, default: float = float("nan")) -> float:
    if value in ("", None):
        return default
    try:
        return float(value)
    except Exception:
        return default


def read_windows(scan_csv: Path, top_csv: Path, mode: str, top_n: int, min_count: int) -> list[Window]:
    path = top_csv if mode == "top" and top_csv.exists() else scan_csv
    df = pd.read_csv(path)
    df = df[df["tag"].astype(str).str.startswith("z")].copy()
    df["N"] = pd.to_numeric(df["N"], errors="coerce")
    df = df[df["N"] >= min_count]
    if mode == "top":
        if "su2_quaia_priority_score" in df.columns:
            df["sort_score"] = pd.to_numeric(df["su2_quaia_priority_score"], errors="coerce").fillna(0.0)
            df = df.sort_values("sort_score", ascending=False)
        df = df.drop_duplicates(["bcut_deg", "tag"]).head(top_n)
    else:
        df = df.sort_values(["bcut_deg", "z_lo", "z_hi", "tag"])

    windows: list[Window] = []
    for _, row in df.iterrows():
        z_lo = parse_float(row.get("z_lo"))
        z_hi = parse_float(row.get("z_hi"))
        if not math.isfinite(z_lo) or not math.isfinite(z_hi):
            continue
        windows.append(
            Window(
                bcut_deg=parse_float(row.get("bcut_deg")),
                tag=str(row.get("tag")),
                z_lo=z_lo,
                z_hi=z_hi,
                z_center=parse_float(row.get("z_center")),
                observed_N=int(parse_float(row.get("N"), 0.0)),
                observed_amp=parse_float(row.get("amp")),
                observed_l_deg=parse_float(row.get("l_deg")),
                observed_b_deg=parse_float(row.get("b_deg")),
                observed_sep_cmb_deg=parse_float(row.get("sep_cmb_deg")),
                observed_amp_par_cmb=parse_float(row.get("amp_par_cmb")),
                observed_frac_par_cmb=parse_float(row.get("frac_par_cmb")),
                observed_formal_snr=parse_float(row.get("formal_snr")),
                priority_score=parse_float(row.get("su2_quaia_priority_score"), 0.0),
            )
        )
    return windows


def fit_redshift_dipole(z: np.ndarray, vec: np.ndarray) -> dict[str, float]:
    n = len(z)
    ones_sum = float(n)
    sx = np.sum(vec, axis=0)
    xtx = np.empty((4, 4), dtype=float)
    xtx[0, 0] = ones_sum
    xtx[0, 1:] = sx
    xtx[1:, 0] = sx
    xtx[1:, 1:] = vec.T @ vec
    xty = np.empty(4, dtype=float)
    xty[0] = float(np.sum(z))
    xty[1:] = vec.T @ z
    beta = np.linalg.pinv(xtx) @ xty
    yty = float(np.dot(z, z))
    sse = max(yty - float(beta @ xty), 0.0)
    dof = max(n - 4, 1)
    sigma2 = sse / dof
    dvec = beta[1:4]
    amp = float(np.linalg.norm(dvec))
    l_deg, b_deg = unit_to_lb(dvec)
    cmb = gal_unit_from_lb(CMB_L_DEG, CMB_B_DEG)
    sep = angular_sep_deg(dvec, cmb)
    amp_par = float(np.dot(dvec, cmb))
    frac_par = amp_par / amp if amp > 0 else float("nan")
    try:
        cov_d = sigma2 * np.linalg.pinv(xtx)[1:, 1:]
        snr = float(math.sqrt(max(float(dvec @ np.linalg.pinv(cov_d) @ dvec), 0.0)))
    except Exception:
        snr = float("nan")
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


def p_ge(values: np.ndarray, observed: float) -> float:
    values = values[np.isfinite(values)]
    if len(values) == 0 or not math.isfinite(observed):
        return float("nan")
    return float((1 + np.sum(values >= observed)) / (1 + len(values)))


def p_le(values: np.ndarray, observed: float) -> float:
    values = values[np.isfinite(values)]
    if len(values) == 0 or not math.isfinite(observed):
        return float("nan")
    return float((1 + np.sum(values <= observed)) / (1 + len(values)))


def q(values: np.ndarray, prob: float) -> float:
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return float("nan")
    return float(np.quantile(values, prob))


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


def save_mock_catalog(path: Path, randoms: dict[str, np.ndarray], rand_idx: np.ndarray, z_mock: np.ndarray, window: Window, mock_id: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        mock_id=np.asarray(mock_id),
        bcut_deg=np.asarray(window.bcut_deg),
        tag=np.asarray(window.tag),
        z_lo=np.asarray(window.z_lo),
        z_hi=np.asarray(window.z_hi),
        ra=randoms["ra"][rand_idx],
        dec=randoms["dec"][rand_idx],
        l_gal=randoms["l_gal"][rand_idx],
        b_gal=randoms["b_gal"][rand_idx],
        z=z_mock,
    )


def run_mocks(args: argparse.Namespace) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    args.out.mkdir(parents=True, exist_ok=True)
    quaia = load_quaia(args.quaia_csv)
    cache_path = args.out / "quaia_randoms_galactic_cache.npz" if args.cache_randoms else None
    randoms = load_randoms(args.randoms_fits, cache_path)
    windows = read_windows(args.scan_csv, args.top_csv, args.window_mode, args.top_n, args.min_count)
    if not windows:
        raise RuntimeError("No redshift windows selected.")

    rng = np.random.default_rng(args.seed)
    summary_rows: list[dict[str, Any]] = []
    pvalue_rows: list[dict[str, Any]] = []
    manifest_rows: list[dict[str, Any]] = []

    for w_i, window in enumerate(windows, start=1):
        obs_mask = (
            (quaia["z"] >= window.z_lo)
            & (quaia["z"] < window.z_hi)
            & (np.abs(quaia["b_gal"]) >= window.bcut_deg)
        )
        obs_idx = np.flatnonzero(obs_mask)
        if len(obs_idx) < args.min_count:
            continue
        obs_z = quaia["z"][obs_idx]
        rand_pool = np.flatnonzero(np.abs(randoms["b_gal"]) >= window.bcut_deg)
        if len(rand_pool) == 0:
            continue
        replace_randoms = len(rand_pool) < len(obs_idx)
        catalog_paths: list[str] = []
        mock_stats: list[dict[str, float]] = []
        for mock_id in range(args.n_mocks):
            rand_idx = rng.choice(rand_pool, size=len(obs_idx), replace=replace_randoms)
            z_mock = rng.permutation(obs_z)
            stats = fit_redshift_dipole(z_mock, randoms["vec"][rand_idx])
            mock_stats.append(stats)
            row = {
                "window_index": w_i,
                "mock_id": mock_id,
                "bcut_deg": window.bcut_deg,
                "tag": window.tag,
                "z_lo": window.z_lo,
                "z_hi": window.z_hi,
                "z_center": window.z_center,
                "N": len(obs_idx),
                **stats,
            }
            summary_rows.append(row)
            if mock_id < args.save_first_catalogs:
                safe_tag = window.tag.replace(".", "p").replace("-", "m")
                catalog_path = args.out / "mock_catalogs" / f"mock_window{w_i:03d}_{safe_tag}_bcut{window.bcut_deg:g}_mock{mock_id:03d}.npz"
                save_mock_catalog(catalog_path, randoms, rand_idx, z_mock, window, mock_id)
                catalog_paths.append(str(catalog_path))
        amp = np.array([row["amp"] for row in mock_stats], dtype=float)
        abs_par = np.abs(np.array([row["amp_par_cmb"] for row in mock_stats], dtype=float))
        snr = np.array([row["formal_snr"] for row in mock_stats], dtype=float)
        sep = np.array([row["sep_cmb_deg"] for row in mock_stats], dtype=float)
        obs_abs_par = abs(window.observed_amp_par_cmb)
        p_row = {
            "window_index": w_i,
            "bcut_deg": window.bcut_deg,
            "tag": window.tag,
            "z_lo": window.z_lo,
            "z_hi": window.z_hi,
            "z_center": window.z_center,
            "N": len(obs_idx),
            "n_mocks": args.n_mocks,
            "observed_amp": window.observed_amp,
            "mock_amp_mean": float(np.nanmean(amp)),
            "mock_amp_p50": q(amp, 0.50),
            "mock_amp_p95": q(amp, 0.95),
            "mock_amp_p99": q(amp, 0.99),
            "p_amp_ge_observed": p_ge(amp, window.observed_amp),
            "observed_abs_amp_par_cmb": obs_abs_par,
            "mock_abs_par_p95": q(abs_par, 0.95),
            "mock_abs_par_p99": q(abs_par, 0.99),
            "p_abs_par_ge_observed": p_ge(abs_par, obs_abs_par),
            "observed_formal_snr": window.observed_formal_snr,
            "mock_snr_p95": q(snr, 0.95),
            "mock_snr_p99": q(snr, 0.99),
            "p_snr_ge_observed": p_ge(snr, window.observed_formal_snr),
            "observed_sep_cmb_deg": window.observed_sep_cmb_deg,
            "mock_sep_p05": q(sep, 0.05),
            "p_sep_le_observed": p_le(sep, window.observed_sep_cmb_deg),
            "priority_score": window.priority_score,
            "saved_catalog_paths": ";".join(catalog_paths),
        }
        pvalue_rows.append(p_row)
        manifest_rows.append(
            {
                "window_index": w_i,
                "bcut_deg": window.bcut_deg,
                "tag": window.tag,
                "z_lo": window.z_lo,
                "z_hi": window.z_hi,
                "N_observed": len(obs_idx),
                "random_pool_N": len(rand_pool),
                "n_mocks": args.n_mocks,
                "mock_catalog_format": "npz with RA, DEC, L_GAL, B_GAL, Z for saved first catalogs",
                "saved_catalog_count": len(catalog_paths),
                "saved_catalog_paths": ";".join(catalog_paths),
            }
        )
        print(
            f"[{w_i}/{len(windows)}] {window.tag} bcut={window.bcut_deg:g} N={len(obs_idx)} "
            f"p_amp={pvalue_rows[-1]['p_amp_ge_observed']:.4g} p_abs_par={pvalue_rows[-1]['p_abs_par_ge_observed']:.4g}",
            flush=True,
        )
    return summary_rows, pvalue_rows, manifest_rows


def markdown_table(rows: list[dict[str, Any]], columns: list[str], max_rows: int = 12) -> str:
    rows = rows[:max_rows]
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


def write_report(path: Path, pvalues: list[dict[str, Any]], manifest: list[dict[str, Any]], args: argparse.Namespace) -> None:
    p_sorted = sorted(pvalues, key=lambda row: float(row.get("p_amp_ge_observed", 1.0)))
    par_sorted = sorted(pvalues, key=lambda row: float(row.get("p_abs_par_ge_observed", 1.0)))
    lines = [
        "# Quaia Redshift-Window Full Mock Catalogs",
        "",
        f"Date: July 14, 2026",
        "",
        "## Purpose",
        "",
        "This run replaces the earlier summary-only Quaia mock control for redshift-window anisotropy. For each selected observed redshift window, it generates object-level mock catalogs by putting the observed window redshift distribution onto random Quaia sky positions passing the same Galactic latitude cut, then refits the same redshift dipole model `z = a0 + D.n`.",
        "",
        "The saved mock catalogs are compact `.npz` object catalogs for the first requested mock(s) in each window. All mocks are nevertheless generated and fitted object-by-object; the per-mock fit summaries are retained in CSV.",
        "",
        "## Configuration",
        "",
        f"- observed Quaia table: `{args.quaia_csv}`",
        f"- random sky catalog: `{args.randoms_fits}`",
        f"- scan CSV: `{args.scan_csv}`",
        f"- top window CSV: `{args.top_csv}`",
        f"- window mode: `{args.window_mode}`",
        f"- top_n: `{args.top_n}`",
        f"- n_mocks per window: `{args.n_mocks}`",
        f"- saved full catalogs per window: `{args.save_first_catalogs}`",
        f"- seed: `{args.seed}`",
        "",
        "## Smallest Amplitude p-values",
        "",
        markdown_table(
            p_sorted,
            [
                "window_index",
                "bcut_deg",
                "tag",
                "N",
                "observed_amp",
                "mock_amp_p95",
                "mock_amp_p99",
                "p_amp_ge_observed",
                "priority_score",
            ],
        ),
        "",
        "## Smallest CMB-Projected p-values",
        "",
        markdown_table(
            par_sorted,
            [
                "window_index",
                "bcut_deg",
                "tag",
                "N",
                "observed_abs_amp_par_cmb",
                "mock_abs_par_p95",
                "mock_abs_par_p99",
                "p_abs_par_ge_observed",
                "priority_score",
            ],
        ),
        "",
        "## Interpretation Gate",
        "",
        "- Do not interpret a Quaia anisotropy window as physical unless it remains unusual against these object-level, window-matched mocks.",
        "- A p-value near the finite mock floor means the window needs more mocks, not immediate promotion.",
        "- These mocks preserve the observed redshift distribution inside each window and the random-catalog sky footprint, but they are still not a complete Quaia selection-function likelihood.",
        "",
        "## Outputs",
        "",
        f"- per-mock summary: `{args.out / 'quaia_redshift_window_mock_summary.csv'}`",
        f"- observed-vs-mock p-values: `{args.out / 'quaia_redshift_window_mock_pvalues.csv'}`",
        f"- manifest: `{args.out / 'quaia_redshift_window_mock_manifest.csv'}`",
        f"- config: `{args.out / 'quaia_redshift_window_mock_config.json'}`",
        f"- saved mock catalogs directory: `{args.out / 'mock_catalogs'}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_plot(path: Path, pvalues: list[dict[str, Any]]) -> str:
    if plt is None:
        return f"Plot skipped because matplotlib is unavailable: {MATPLOTLIB_IMPORT_ERROR}"
    rows = sorted(pvalues, key=lambda row: int(row["window_index"]))
    labels = [f"{row['window_index']}:{row['tag']} |b|>{row['bcut_deg']:g}" for row in rows]
    p_amp = np.array([row["p_amp_ge_observed"] for row in rows], dtype=float)
    p_par = np.array([row["p_abs_par_ge_observed"] for row in rows], dtype=float)
    fig_h = max(5.5, 0.35 * len(rows))
    fig, ax = plt.subplots(figsize=(11.0, fig_h), dpi=150)
    y = np.arange(len(rows))
    ax.scatter(p_amp, y, label="amplitude p", color="#2E74B5")
    ax.scatter(p_par, y, label="|CMB projection| p", color="#8E5A3C")
    ax.axvline(0.05, color="#444444", lw=1.0, ls="--")
    ax.set_xscale("log")
    ax.set_xlabel("mock p-value")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.invert_yaxis()
    ax.grid(True, axis="x", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)
    return f"Saved plot: {path}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate object-level Quaia redshift-window mocks and p-values.")
    parser.add_argument("--quaia-csv", type=Path, default=DEFAULT_QUAIA)
    parser.add_argument("--randoms-fits", type=Path, default=DEFAULT_RANDOMS)
    parser.add_argument("--scan-csv", type=Path, default=DEFAULT_SCAN)
    parser.add_argument("--top-csv", type=Path, default=DEFAULT_TOP)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--window-mode", choices=["top", "all"], default="top")
    parser.add_argument("--top-n", type=int, default=25)
    parser.add_argument("--min-count", type=int, default=50)
    parser.add_argument("--n-mocks", type=int, default=200)
    parser.add_argument("--save-first-catalogs", type=int, default=1)
    parser.add_argument("--seed", type=int, default=260714)
    parser.add_argument("--cache-randoms", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    summary, pvalues, manifest = run_mocks(args)
    write_csv(args.out / "quaia_redshift_window_mock_summary.csv", summary)
    write_csv(args.out / "quaia_redshift_window_mock_pvalues.csv", pvalues)
    write_csv(args.out / "quaia_redshift_window_mock_manifest.csv", manifest)
    plot_status = write_plot(args.out / "quaia_redshift_window_mock_pvalues.png", pvalues)
    (args.out / "quaia_redshift_window_mock_config.json").write_text(json.dumps(vars(args), indent=2, default=str), encoding="utf-8")
    write_report(args.out / "quaia_redshift_window_mock_report.md", pvalues, manifest, args)
    print(f"Saved per-mock summary: {args.out / 'quaia_redshift_window_mock_summary.csv'}")
    print(f"Saved p-values: {args.out / 'quaia_redshift_window_mock_pvalues.csv'}")
    print(f"Saved report: {args.out / 'quaia_redshift_window_mock_report.md'}")
    print(plot_status)
    return 0


if __name__ == "__main__":
    sys.exit(main())
