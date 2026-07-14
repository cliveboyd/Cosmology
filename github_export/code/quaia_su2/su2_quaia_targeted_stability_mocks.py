from __future__ import annotations

import argparse
import csv
import json
import math
import sys
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


ROOT = Path(r"C:\Users\clive\Documents\Cosmology")
DEFAULT_QUAIA = ROOT / "quaia_G20p0_basic_gal.csv"
DEFAULT_RANDOMS = ROOT / "quaia_G20p0_randoms_simple.fits"
DEFAULT_OUT = ROOT / "plamb_runs" / "diagnostics" / "su2_quaia_targeted_stability_mocks_20260714"


EQ_TO_GAL = np.array(
    [
        [-0.0548755604, -0.8734370902, -0.4838350155],
        [0.4941094279, -0.4448296300, 0.7469822445],
        [-0.8676661490, -0.1980763734, 0.4559837762],
    ],
    dtype=float,
)


def gal_unit_array(l_deg: np.ndarray, b_deg: np.ndarray) -> np.ndarray:
    l = np.deg2rad(l_deg)
    b = np.deg2rad(b_deg)
    cb = np.cos(b)
    return np.column_stack([cb * np.cos(l), cb * np.sin(l), np.sin(b)])


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


def unit_from_lb(l_deg: float, b_deg: float) -> np.ndarray:
    l = math.radians(l_deg)
    b = math.radians(b_deg)
    cb = math.cos(b)
    return np.array([cb * math.cos(l), cb * math.sin(l), math.sin(b)], dtype=float)


def unit_to_lb(vec: np.ndarray) -> tuple[float, float]:
    norm = float(np.linalg.norm(vec))
    if norm <= 0.0 or not math.isfinite(norm):
        return float("nan"), float("nan")
    unit = vec / norm
    return math.degrees(math.atan2(float(unit[1]), float(unit[0]))) % 360.0, math.degrees(math.asin(float(np.clip(unit[2], -1.0, 1.0))))


def angular_sep_deg(a: np.ndarray, b: np.ndarray) -> float:
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denom <= 0.0:
        return float("nan")
    dot = float(np.clip(np.dot(a, b) / denom, -1.0, 1.0))
    return math.degrees(math.acos(dot))


def fit_redshift_dipole(z: np.ndarray, vec: np.ndarray) -> dict[str, float]:
    n = len(z)
    sx = np.sum(vec, axis=0)
    xtx = np.empty((4, 4), dtype=float)
    xtx[0, 0] = float(n)
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
        "formal_snr": snr,
        "rms_resid": float(math.sqrt(sigma2)),
    }


def load_quaia(path: Path) -> dict[str, np.ndarray]:
    df = pd.read_csv(path)
    cols = {c.upper(): c for c in df.columns}
    z = df[cols["Z"]].to_numpy(float)
    l_gal = df[cols["L_GAL"]].to_numpy(float)
    b_gal = df[cols["B_GAL"]].to_numpy(float)
    ok = np.isfinite(z) & np.isfinite(l_gal) & np.isfinite(b_gal) & (z > 0.0)
    return {"z": z[ok], "l_gal": l_gal[ok], "b_gal": b_gal[ok], "vec": gal_unit_array(l_gal[ok], b_gal[ok])}


def load_randoms(path: Path, cache_path: Path | None) -> dict[str, np.ndarray]:
    if cache_path is not None and cache_path.exists():
        data = np.load(cache_path)
        return {key: np.asarray(data[key]) for key in ["ra", "dec", "l_gal", "b_gal", "vec"]}
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


def target_windows() -> list[dict[str, Any]]:
    return [
        {"target": "broad_0p95_1p45", "z_lo": 0.95, "z_hi": 1.45, "note": "strongest stable broad control"},
        {"target": "narrow_1p15_1p35", "z_lo": 1.15, "z_hi": 1.35, "note": "stable narrow follow-up around z~1.25"},
    ]


def summarize_fits(target: dict[str, Any], fit_rows: list[dict[str, Any]], prefix: str = "") -> dict[str, Any]:
    amps = np.array([row["amp"] for row in fit_rows], dtype=float)
    snrs = np.array([row["formal_snr"] for row in fit_rows], dtype=float)
    dirs = np.array([unit_from_lb(row["l_deg"], row["b_deg"]) for row in fit_rows], dtype=float)
    max_pair_sep = 0.0
    for i in range(len(dirs)):
        for j in range(i + 1, len(dirs)):
            max_pair_sep = max(max_pair_sep, angular_sep_deg(dirs[i], dirs[j]))
    mean_dir = np.sum(dirs, axis=0)
    mean_resultant = float(np.linalg.norm(mean_dir) / max(len(dirs), 1))
    score = float(np.nanmax(snrs) * mean_resultant / (1.0 + max_pair_sep / 45.0))
    return {
        f"{prefix}target": target["target"],
        f"{prefix}z_lo": target["z_lo"],
        f"{prefix}z_hi": target["z_hi"],
        f"{prefix}n_bcuts": len(fit_rows),
        f"{prefix}N_min": int(min(row["N"] for row in fit_rows)),
        f"{prefix}N_max": int(max(row["N"] for row in fit_rows)),
        f"{prefix}amp_mean": float(np.nanmean(amps)),
        f"{prefix}amp_min": float(np.nanmin(amps)),
        f"{prefix}amp_max": float(np.nanmax(amps)),
        f"{prefix}amp_cv": float(np.nanstd(amps) / np.nanmean(amps)) if np.nanmean(amps) > 0 else float("nan"),
        f"{prefix}snr_max": float(np.nanmax(snrs)),
        f"{prefix}snr_mean": float(np.nanmean(snrs)),
        f"{prefix}max_pair_direction_sep_deg": float(max_pair_sep),
        f"{prefix}mean_resultant_length": mean_resultant,
        f"{prefix}coherence_score": score,
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


def run(args: argparse.Namespace) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    args.out.mkdir(parents=True, exist_ok=True)
    quaia = load_quaia(args.quaia_csv)
    cache_path = args.out / "quaia_randoms_galactic_cache.npz" if args.cache_randoms else None
    randoms = load_randoms(args.randoms_fits, cache_path)
    bcuts = [float(x) for x in args.bcuts]
    min_bcut = min(bcuts)
    rng = np.random.default_rng(args.seed)
    rand_pool = np.flatnonzero(np.abs(randoms["b_gal"]) >= min_bcut)
    targets = target_windows()

    observed_fits: list[dict[str, Any]] = []
    observed_summary: list[dict[str, Any]] = []
    mock_rows: list[dict[str, Any]] = []
    pvalue_rows: list[dict[str, Any]] = []

    prepared: list[dict[str, Any]] = []
    for target in targets:
        base_mask = (quaia["z"] >= target["z_lo"]) & (quaia["z"] < target["z_hi"]) & (np.abs(quaia["b_gal"]) >= min_bcut)
        base_idx = np.flatnonzero(base_mask)
        base_z = quaia["z"][base_idx]
        target_fits: list[dict[str, Any]] = []
        for bcut in bcuts:
            mask = (quaia["z"] >= target["z_lo"]) & (quaia["z"] < target["z_hi"]) & (np.abs(quaia["b_gal"]) >= bcut)
            idx = np.flatnonzero(mask)
            stats = fit_redshift_dipole(quaia["z"][idx], quaia["vec"][idx])
            row = {"target": target["target"], "z_lo": target["z_lo"], "z_hi": target["z_hi"], "bcut_deg": bcut, "N": len(idx), **stats}
            observed_fits.append(row)
            target_fits.append(row)
        obs_sum = summarize_fits(target, target_fits)
        observed_summary.append(obs_sum)
        prepared.append({"target": target, "base_z": base_z, "N_base": len(base_z)})
        print(
            f"observed {target['target']}: N_base={len(base_z)} "
            f"snr_max={obs_sum['snr_max']:.4g} max_pair_sep={obs_sum['max_pair_direction_sep_deg']:.2f} "
            f"coherence={obs_sum['coherence_score']:.4g}",
            flush=True,
        )

    replace_randoms = len(rand_pool) < max(item["N_base"] for item in prepared)
    for mock_id in range(args.n_mocks):
        for item in prepared:
            target = item["target"]
            z_mock = rng.permutation(item["base_z"])
            rand_idx = rng.choice(rand_pool, size=item["N_base"], replace=replace_randoms)
            r_b = randoms["b_gal"][rand_idx]
            r_vec = randoms["vec"][rand_idx]
            fit_rows: list[dict[str, Any]] = []
            for bcut in bcuts:
                sub = np.abs(r_b) >= bcut
                if int(np.sum(sub)) < args.min_count:
                    continue
                stats = fit_redshift_dipole(z_mock[sub], r_vec[sub])
                fit_rows.append({"target": target["target"], "bcut_deg": bcut, "N": int(np.sum(sub)), **stats})
            if len(fit_rows) != len(bcuts):
                continue
            row = {"mock_id": mock_id, **summarize_fits(target, fit_rows)}
            mock_rows.append(row)
        if (mock_id + 1) % args.progress_every == 0 or mock_id == 0:
            print(f"mock {mock_id + 1}/{args.n_mocks}", flush=True)

    mock_df = pd.DataFrame(mock_rows)
    for obs in observed_summary:
        target = obs["target"]
        m = mock_df[mock_df["target"] == target]
        amp_mean = m["amp_mean"].to_numpy(float)
        snr_max = m["snr_max"].to_numpy(float)
        pair_sep = m["max_pair_direction_sep_deg"].to_numpy(float)
        resultant = m["mean_resultant_length"].to_numpy(float)
        coherence = m["coherence_score"].to_numpy(float)
        pvalue_rows.append(
            {
                "target": target,
                "z_lo": obs["z_lo"],
                "z_hi": obs["z_hi"],
                "n_mocks": int(len(m)),
                "observed_amp_mean": obs["amp_mean"],
                "mock_amp_mean_p50": q(amp_mean, 0.50),
                "mock_amp_mean_p95": q(amp_mean, 0.95),
                "p_amp_mean_ge_observed": p_ge(amp_mean, obs["amp_mean"]),
                "observed_snr_max": obs["snr_max"],
                "mock_snr_max_p50": q(snr_max, 0.50),
                "mock_snr_max_p95": q(snr_max, 0.95),
                "p_snr_max_ge_observed": p_ge(snr_max, obs["snr_max"]),
                "observed_max_pair_direction_sep_deg": obs["max_pair_direction_sep_deg"],
                "mock_pair_sep_p05": q(pair_sep, 0.05),
                "mock_pair_sep_p50": q(pair_sep, 0.50),
                "p_pair_sep_le_observed": p_le(pair_sep, obs["max_pair_direction_sep_deg"]),
                "observed_mean_resultant_length": obs["mean_resultant_length"],
                "mock_resultant_p50": q(resultant, 0.50),
                "mock_resultant_p95": q(resultant, 0.95),
                "p_resultant_ge_observed": p_ge(resultant, obs["mean_resultant_length"]),
                "observed_coherence_score": obs["coherence_score"],
                "mock_coherence_p50": q(coherence, 0.50),
                "mock_coherence_p95": q(coherence, 0.95),
                "p_coherence_ge_observed": p_ge(coherence, obs["coherence_score"]),
                "p_joint_snr_and_pair_sep": float((1 + np.sum((snr_max >= obs["snr_max"]) & (pair_sep <= obs["max_pair_direction_sep_deg"]))) / (1 + len(m))) if len(m) else float("nan"),
            }
        )
    return observed_fits, observed_summary, mock_rows, pvalue_rows


def write_report(path: Path, pvalues: list[dict[str, Any]], observed_summary: list[dict[str, Any]], args: argparse.Namespace) -> None:
    lines = [
        "# SU2 / Quaia Targeted Stability Mock Percentiles",
        "",
        "Date: July 14, 2026",
        "",
        "## Purpose",
        "",
        "This run checks whether the observed direction/amplitude stability in the two targeted follow-up windows is unusual against nested random-sky mocks.",
        "",
        "For each target window, the mock draws random Quaia sky positions at the loosest latitude cut and evaluates the stricter latitude cuts as nested subsets. This preserves the latitude-cut correlation structure needed for a stability test.",
        "",
        "## Configuration",
        "",
        f"- n_mocks: `{args.n_mocks}`",
        f"- latitude cuts: `{', '.join(str(x) for x in args.bcuts)}`",
        f"- seed: `{args.seed}`",
        f"- Quaia table: `{args.quaia_csv}`",
        f"- random sky catalog: `{args.randoms_fits}`",
        "",
        "## Observed Stability Metrics",
        "",
        markdown_table(
            observed_summary,
            [
                "target",
                "z_lo",
                "z_hi",
                "amp_mean",
                "snr_max",
                "max_pair_direction_sep_deg",
                "mean_resultant_length",
                "coherence_score",
            ],
        ),
        "",
        "## Mock Percentiles",
        "",
        markdown_table(
            pvalues,
            [
                "target",
                "n_mocks",
                "observed_snr_max",
                "p_snr_max_ge_observed",
                "observed_max_pair_direction_sep_deg",
                "p_pair_sep_le_observed",
                "observed_coherence_score",
                "p_coherence_ge_observed",
                "p_joint_snr_and_pair_sep",
            ],
        ),
        "",
        "## Interpretation Gate",
        "",
        "- A targeted window is interesting only if its amplitude/SNR and direction clustering are unusual together.",
        "- Small `p_pair_sep_le_observed` means the observed directions are unusually tightly clustered across latitude cuts.",
        "- Small `p_joint_snr_and_pair_sep` means the mock rarely gets both comparable SNR and comparable direction clustering.",
        "- This remains exploratory because the full global look-elsewhere gate did not pass.",
        "",
        "## Outputs",
        "",
        f"- observed fits: `{args.out / 'su2_quaia_targeted_stability_observed_fits.csv'}`",
        f"- observed summary: `{args.out / 'su2_quaia_targeted_stability_observed_summary.csv'}`",
        f"- mock summary: `{args.out / 'su2_quaia_targeted_stability_mock_summary.csv'}`",
        f"- p-values: `{args.out / 'su2_quaia_targeted_stability_pvalues.csv'}`",
        f"- config: `{args.out / 'su2_quaia_targeted_stability_config.json'}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Targeted nested random-sky stability mock percentiles for SU2/Quaia windows.")
    parser.add_argument("--quaia-csv", type=Path, default=DEFAULT_QUAIA)
    parser.add_argument("--randoms-fits", type=Path, default=DEFAULT_RANDOMS)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--n-mocks", type=int, default=1500)
    parser.add_argument("--bcuts", type=float, nargs="+", default=[10, 15, 20, 25, 30, 35])
    parser.add_argument("--min-count", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=290714)
    parser.add_argument("--progress-every", type=int, default=50)
    parser.add_argument("--cache-randoms", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    observed_fits, observed_summary, mock_rows, pvalue_rows = run(args)
    write_csv(args.out / "su2_quaia_targeted_stability_observed_fits.csv", observed_fits)
    write_csv(args.out / "su2_quaia_targeted_stability_observed_summary.csv", observed_summary)
    write_csv(args.out / "su2_quaia_targeted_stability_mock_summary.csv", mock_rows)
    write_csv(args.out / "su2_quaia_targeted_stability_pvalues.csv", pvalue_rows)
    (args.out / "su2_quaia_targeted_stability_config.json").write_text(json.dumps(vars(args), indent=2, default=str), encoding="utf-8")
    write_report(args.out / "su2_quaia_targeted_stability_report.md", pvalue_rows, observed_summary, args)
    print(f"Saved observed fits: {args.out / 'su2_quaia_targeted_stability_observed_fits.csv'}", flush=True)
    print(f"Saved observed summary: {args.out / 'su2_quaia_targeted_stability_observed_summary.csv'}", flush=True)
    print(f"Saved mock summary: {args.out / 'su2_quaia_targeted_stability_mock_summary.csv'}", flush=True)
    print(f"Saved p-values: {args.out / 'su2_quaia_targeted_stability_pvalues.csv'}", flush=True)
    print(f"Saved report: {args.out / 'su2_quaia_targeted_stability_report.md'}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
