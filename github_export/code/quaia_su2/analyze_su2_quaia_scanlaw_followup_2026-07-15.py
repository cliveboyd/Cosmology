from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

import analyze_su2_quaia_external_dust_gate as ext


ROOT = Path(r"C:\Users\clive\Documents\Cosmology")
DEFAULT_OUT = ROOT / "plamb_runs" / "diagnostics" / "su2_quaia_scanlaw_followup_20260715"


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


def precompute_baseline(vec: np.ndarray) -> dict[str, np.ndarray]:
    n = len(vec)
    sx = np.sum(vec, axis=0)
    xtx = np.empty((4, 4), dtype=float)
    xtx[0, 0] = float(n)
    xtx[0, 1:] = sx
    xtx[1:, 0] = sx
    xtx[1:, 1:] = vec.T @ vec
    return {"pinv_xtx": np.linalg.pinv(xtx)}


def fit_fast(z: np.ndarray, vec: np.ndarray, geom: dict[str, np.ndarray]) -> dict[str, float]:
    xty = np.empty(4, dtype=float)
    xty[0] = float(np.sum(z))
    xty[1:] = vec.T @ z
    pinv_xtx = geom["pinv_xtx"]
    beta = pinv_xtx @ xty
    yty = float(np.dot(z, z))
    sse = max(yty - float(beta @ xty), 0.0)
    dof = max(len(z) - 4, 1)
    sigma2 = sse / dof
    dvec = beta[1:4]
    amp = float(np.linalg.norm(dvec))
    l_deg, b_deg = ext.unit_to_lb(dvec)
    try:
        cov_d = sigma2 * pinv_xtx[1:, 1:]
        snr = float(math.sqrt(max(float(dvec @ np.linalg.pinv(cov_d) @ dvec), 0.0)))
    except Exception:
        snr = float("nan")
    return {
        "amp": amp,
        "l_deg": l_deg,
        "b_deg": b_deg,
        "formal_snr": snr,
        "rms_resid": float(math.sqrt(max(sigma2, 0.0))),
    }


def summarize_fits(label: str, fits: list[dict[str, Any]]) -> dict[str, Any]:
    amps = np.array([row["amp"] for row in fits], dtype=float)
    snrs = np.array([row["formal_snr"] for row in fits], dtype=float)
    dirs = np.array([unit_from_lb(row["l_deg"], row["b_deg"]) for row in fits], dtype=float)
    max_pair_sep = 0.0
    for i in range(len(dirs)):
        for j in range(i + 1, len(dirs)):
            max_pair_sep = max(max_pair_sep, angular_sep_deg(dirs[i], dirs[j]))
    mean_dir = np.sum(dirs, axis=0)
    mean_resultant = float(np.linalg.norm(mean_dir) / max(len(dirs), 1))
    coherence_score = float(np.nanmax(snrs) * mean_resultant / (1.0 + max_pair_sep / 45.0))
    return {
        "label": label,
        "n_bcuts": len(fits),
        "N_min": int(min(row["N"] for row in fits)),
        "N_max": int(max(row["N"] for row in fits)),
        "amp_mean": float(np.nanmean(amps)),
        "amp_min": float(np.nanmin(amps)),
        "amp_max": float(np.nanmax(amps)),
        "snr_max": float(np.nanmax(snrs)),
        "snr_mean": float(np.nanmean(snrs)),
        "max_pair_direction_sep_deg": float(max_pair_sep),
        "mean_resultant_length": mean_resultant,
        "coherence_score": coherence_score,
    }


def p_ge(values: np.ndarray, observed: float) -> float:
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return float("nan")
    return float((1 + np.sum(values >= observed)) / (1 + len(values)))


def p_le(values: np.ndarray, observed: float) -> float:
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return float("nan")
    return float((1 + np.sum(values <= observed)) / (1 + len(values)))


def qtile(values: np.ndarray, prob: float) -> float:
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return float("nan")
    return float(np.quantile(values, prob))


def quartile_edges(values: np.ndarray) -> np.ndarray:
    finite = values[np.isfinite(values)]
    return np.quantile(finite, [0.0, 0.25, 0.5, 0.75, 1.0])


def quartile_index(values: np.ndarray, edges: np.ndarray) -> np.ndarray:
    qidx = np.searchsorted(edges[1:-1], values, side="right") + 1
    return np.clip(qidx, 1, 4)


def phase_quartile(cos2: np.ndarray, sin2: np.ndarray) -> np.ndarray:
    phase = (np.arctan2(sin2, cos2) + 2.0 * np.pi) % (2.0 * np.pi)
    return np.floor(phase / (0.5 * np.pi)).astype(np.int64) + 1


def standardize(values: np.ndarray) -> np.ndarray:
    med = np.nanmedian(values)
    scale = np.nanstd(values)
    if not math.isfinite(scale) or scale <= 0.0:
        scale = 1.0
    return (values - med) / scale


def fit_template_residual_dipole(z: np.ndarray, vec: np.ndarray, controls: list[np.ndarray]) -> dict[str, float]:
    columns = [np.ones(len(z))]
    columns.extend(standardize(control) for control in controls)
    x = np.column_stack(columns)
    beta, *_ = np.linalg.lstsq(x, z, rcond=None)
    resid = z - x @ beta
    return fit_fast(resid, vec, precompute_baseline(vec))


def build_scan_bins(q: dict[str, np.ndarray], mask: np.ndarray) -> tuple[np.ndarray, list[dict[str, Any]]]:
    scan_count = q["gaia_scan_count_log1p_dr3"][mask]
    scan_resultant = q["gaia_scan_angle_resultant_dr3"][mask]
    scan_cos = q["gaia_scan_angle_cos2_mean_dr3"][mask]
    scan_sin = q["gaia_scan_angle_sin2_mean_dr3"][mask]
    count_edges = quartile_edges(scan_count)
    result_edges = quartile_edges(scan_resultant)
    count_q = quartile_index(scan_count, count_edges)
    result_q = quartile_index(scan_resultant, result_edges)
    phase_q = phase_quartile(scan_cos, scan_sin)
    bins = (count_q - 1) * 16 + (result_q - 1) * 4 + phase_q
    rows = []
    for bin_id in sorted(np.unique(bins)):
        sub = bins == bin_id
        rows.append(
            {
                "scan_bin": int(bin_id),
                "N": int(np.sum(sub)),
                "scan_count_q": int(((bin_id - 1) // 16) + 1),
                "scan_resultant_q": int((((bin_id - 1) % 16) // 4) + 1),
                "scan_phase_q": int(((bin_id - 1) % 4) + 1),
                "mean_scan_count_log1p": float(np.mean(scan_count[sub])),
                "mean_scan_resultant": float(np.mean(scan_resultant[sub])),
            }
        )
    return bins.astype(np.int64), rows


def fit_across_bcuts(
    label: str,
    z: np.ndarray,
    b: np.ndarray,
    vec: np.ndarray,
    bcuts: list[float],
    geoms: dict[float, dict[str, np.ndarray]] | None = None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows = []
    for bcut in bcuts:
        mask = np.abs(b) >= bcut
        geom = geoms[bcut] if geoms is not None else precompute_baseline(vec[mask])
        fit = fit_fast(z[mask], vec[mask], geom)
        rows.append({"label": label, "bcut_deg": bcut, "N": int(np.sum(mask)), **fit})
    return rows, summarize_fits(label, rows)


def scanlaw_strata(q: dict[str, np.ndarray], base_mask: np.ndarray, bcuts: list[float]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    summary: list[dict[str, Any]] = []
    variables = [
        ("gaia_scan_count_log1p_dr3", "Gaia scan count"),
        ("gaia_scan_angle_resultant_dr3", "Gaia scan-angle resultant"),
    ]
    for name, label in variables:
        values = q[name][base_mask]
        edges = quartile_edges(values)
        qidx = quartile_index(values, edges)
        quartile_summaries = []
        for qi in range(1, 5):
            sub_global = np.flatnonzero(base_mask)[qidx == qi]
            z = q["z"][sub_global]
            b = q["b"][sub_global]
            vec = q["vec"][sub_global]
            fit_rows, fit_summary = fit_across_bcuts(f"{name}_q{qi}", z, b, vec, bcuts)
            for row in fit_rows:
                row.update(
                    {
                        "strata_variable": name,
                        "strata_label": label,
                        "quartile": qi,
                        "quartile_min": float(edges[qi - 1]),
                        "quartile_max": float(edges[qi]),
                    }
                )
            rows.extend(fit_rows)
            fit_summary.update(
                {
                    "strata_variable": name,
                    "strata_label": label,
                    "quartile": qi,
                    "quartile_min": float(edges[qi - 1]),
                    "quartile_max": float(edges[qi]),
                }
            )
            quartile_summaries.append(fit_summary)
            summary.append(fit_summary)

        q1 = quartile_summaries[0]
        q4 = quartile_summaries[-1]
        dirs = [unit_from_lb(fit["l_deg"], fit["b_deg"]) for fit in rows if fit.get("strata_variable") == name]
        max_sep = 0.0
        for i in range(len(dirs)):
            for j in range(i + 1, len(dirs)):
                max_sep = max(max_sep, angular_sep_deg(dirs[i], dirs[j]))
        summary.append(
            {
                "label": f"{name}_quartile_gradient",
                "strata_variable": name,
                "strata_label": label,
                "quartile": "Q4/Q1",
                "amp_q4_over_q1": q4["amp_mean"] / q1["amp_mean"] if q1["amp_mean"] else float("nan"),
                "snr_q4_minus_q1": q4["snr_max"] - q1["snr_max"],
                "max_fit_direction_sep_deg": max_sep,
                "q1_N_min": q1["N_min"],
                "q4_N_min": q4["N_min"],
            }
        )
    return rows, summary


def residual_tests(
    q: dict[str, np.ndarray],
    base_mask: np.ndarray,
    bcuts: list[float],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    indices = np.flatnonzero(base_mask)
    z = q["z"][indices]
    b = q["b"][indices]
    vec = q["vec"][indices]
    external_controls = [
        q["ebv_sfd"][indices],
        q["sfd_log1p"][indices],
        q["sfd_sq"][indices],
        q["selection_T"][indices],
        q["random_density_log1p"][indices],
        q["gaia_scan_count_log1p_dr3"][indices],
        q["gaia_scan_angle_cos2_mean_dr3"][indices],
        q["gaia_scan_angle_sin2_mean_dr3"][indices],
        q["gaia_scan_angle_resultant_dr3"][indices],
    ]
    all_proxy_controls = external_controls + [
        q["zerr"][indices],
        q["g"][indices],
        q["w1"][indices],
        q["w2"][indices],
        q["w1_w2"][indices],
        q["pmra_error"][indices],
        q["pmdec_error"][indices],
    ]
    test_groups = [
        ("baseline_no_residualization", []),
        ("residual_after_sfd_selection_gaia", external_controls),
        ("residual_after_all_available_proxy_stack", all_proxy_controls),
    ]
    rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    baseline_summary: dict[str, Any] | None = None
    for label, controls in test_groups:
        fit_rows = []
        for bcut in bcuts:
            sub = np.abs(b) >= bcut
            if controls:
                fit = fit_template_residual_dipole(z[sub], vec[sub], [control[sub] for control in controls])
            else:
                fit = fit_fast(z[sub], vec[sub], precompute_baseline(vec[sub]))
            fit_rows.append({"label": label, "bcut_deg": bcut, "N": int(np.sum(sub)), **fit})
        summary = summarize_fits(label, fit_rows)
        if baseline_summary is None:
            baseline_summary = summary
            summary["amp_ratio_vs_baseline"] = 1.0
            summary["snr_delta_vs_baseline"] = 0.0
            summary["pair_sep_delta_vs_baseline"] = 0.0
        else:
            summary["amp_ratio_vs_baseline"] = summary["amp_mean"] / baseline_summary["amp_mean"]
            summary["snr_delta_vs_baseline"] = summary["snr_max"] - baseline_summary["snr_max"]
            summary["pair_sep_delta_vs_baseline"] = summary["max_pair_direction_sep_deg"] - baseline_summary["max_pair_direction_sep_deg"]
        rows.extend(fit_rows)
        summary_rows.append(summary)
    return rows, summary_rows


def scanlaw_preserving_null(
    q: dict[str, np.ndarray],
    base_mask: np.ndarray,
    bcuts: list[float],
    n_mocks: int,
    seed: int,
    progress_every: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    indices = np.flatnonzero(base_mask)
    z = q["z"][indices]
    b = q["b"][indices]
    vec = q["vec"][indices]
    masks = {bcut: np.abs(b) >= bcut for bcut in bcuts}
    geoms = {bcut: precompute_baseline(vec[masks[bcut]]) for bcut in bcuts}
    obs_rows, obs_summary = fit_across_bcuts("observed", z, b, vec, bcuts, geoms)
    bins, bin_rows = build_scan_bins(q, base_mask)
    bin_indices = [np.flatnonzero(bins == bin_id) for bin_id in sorted(np.unique(bins))]
    rng = np.random.default_rng(seed)
    mock_rows: list[dict[str, Any]] = []
    for mock_id in range(n_mocks):
        z_mock = z.copy()
        for sub in bin_indices:
            if len(sub) > 1:
                z_mock[sub] = rng.permutation(z_mock[sub])
        _, summary = fit_across_bcuts("mock", z_mock, b, vec, bcuts, geoms)
        summary["mock_id"] = mock_id
        mock_rows.append(summary)
        if progress_every and ((mock_id + 1) % progress_every == 0 or mock_id == 0):
            print(f"scanlaw-preserving shuffle {mock_id + 1}/{n_mocks}", flush=True)

    amp = np.array([row["amp_mean"] for row in mock_rows], dtype=float)
    snr = np.array([row["snr_max"] for row in mock_rows], dtype=float)
    sep = np.array([row["max_pair_direction_sep_deg"] for row in mock_rows], dtype=float)
    coh = np.array([row["coherence_score"] for row in mock_rows], dtype=float)
    pvalue = {
        "null_kind": "scanlaw_count_resultant_phase_zshuffle",
        "n_mocks": len(mock_rows),
        "n_scan_bins": len(bin_indices),
        "observed_amp_mean": obs_summary["amp_mean"],
        "mock_amp_mean_p50": qtile(amp, 0.50),
        "mock_amp_mean_p95": qtile(amp, 0.95),
        "mock_amp_mean_p99": qtile(amp, 0.99),
        "p_amp_mean_ge_observed": p_ge(amp, obs_summary["amp_mean"]),
        "observed_snr_max": obs_summary["snr_max"],
        "mock_snr_max_p50": qtile(snr, 0.50),
        "mock_snr_max_p95": qtile(snr, 0.95),
        "mock_snr_max_p99": qtile(snr, 0.99),
        "p_snr_max_ge_observed": p_ge(snr, obs_summary["snr_max"]),
        "observed_max_pair_direction_sep_deg": obs_summary["max_pair_direction_sep_deg"],
        "mock_pair_sep_p05": qtile(sep, 0.05),
        "mock_pair_sep_p50": qtile(sep, 0.50),
        "p_pair_sep_le_observed": p_le(sep, obs_summary["max_pair_direction_sep_deg"]),
        "observed_coherence_score": obs_summary["coherence_score"],
        "mock_coherence_p50": qtile(coh, 0.50),
        "mock_coherence_p95": qtile(coh, 0.95),
        "mock_coherence_p99": qtile(coh, 0.99),
        "p_coherence_ge_observed": p_ge(coh, obs_summary["coherence_score"]),
        "p_joint_snr_and_pair_sep": float((1 + np.sum((snr >= obs_summary["snr_max"]) & (sep <= obs_summary["max_pair_direction_sep_deg"]))) / (1 + len(mock_rows))),
    }
    return obs_rows, mock_rows, [pvalue], bin_rows


def write_report(
    out_dir: Path,
    pvalues: list[dict[str, Any]],
    strata_summary: list[dict[str, Any]],
    residual_summary: list[dict[str, Any]],
    args: argparse.Namespace,
) -> None:
    p = pvalues[0]
    residual_external = next(row for row in residual_summary if row["label"] == "residual_after_sfd_selection_gaia")
    if p["p_joint_snr_and_pair_sep"] <= 0.01:
        null_line = "The locked mode remains rare under scan-law-preserving redshift shuffles."
    else:
        null_line = "The locked mode is not rare under scan-law-preserving redshift shuffles."
    if residual_external["amp_ratio_vs_baseline"] >= 0.8:
        residual_line = "The residual-mode test still retains most of the baseline amplitude after SFD + selection/depth + Gaia scan-law regression."
    else:
        residual_line = "The residual-mode test loses a material fraction of the baseline amplitude after SFD + selection/depth + Gaia scan-law regression."

    report = [
        "# SU2 / Quaia Gaia Scan-Law Follow-Up",
        "",
        "Date: July 15, 2026",
        "",
        "## Purpose",
        "",
        "This follows the external-template gate where the Gaia scan-law controls weakened and rotated the locked `0.95 <= z < 1.45` angular mode.",
        "",
        "## Tests",
        "",
        "1. Scan-law preserving null mocks: shuffle redshifts inside bins of Gaia scan count, scan-angle resultant, and scan-angle phase.",
        "2. Scan-law strata: fit the locked mode separately in scan-count and scan-angle-resultant quartiles.",
        "3. Residual mode: regress out SFD dust, Quaia selection/depth, and Gaia scan-law templates before refitting the dipole.",
        "",
        "## Scan-Law Preserving Null",
        "",
        markdown_table(
            pvalues,
            [
                "null_kind",
                "n_mocks",
                "n_scan_bins",
                "observed_snr_max",
                "p_snr_max_ge_observed",
                "p_pair_sep_le_observed",
                "p_coherence_ge_observed",
                "p_joint_snr_and_pair_sep",
                "mock_snr_max_p99",
            ],
        ),
        "",
        f"- {null_line}",
        "",
        "## Scan-Law Strata Summary",
        "",
        markdown_table(
            strata_summary,
            [
                "strata_variable",
                "quartile",
                "amp_mean",
                "snr_max",
                "max_pair_direction_sep_deg",
                "N_min",
                "amp_q4_over_q1",
                "max_fit_direction_sep_deg",
            ],
        ),
        "",
        "## Residual Mode Summary",
        "",
        markdown_table(
            residual_summary,
            [
                "label",
                "amp_ratio_vs_baseline",
                "snr_max",
                "snr_delta_vs_baseline",
                "max_pair_direction_sep_deg",
                "mean_resultant_length",
                "N_min",
                "N_max",
            ],
        ),
        "",
        f"- {residual_line}",
        "",
        "## Interpretation",
        "",
        "The promotion gate is not simply whether the signal survives dust and selection/depth. It must also survive an explicit Gaia scan-law null and retain a coherent residual mode after scan-law regression.",
        "",
        "## Configuration",
        "",
        f"- mocks: `{args.n_mocks}`",
        f"- redshift window: `{args.z_min} <= z < {args.z_max}`",
        f"- latitude cuts: `{', '.join(str(x) for x in args.b_cuts)}`",
        f"- seed: `{args.seed}`",
    ]
    (out_dir / "su2_quaia_scanlaw_followup_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    readout = [
        "# SU2 / Quaia Gaia Scan-Law Follow-Up Readout",
        "",
        "Date: July 15, 2026",
        "",
        "## Bottom Line",
        "",
        null_line,
        "",
        residual_line,
        "",
        "## Key P-Values",
        "",
        markdown_table(
            pvalues,
            [
                "n_mocks",
                "n_scan_bins",
                "p_snr_max_ge_observed",
                "p_pair_sep_le_observed",
                "p_coherence_ge_observed",
                "p_joint_snr_and_pair_sep",
            ],
        ),
        "",
        "## Residual Summary",
        "",
        markdown_table(
            residual_summary,
            [
                "label",
                "amp_ratio_vs_baseline",
                "snr_max",
                "max_pair_direction_sep_deg",
            ],
        ),
    ]
    (out_dir / "su2_quaia_scanlaw_followup_readout.md").write_text("\n".join(readout) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Gaia scan-law follow-up gates for SU2/Quaia.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--n-mocks", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=150715)
    parser.add_argument("--z-min", type=float, default=0.95)
    parser.add_argument("--z-max", type=float, default=1.45)
    parser.add_argument("--b-cuts", type=float, nargs="+", default=[15.0, 20.0, 25.0, 30.0, 35.0])
    parser.add_argument("--progress-every", type=int, default=250)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    gaiascanlaw_data_dir = ext.resolve_gaiascanlaw_data_dir()
    q = ext.load_data(ext.DEFAULT_QUAIA, ext.DEFAULT_SELECTION, ext.DEFAULT_RANDOMS, ext.DEFAULT_SFD_DIR, gaiascanlaw_data_dir)
    q["vec"] = q["vec"] if "vec" in q else q["vec_gal"]

    base_mask = (
        (q["z"] >= args.z_min)
        & (q["z"] < args.z_max)
        & (np.abs(q["b"]) >= min(args.b_cuts))
        & np.isfinite(q["z"])
        & np.isfinite(q["b"])
        & np.isfinite(q["gaia_scan_count_log1p_dr3"])
        & np.isfinite(q["gaia_scan_angle_resultant_dr3"])
        & np.isfinite(q["gaia_scan_angle_cos2_mean_dr3"])
        & np.isfinite(q["gaia_scan_angle_sin2_mean_dr3"])
    )

    obs_rows, mock_rows, pvalues, bin_rows = scanlaw_preserving_null(
        q,
        base_mask,
        args.b_cuts,
        args.n_mocks,
        args.seed,
        args.progress_every,
    )
    strata_rows, strata_summary = scanlaw_strata(q, base_mask, args.b_cuts)
    residual_rows, residual_summary = residual_tests(q, base_mask, args.b_cuts)

    write_csv(args.out_dir / "su2_quaia_scanlaw_followup_observed_fits.csv", obs_rows)
    write_csv(args.out_dir / "su2_quaia_scanlaw_followup_mock_summary.csv", mock_rows)
    write_csv(args.out_dir / "su2_quaia_scanlaw_followup_pvalues.csv", pvalues)
    write_csv(args.out_dir / "su2_quaia_scanlaw_followup_scan_bins.csv", bin_rows)
    write_csv(args.out_dir / "su2_quaia_scanlaw_followup_strata_fit_rows.csv", strata_rows)
    write_csv(args.out_dir / "su2_quaia_scanlaw_followup_strata_summary.csv", strata_summary)
    write_csv(args.out_dir / "su2_quaia_scanlaw_followup_residual_fit_rows.csv", residual_rows)
    write_csv(args.out_dir / "su2_quaia_scanlaw_followup_residual_summary.csv", residual_summary)

    config = {
        "date": "2026-07-15",
        "z_min": args.z_min,
        "z_max": args.z_max,
        "b_cuts": args.b_cuts,
        "n_mocks": args.n_mocks,
        "seed": args.seed,
        "scanlaw_data_dir": str(gaiascanlaw_data_dir) if gaiascanlaw_data_dir else None,
        "scan_bins": "quartiles of Gaia scan count log1p, quartiles of scan-angle resultant, and four scan-angle phase quadrants",
    }
    (args.out_dir / "su2_quaia_scanlaw_followup_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")
    write_report(args.out_dir, pvalues, strata_summary, residual_summary, args)

    print(f"Saved readout: {args.out_dir / 'su2_quaia_scanlaw_followup_readout.md'}", flush=True)
    print(f"Saved report: {args.out_dir / 'su2_quaia_scanlaw_followup_report.md'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
