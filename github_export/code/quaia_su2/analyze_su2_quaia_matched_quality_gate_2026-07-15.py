from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures, StandardScaler

try:
    from astropy.io import fits
except Exception as exc:  # pragma: no cover
    fits = None
    FITS_IMPORT_ERROR = exc
else:
    FITS_IMPORT_ERROR = None


ROOT = Path(r"C:\Users\clive\Documents\Cosmology")
DEFAULT_QUAIA = ROOT / "quaia_G20.0.fits"
DEFAULT_OUT = ROOT / "plamb_runs" / "diagnostics" / "su2_quaia_matched_quality_gate_20260715"


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


def gal_unit_array(l_deg: np.ndarray, b_deg: np.ndarray) -> np.ndarray:
    l = np.deg2rad(l_deg)
    b = np.deg2rad(b_deg)
    cb = np.cos(b)
    return np.column_stack([cb * np.cos(l), cb * np.sin(l), np.sin(b)])


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
    l_deg = math.degrees(math.atan2(float(unit[1]), float(unit[0]))) % 360.0
    b_deg = math.degrees(math.asin(float(np.clip(unit[2], -1.0, 1.0))))
    return l_deg, b_deg


def angular_sep_deg(a: np.ndarray, b: np.ndarray) -> float:
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denom <= 0.0:
        return float("nan")
    dot = float(np.clip(np.dot(a, b) / denom, -1.0, 1.0))
    return math.degrees(math.acos(dot))


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    return float(np.sum(weights * values) / np.sum(weights))


def weighted_var(values: np.ndarray, weights: np.ndarray) -> float:
    mu = weighted_mean(values, weights)
    return float(np.sum(weights * (values - mu) ** 2) / np.sum(weights))


def weighted_ks(x0: np.ndarray, w0: np.ndarray, x1: np.ndarray, w1: np.ndarray) -> float:
    x0 = np.asarray(x0, dtype=float)
    x1 = np.asarray(x1, dtype=float)
    w0 = np.asarray(w0, dtype=float)
    w1 = np.asarray(w1, dtype=float)
    mask0 = np.isfinite(x0) & np.isfinite(w0) & (w0 > 0.0)
    mask1 = np.isfinite(x1) & np.isfinite(w1) & (w1 > 0.0)
    x0, w0 = x0[mask0], w0[mask0]
    x1, w1 = x1[mask1], w1[mask1]
    if len(x0) == 0 or len(x1) == 0:
        return float("nan")
    order0 = np.argsort(x0)
    order1 = np.argsort(x1)
    x0, w0 = x0[order0], w0[order0] / np.sum(w0)
    x1, w1 = x1[order1], w1[order1] / np.sum(w1)
    grid = np.unique(np.concatenate([x0, x1]))
    c0 = np.cumsum(w0)
    c1 = np.cumsum(w1)
    idx0 = np.searchsorted(x0, grid, side="right") - 1
    idx1 = np.searchsorted(x1, grid, side="right") - 1
    f0 = np.where(idx0 >= 0, c0[np.clip(idx0, 0, len(c0) - 1)], 0.0)
    f1 = np.where(idx1 >= 0, c1[np.clip(idx1, 0, len(c1) - 1)], 0.0)
    return float(np.max(np.abs(f0 - f1)))


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


def qtile(values: np.ndarray, prob: float) -> float:
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return float("nan")
    return float(np.quantile(values, prob))


def load_quaia(path: Path) -> dict[str, np.ndarray]:
    if fits is None:
        raise RuntimeError(f"astropy.io.fits unavailable: {FITS_IMPORT_ERROR}")
    with fits.open(path, memmap=True) as hdul:
        data = hdul[1].data
        q = {
            "z": np.asarray(data["redshift_quaia"], dtype=float),
            "zerr": np.asarray(data["redshift_quaia_err"], dtype=float),
            "l": np.asarray(data["l"], dtype=float),
            "b": np.asarray(data["b"], dtype=float),
            "g": np.asarray(data["phot_g_mean_mag"], dtype=float),
            "bp_rp": np.asarray(data["phot_bp_mean_mag"], dtype=float) - np.asarray(data["phot_rp_mean_mag"], dtype=float),
            "w1": np.asarray(data["mag_w1_vg"], dtype=float),
            "w2": np.asarray(data["mag_w2_vg"], dtype=float),
            "w1_w2": np.asarray(data["mag_w1_vg"], dtype=float) - np.asarray(data["mag_w2_vg"], dtype=float),
            "pmra_error": np.asarray(data["pmra_error"], dtype=float),
            "pmdec_error": np.asarray(data["pmdec_error"], dtype=float),
        }
    q["vec"] = gal_unit_array(q["l"], q["b"])
    return q


def quality_groups() -> list[dict[str, Any]]:
    return [
        {"group": "zerr_only", "variables": ["zerr"]},
        {"group": "gaia_colour_quality", "variables": ["g", "bp_rp"]},
        {"group": "wise_colour_quality", "variables": ["w1", "w2", "w1_w2"]},
        {"group": "pm_error_quality", "variables": ["pmra_error", "pmdec_error"]},
        {"group": "all_catalogue_quality", "variables": ["zerr", "g", "bp_rp", "w1", "w2", "w1_w2", "pmra_error", "pmdec_error"]},
    ]


def fit_weighted_dipole(z: np.ndarray, vec: np.ndarray, weights: np.ndarray | None = None) -> dict[str, float]:
    n = len(z)
    x = np.column_stack([np.ones(n), vec[:, 0], vec[:, 1], vec[:, 2]])
    if weights is None:
        weights = np.ones(n, dtype=float)
    weights = np.asarray(weights, dtype=float)
    weights = np.where(np.isfinite(weights) & (weights > 0.0), weights, 0.0)
    if np.sum(weights > 0.0) < 5:
        return {"a0": float("nan"), "amp": float("nan"), "l_deg": float("nan"), "b_deg": float("nan"), "formal_snr": float("nan"), "rms_resid": float("nan")}
    w = weights / np.mean(weights[weights > 0.0])
    xw = x * w[:, None]
    xtwx = x.T @ xw
    xtwz = x.T @ (w * z)
    pinv = np.linalg.pinv(xtwx)
    beta = pinv @ xtwz
    resid = z - x @ beta
    ess = float((np.sum(w) ** 2) / np.sum(w * w))
    dof = max(ess - x.shape[1], 1.0)
    sigma2 = float(np.sum(w * resid * resid) / dof)
    dvec = beta[1:4]
    amp = float(np.linalg.norm(dvec))
    l_deg, b_deg = unit_to_lb(dvec)
    try:
        cov_d = sigma2 * pinv[1:4, 1:4]
        snr = float(math.sqrt(max(float(dvec @ np.linalg.pinv(cov_d) @ dvec), 0.0)))
    except Exception:
        snr = float("nan")
    return {
        "a0": float(beta[0]),
        "amp": amp,
        "l_deg": l_deg,
        "b_deg": b_deg,
        "formal_snr": snr,
        "rms_resid": float(math.sqrt(max(sigma2, 0.0))),
    }


def summarize_fits(label: str, fits: list[dict[str, Any]], baseline_fits: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    amps = np.array([row["amp"] for row in fits], dtype=float)
    snrs = np.array([row["formal_snr"] for row in fits], dtype=float)
    dirs = np.array([unit_from_lb(row["l_deg"], row["b_deg"]) for row in fits], dtype=float)
    max_pair_sep = 0.0
    for i in range(len(dirs)):
        for j in range(i + 1, len(dirs)):
            max_pair_sep = max(max_pair_sep, angular_sep_deg(dirs[i], dirs[j]))
    mean_dir = np.sum(dirs, axis=0)
    mean_resultant = float(np.linalg.norm(mean_dir) / max(len(dirs), 1))
    row = {
        "fit_group": label,
        "n_bcuts": len(fits),
        "N_min": int(min(fit["N"] for fit in fits)),
        "N_max": int(max(fit["N"] for fit in fits)),
        "ess_min": float(min(fit.get("ess", fit["N"]) for fit in fits)),
        "ess_max": float(max(fit.get("ess", fit["N"]) for fit in fits)),
        "amp_mean": float(np.nanmean(amps)),
        "amp_min": float(np.nanmin(amps)),
        "amp_max": float(np.nanmax(amps)),
        "snr_max": float(np.nanmax(snrs)),
        "snr_mean": float(np.nanmean(snrs)),
        "max_pair_direction_sep_deg": float(max_pair_sep),
        "mean_resultant_length": mean_resultant,
        "coherence_score": float(np.nanmax(snrs) * mean_resultant / (1.0 + max_pair_sep / 45.0)),
    }
    if baseline_fits is None:
        row["amp_ratio_vs_baseline"] = 1.0
        row["max_direction_sep_vs_baseline_deg"] = 0.0
        row["snr_delta_vs_baseline"] = 0.0
    else:
        base_amps = np.array([fit["amp"] for fit in baseline_fits], dtype=float)
        base_dirs = np.array([unit_from_lb(fit["l_deg"], fit["b_deg"]) for fit in baseline_fits], dtype=float)
        row["amp_ratio_vs_baseline"] = row["amp_mean"] / float(np.nanmean(base_amps))
        row["snr_delta_vs_baseline"] = row["snr_max"] - float(np.nanmax([fit["formal_snr"] for fit in baseline_fits]))
        seps = [angular_sep_deg(dirs[i], base_dirs[i]) for i in range(min(len(dirs), len(base_dirs)))]
        row["max_direction_sep_vs_baseline_deg"] = float(np.nanmax(seps))
    return row


def build_features(q: dict[str, np.ndarray], idx: np.ndarray, variables: list[str], poly_degree: int) -> np.ndarray:
    raw = np.column_stack([q[name][idx] for name in variables]).astype(float)
    if poly_degree <= 1:
        return raw
    # Keep interactions out for speed and stability; squares capture non-linear quality tails.
    poly = PolynomialFeatures(degree=poly_degree, include_bias=False, interaction_only=False)
    return poly.fit_transform(raw)


def overlap_weights(features: np.ndarray, hemisphere: np.ndarray, seed: int) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    model = make_pipeline(
        StandardScaler(),
        LogisticRegression(max_iter=1000, C=1.0, solver="lbfgs", random_state=seed),
    )
    model.fit(features, hemisphere.astype(int))
    p = model.predict_proba(features)[:, 1]
    p = np.clip(p, 0.02, 0.98)
    weights = np.where(hemisphere, 1.0 - p, p)
    weights /= np.mean(weights)
    info = {
        "propensity_min": float(np.min(p)),
        "propensity_p01": float(np.quantile(p, 0.01)),
        "propensity_p50": float(np.quantile(p, 0.50)),
        "propensity_p99": float(np.quantile(p, 0.99)),
        "propensity_max": float(np.max(p)),
        "weight_min": float(np.min(weights)),
        "weight_p99": float(np.quantile(weights, 0.99)),
        "weight_max": float(np.max(weights)),
        "ess": float((np.sum(weights) ** 2) / np.sum(weights * weights)),
    }
    return weights, p, info


def balance_rows(
    q: dict[str, np.ndarray],
    idx: np.ndarray,
    variables: list[str],
    hemisphere: np.ndarray,
    weights: np.ndarray,
    fit_group: str,
    bcut: float,
) -> list[dict[str, Any]]:
    rows = []
    base_weights = np.ones(len(idx), dtype=float)
    for name in variables:
        values = q[name][idx]
        h0 = ~hemisphere
        h1 = hemisphere
        pooled_sd = math.sqrt(max((np.nanvar(values[h0]) + np.nanvar(values[h1])) / 2.0, 1e-30))
        before_smd = (float(np.nanmean(values[h1])) - float(np.nanmean(values[h0]))) / pooled_sd
        w0 = weights[h0]
        w1 = weights[h1]
        after_var = (weighted_var(values[h0], w0) + weighted_var(values[h1], w1)) / 2.0
        after_sd = math.sqrt(max(after_var, 1e-30))
        after_smd = (weighted_mean(values[h1], w1) - weighted_mean(values[h0], w0)) / after_sd
        rows.append(
            {
                "fit_group": fit_group,
                "bcut_deg": bcut,
                "quality_variable": name,
                "smd_before": before_smd,
                "smd_after": after_smd,
                "abs_smd_before": abs(before_smd),
                "abs_smd_after": abs(after_smd),
                "weighted_ks_before": weighted_ks(values[h0], base_weights[h0], values[h1], base_weights[h1]),
                "weighted_ks_after": weighted_ks(values[h0], w0, values[h1], w1),
            }
        )
    return rows


def propensity_strata(propensity: np.ndarray, n_strata: int) -> np.ndarray:
    edges = np.unique(np.quantile(propensity[np.isfinite(propensity)], np.linspace(0, 1, n_strata + 1)))
    if len(edges) <= 2:
        return np.zeros(len(propensity), dtype=int)
    return np.searchsorted(edges[1:-1], propensity, side="right")


def run(args: argparse.Namespace) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    q = load_quaia(args.quaia_fits)
    bcuts = [float(x) for x in args.bcuts]
    base = (
        (q["z"] >= args.z_min)
        & (q["z"] < args.z_max)
        & np.isfinite(q["z"])
        & np.isfinite(q["l"])
        & np.isfinite(q["b"])
    )
    fit_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    balance: list[dict[str, Any]] = []
    mock_rows: list[dict[str, Any]] = []
    pvalue_rows: list[dict[str, Any]] = []
    rng = np.random.default_rng(args.seed)

    baseline_fits: list[dict[str, Any]] = []
    for bcut in bcuts:
        idx = np.flatnonzero(base & (np.abs(q["b"]) >= bcut))
        fit = fit_weighted_dipole(q["z"][idx], q["vec"][idx])
        row = {
            "fit_group": "baseline_unweighted",
            "bcut_deg": bcut,
            "N": int(len(idx)),
            "ess": float(len(idx)),
            "variables": "",
            **fit,
            "amp_ratio_vs_baseline_same_bcut": 1.0,
            "direction_sep_vs_baseline_same_bcut_deg": 0.0,
        }
        fit_rows.append(row)
        baseline_fits.append(row)
    baseline_summary = summarize_fits("baseline_unweighted", baseline_fits)
    summary_rows.append(baseline_summary)

    group_state: dict[str, dict[str, Any]] = {}
    for group in quality_groups():
        group_name = group["group"]
        variables = group["variables"]
        weighted_fits: list[dict[str, Any]] = []
        group_state[group_name] = {"bcuts": {}}
        for bcut, base_fit in zip(bcuts, baseline_fits):
            idx = np.flatnonzero(base & (np.abs(q["b"]) >= bcut))
            finite = np.ones(len(idx), dtype=bool)
            for name in variables:
                finite &= np.isfinite(q[name][idx])
            idx = idx[finite]
            z = q["z"][idx]
            vec = q["vec"][idx]
            base_dir = unit_from_lb(base_fit["l_deg"], base_fit["b_deg"])
            projection = vec @ base_dir
            hemisphere = projection >= 0.0
            if np.sum(hemisphere) < 100 or np.sum(~hemisphere) < 100:
                continue
            features = build_features(q, idx, variables, args.poly_degree)
            weights, propensity, info = overlap_weights(features, hemisphere, args.seed)
            fit = fit_weighted_dipole(z, vec, weights)
            direction_sep = angular_sep_deg(unit_from_lb(fit["l_deg"], fit["b_deg"]), base_dir)
            row = {
                "fit_group": group_name,
                "bcut_deg": bcut,
                "N": int(len(idx)),
                "ess": info["ess"],
                "variables": ",".join(variables),
                **fit,
                "amp_ratio_vs_baseline_same_bcut": fit["amp"] / base_fit["amp"],
                "direction_sep_vs_baseline_same_bcut_deg": direction_sep,
                **info,
            }
            fit_rows.append(row)
            weighted_fits.append(row)
            bal = balance_rows(q, idx, variables, hemisphere, weights, group_name, bcut)
            balance.extend(bal)
            group_state[group_name]["bcuts"][bcut] = {
                "idx": idx,
                "weights": weights,
                "propensity": propensity,
                "baseline_dir": base_dir,
                "base_fit": base_fit,
                "strata": propensity_strata(propensity, args.n_strata),
            }
        summary = summarize_fits(group_name, weighted_fits, baseline_fits)
        group_balance = [row for row in balance if row["fit_group"] == group_name]
        summary["max_abs_smd_before"] = float(np.nanmax([row["abs_smd_before"] for row in group_balance]))
        summary["max_abs_smd_after"] = float(np.nanmax([row["abs_smd_after"] for row in group_balance]))
        summary["max_weighted_ks_before"] = float(np.nanmax([row["weighted_ks_before"] for row in group_balance]))
        summary["max_weighted_ks_after"] = float(np.nanmax([row["weighted_ks_after"] for row in group_balance]))
        summary["passes_amp_dir_gate"] = bool(summary["amp_ratio_vs_baseline"] >= args.pass_amp_ratio and summary["max_direction_sep_vs_baseline_deg"] <= args.pass_direction_deg)
        summary["passes_balance_gate"] = bool(summary["max_abs_smd_after"] <= args.pass_smd and summary["max_weighted_ks_after"] <= args.pass_ks)
        summary["readout"] = (
            "survives matched-quality gate"
            if summary["passes_amp_dir_gate"] and summary["passes_balance_gate"]
            else "does not pass matched-quality gate"
        )
        summary_rows.append(summary)

    # Matched-quality redshift shuffles. Keep the fitted weights fixed and shuffle redshift
    # within propensity strata for each bcut/group. This preserves the quality-overlap
    # structure used by the weights.
    for group in quality_groups():
        group_name = group["group"]
        observed_summary = next(row for row in summary_rows if row["fit_group"] == group_name)
        state = group_state[group_name]["bcuts"]
        if args.n_mocks <= 0:
            continue
        for mock_id in range(args.n_mocks):
            fits = []
            for bcut in bcuts:
                st = state.get(bcut)
                if st is None:
                    continue
                idx = st["idx"]
                z_mock = q["z"][idx].copy()
                strata = st["strata"]
                for sid in np.unique(strata):
                    sub = np.flatnonzero(strata == sid)
                    if len(sub) > 1:
                        z_mock[sub] = rng.permutation(z_mock[sub])
                fit = fit_weighted_dipole(z_mock, q["vec"][idx], st["weights"])
                fits.append({"fit_group": group_name, "bcut_deg": bcut, "N": int(len(idx)), "ess": float((np.sum(st["weights"]) ** 2) / np.sum(st["weights"] * st["weights"])), **fit})
            if len(fits) == len(bcuts):
                summary = summarize_fits(group_name, fits, baseline_fits)
                summary["mock_id"] = mock_id
                mock_rows.append(summary)
            if args.progress_every and ((mock_id + 1) % args.progress_every == 0 or mock_id == 0):
                print(f"{group_name} matched-quality shuffle {mock_id + 1}/{args.n_mocks}", flush=True)
        sub = [row for row in mock_rows if row["fit_group"] == group_name]
        amp = np.array([row["amp_mean"] for row in sub], dtype=float)
        snr = np.array([row["snr_max"] for row in sub], dtype=float)
        sep = np.array([row["max_pair_direction_sep_deg"] for row in sub], dtype=float)
        coh = np.array([row["coherence_score"] for row in sub], dtype=float)
        joint = float((1 + np.sum((snr >= observed_summary["snr_max"]) & (sep <= observed_summary["max_pair_direction_sep_deg"]))) / (1 + len(sub))) if sub else float("nan")
        pvalue_rows.append(
            {
                "fit_group": group_name,
                "null_kind": "matched_quality_propensity_strata_zshuffle",
                "n_mocks": len(sub),
                "observed_amp_mean": observed_summary["amp_mean"],
                "mock_amp_mean_p50": qtile(amp, 0.50),
                "mock_amp_mean_p95": qtile(amp, 0.95),
                "p_amp_mean_ge_observed": p_ge(amp, observed_summary["amp_mean"]),
                "observed_snr_max": observed_summary["snr_max"],
                "mock_snr_max_p50": qtile(snr, 0.50),
                "mock_snr_max_p95": qtile(snr, 0.95),
                "mock_snr_max_p99": qtile(snr, 0.99),
                "p_snr_max_ge_observed": p_ge(snr, observed_summary["snr_max"]),
                "observed_max_pair_direction_sep_deg": observed_summary["max_pair_direction_sep_deg"],
                "mock_pair_sep_p05": qtile(sep, 0.05),
                "mock_pair_sep_p50": qtile(sep, 0.50),
                "p_pair_sep_le_observed": p_le(sep, observed_summary["max_pair_direction_sep_deg"]),
                "observed_coherence_score": observed_summary["coherence_score"],
                "mock_coherence_p50": qtile(coh, 0.50),
                "mock_coherence_p95": qtile(coh, 0.95),
                "p_coherence_ge_observed": p_ge(coh, observed_summary["coherence_score"]),
                "p_joint_snr_and_pair_sep": joint,
            }
        )

    config = {
        "date": "2026-07-15",
        "target": "broad_primary_0p95_1p45",
        "z_min": args.z_min,
        "z_max": args.z_max,
        "bcuts": bcuts,
        "quality_groups": quality_groups(),
        "method": "propensity-overlap weights balancing quality variables between hemispheres defined by each bcut baseline dipole axis",
        "equation_style": [
            "z_i      = c0 + d . n_i + epsilon_i",
            "A_dipole = ||d||",
            "SNR_d    = sqrt(d^T Cov(d)^(-1) d)",
        ],
        "n_mocks": args.n_mocks,
        "n_strata": args.n_strata,
        "seed": args.seed,
        "pass_amp_ratio": args.pass_amp_ratio,
        "pass_direction_deg": args.pass_direction_deg,
        "pass_smd": args.pass_smd,
        "pass_ks": args.pass_ks,
    }
    return fit_rows, summary_rows, balance, mock_rows, pvalue_rows, config


def write_report(path: Path, summaries: list[dict[str, Any]], pvalues: list[dict[str, Any]], args: argparse.Namespace) -> None:
    nonbase = [row for row in summaries if row["fit_group"] != "baseline_unweighted"]
    lines = [
        "# SU2 / Quaia Matched Catalogue-Quality Gate",
        "",
        "Date: July 15, 2026",
        "",
        "## Purpose",
        "",
        "This gate tests whether the locked Quaia `0.95 < z < 1.45` angular mode survives when catalogue-quality variables are explicitly balanced between the two hemispheres defined by the fitted dipole axis.",
        "",
        "The analysis uses propensity-overlap weights. For each latitude cut, a logistic model predicts the positive dipole-axis hemisphere from catalogue-quality variables. Objects are then overlap-weighted so the two hemispheres have matched quality distributions before the redshift dipole is refit.",
        "",
        "## Model",
        "",
        "```text",
        "z_i       = c0 + d . n_i + epsilon_i",
        "A_dipole  = ||d||",
        "SNR_d     = sqrt(d^T Cov(d)^(-1) d)",
        "```",
        "",
        "## Locked Configuration",
        "",
        f"- redshift window: `{args.z_min} <= z < {args.z_max}`",
        f"- latitude cuts: `{', '.join(str(x) for x in args.bcuts)}`",
        f"- matched-quality shuffle mocks per group: `{args.n_mocks}`",
        "",
        "## Matched-Quality Summary",
        "",
        markdown_table(
            nonbase,
            [
                "fit_group",
                "amp_ratio_vs_baseline",
                "max_direction_sep_vs_baseline_deg",
                "snr_max",
                "max_abs_smd_after",
                "max_weighted_ks_after",
                "readout",
            ],
        ),
        "",
        "## Matched-Quality Shuffle P-Values",
        "",
        markdown_table(
            pvalues,
            [
                "fit_group",
                "n_mocks",
                "p_snr_max_ge_observed",
                "p_pair_sep_le_observed",
                "p_coherence_ge_observed",
                "p_joint_snr_and_pair_sep",
            ],
        ),
        "",
        "## Interpretation Gate",
        "",
        f"- Amp/direction gate: amplitude ratio >= `{args.pass_amp_ratio}` and maximum direction rotation <= `{args.pass_direction_deg}` degrees.",
        f"- Balance gate: post-weight max absolute SMD <= `{args.pass_smd}` and post-weight weighted KS <= `{args.pass_ks}`.",
        "- If all-catalogue quality passes both gates, the catalogue-quality explanation is weakened.",
        "- If all-catalogue quality fails the amp/direction gate after good balance, the catalogue-quality explanation remains the leading caveat.",
        "",
        "## Outputs",
        "",
        f"- summary: `{args.out / 'su2_quaia_matched_quality_summary.csv'}`",
        f"- fit rows: `{args.out / 'su2_quaia_matched_quality_fit_rows.csv'}`",
        f"- balance diagnostics: `{args.out / 'su2_quaia_matched_quality_balance.csv'}`",
        f"- p-values: `{args.out / 'su2_quaia_matched_quality_pvalues.csv'}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_readout(path: Path, summaries: list[dict[str, Any]], pvalues: list[dict[str, Any]]) -> None:
    nonbase = [row for row in summaries if row["fit_group"] != "baseline_unweighted"]
    all_quality = next((row for row in nonbase if row["fit_group"] == "all_catalogue_quality"), None)
    all_p = next((row for row in pvalues if row["fit_group"] == "all_catalogue_quality"), None)
    if all_quality is None:
        bottom = "All-catalogue-quality gate did not run."
    elif all_quality.get("readout") == "survives matched-quality gate":
        bottom = "The locked mode survives the all-catalogue-quality matching gate; catalogue-quality alone is weakened as the explanation."
    else:
        bottom = "The locked mode does not pass the all-catalogue-quality matching gate; catalogue-quality remains the leading promotion-blocking caveat."
    lines = [
        "# SU2 / Quaia Matched Catalogue-Quality Gate Readout",
        "",
        "Date: July 15, 2026",
        "",
        "## Bottom Line",
        "",
        bottom,
        "",
        "## Gate Summary",
        "",
        markdown_table(
            nonbase,
            [
                "fit_group",
                "amp_ratio_vs_baseline",
                "max_direction_sep_vs_baseline_deg",
                "snr_max",
                "max_abs_smd_after",
                "max_weighted_ks_after",
                "readout",
            ],
        ),
        "",
        "## Shuffle P-Values",
        "",
        markdown_table(
            pvalues,
            [
                "fit_group",
                "n_mocks",
                "p_snr_max_ge_observed",
                "p_pair_sep_le_observed",
                "p_coherence_ge_observed",
                "p_joint_snr_and_pair_sep",
            ],
        ),
    ]
    if all_quality is not None and all_p is not None:
        lines.extend(
            [
                "",
                "## All-Catalogue Quality Detail",
                "",
                f"- amplitude ratio vs baseline: `{fmt(all_quality['amp_ratio_vs_baseline'])}`",
                f"- maximum direction rotation vs baseline: `{fmt(all_quality['max_direction_sep_vs_baseline_deg'])}` deg",
                f"- max post-weight absolute SMD: `{fmt(all_quality['max_abs_smd_after'])}`",
                f"- max post-weight weighted KS: `{fmt(all_quality['max_weighted_ks_after'])}`",
                f"- matched-quality joint p-value: `{fmt(all_p['p_joint_snr_and_pair_sep'])}`",
            ]
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Matched catalogue-quality gate for the locked SU2/Quaia angular mode.")
    parser.add_argument("--quaia-fits", type=Path, default=DEFAULT_QUAIA)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--z-min", type=float, default=0.95)
    parser.add_argument("--z-max", type=float, default=1.45)
    parser.add_argument("--bcuts", type=float, nargs="+", default=[10, 15, 20, 25, 30, 35])
    parser.add_argument("--n-mocks", type=int, default=500)
    parser.add_argument("--n-strata", type=int, default=10)
    parser.add_argument("--seed", type=int, default=170715)
    parser.add_argument("--poly-degree", type=int, default=2)
    parser.add_argument("--progress-every", type=int, default=100)
    parser.add_argument("--pass-amp-ratio", type=float, default=0.8)
    parser.add_argument("--pass-direction-deg", type=float, default=30.0)
    parser.add_argument("--pass-smd", type=float, default=0.1)
    parser.add_argument("--pass-ks", type=float, default=0.1)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    fit_rows, summaries, balance, mocks, pvalues, config = run(args)
    write_csv(args.out / "su2_quaia_matched_quality_fit_rows.csv", fit_rows)
    write_csv(args.out / "su2_quaia_matched_quality_summary.csv", summaries)
    write_csv(args.out / "su2_quaia_matched_quality_balance.csv", balance)
    write_csv(args.out / "su2_quaia_matched_quality_mock_summary.csv", mocks)
    write_csv(args.out / "su2_quaia_matched_quality_pvalues.csv", pvalues)
    (args.out / "su2_quaia_matched_quality_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")
    write_report(args.out / "su2_quaia_matched_quality_report.md", summaries, pvalues, args)
    write_readout(args.out / "su2_quaia_matched_quality_readout.md", summaries, pvalues)
    print(f"Saved summary: {args.out / 'su2_quaia_matched_quality_summary.csv'}", flush=True)
    print(f"Saved p-values: {args.out / 'su2_quaia_matched_quality_pvalues.csv'}", flush=True)
    print(f"Saved report: {args.out / 'su2_quaia_matched_quality_report.md'}", flush=True)
    print(f"Saved readout: {args.out / 'su2_quaia_matched_quality_readout.md'}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
