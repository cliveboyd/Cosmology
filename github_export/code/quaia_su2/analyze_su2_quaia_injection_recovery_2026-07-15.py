from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits
from scipy.stats import chi2
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(r"C:\Users\clive\Documents\Cosmology")
CODE_DIR = ROOT / "github_export" / "code" / "quaia_su2"
OUTPUTS = Path(r"C:\Users\clive\Documents\Codex\2026-07-14\ok-continu\outputs")
DEFAULT_OUT = ROOT / "plamb_runs" / "diagnostics" / "su2_quaia_injection_recovery_20260715"
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


def parse_csv_floats(text: str) -> list[float]:
    return [float(part.strip()) for part in text.split(",") if part.strip()]


def finite_mask(*arrays: np.ndarray) -> np.ndarray:
    mask = np.ones(len(arrays[0]), dtype=bool)
    for arr in arrays:
        mask &= np.isfinite(arr)
    return mask


def standardize_from_mask(values: np.ndarray, mask: np.ndarray) -> np.ndarray:
    med = float(np.nanmedian(values[mask]))
    scale = float(np.nanstd(values[mask]))
    if not math.isfinite(scale) or scale <= 0.0:
        scale = 1.0
    return (values - med) / scale


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


def random_unit(rng: np.random.Generator) -> np.ndarray:
    z = rng.uniform(-1.0, 1.0)
    phi = rng.uniform(0.0, 2.0 * math.pi)
    r = math.sqrt(max(1.0 - z * z, 0.0))
    return np.array([r * math.cos(phi), r * math.sin(phi), z], dtype=float)


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
    q["bp_rp_z"] = standardize_from_mask(q["bp_rp"], base_mask)
    q["w1_w2_z"] = standardize_from_mask(q["w1_w2"], base_mask)
    q["bp_rp_x_w1_w2_z"] = standardize_from_mask(q["bp_rp_z"] * q["w1_w2_z"], base_mask)
    meta = {
        "gaiascanlaw_data_dir": str(gaiascanlaw_data_dir) if gaiascanlaw_data_dir else None,
        "scanlaw_templates_present": bool(np.all(np.isfinite(q["gaia_scan_count_log1p_dr3"][base_mask]))),
    }
    return q, meta


def all_quality_variables(q: dict[str, np.ndarray]) -> list[tuple[str, np.ndarray]]:
    return [
        ("zerr", q["zerr"]),
        ("g", q["g"]),
        ("bp_rp", q["bp_rp"]),
        ("w1", q["w1"]),
        ("w2", q["w2"]),
        ("w1_w2", q["w1_w2"]),
        ("pmra_error", q["pmra_error"]),
        ("pmdec_error", q["pmdec_error"]),
    ]


def scanlaw_colour_controls(q: dict[str, np.ndarray]) -> list[tuple[str, np.ndarray]]:
    return [
        ("gaia_scan_count_log1p_dr3", q["gaia_scan_count_log1p_dr3"]),
        ("gaia_scan_angle_cos2_mean_dr3", q["gaia_scan_angle_cos2_mean_dr3"]),
        ("gaia_scan_angle_sin2_mean_dr3", q["gaia_scan_angle_sin2_mean_dr3"]),
        ("gaia_scan_angle_resultant_dr3", q["gaia_scan_angle_resultant_dr3"]),
        ("bp_rp_z", q["bp_rp_z"]),
        ("w1_w2_z", q["w1_w2_z"]),
        ("bp_rp_x_w1_w2_z", q["bp_rp_x_w1_w2_z"]),
    ]


def design_matrix(vec: np.ndarray, controls: list[np.ndarray] | None = None, include_dipole: bool = True) -> tuple[np.ndarray, slice | None]:
    columns = [np.ones(len(vec))]
    dip_slice: slice | None = None
    if include_dipole:
        start = len(columns)
        columns.extend([vec[:, 0], vec[:, 1], vec[:, 2]])
        dip_slice = slice(start, start + 3)
    if controls:
        for control in controls:
            c = np.asarray(control, dtype=float)
            med = float(np.nanmedian(c))
            scale = float(np.nanstd(c))
            if not math.isfinite(scale) or scale <= 0.0:
                scale = 1.0
            columns.append((c - med) / scale)
    return np.column_stack(columns).astype(float), dip_slice


class LinearCache:
    def __init__(self, x: np.ndarray, weights: np.ndarray | None = None, dip_slice: slice | None = None, bic_n: int | None = None):
        self.x = x
        self.n, self.k = x.shape
        self.dip_slice = dip_slice
        self.weights = None if weights is None else np.asarray(weights, dtype=float)
        if self.weights is None:
            self.xtw = x.T
            self.pinv = np.linalg.pinv(x.T @ x)
            self.ess = float(self.n)
        else:
            w = np.where(np.isfinite(self.weights) & (self.weights > 0.0), self.weights, 0.0)
            w = w / np.mean(w[w > 0.0])
            self.weights = w
            self.xtw = x.T * w
            self.pinv = np.linalg.pinv(x.T @ (x * w[:, None]))
            self.ess = float((np.sum(w) ** 2) / np.sum(w * w))
        self.bic_n = int(bic_n or self.n)
        if dip_slice is not None:
            self.dip_metric = np.linalg.pinv(self.pinv[dip_slice, dip_slice])
        else:
            self.dip_metric = None

    def fit_batch(self, z: np.ndarray, use_ess_dof: bool = False) -> dict[str, np.ndarray]:
        if z.ndim == 1:
            z = z[:, None]
        xty = self.xtw @ z
        beta = self.pinv @ xty
        if self.weights is None:
            ztz = np.sum(z * z, axis=0)
        else:
            ztz = np.sum(self.weights[:, None] * z * z, axis=0)
        rss = np.maximum(ztz - np.sum(beta * xty, axis=0), 0.0)
        dof_base = self.ess if use_ess_dof else self.n
        dof = np.maximum(dof_base - self.k, 1.0)
        sigma2 = rss / dof
        sigma2_mle = np.maximum(rss / max(self.bic_n, 1), 1e-300)
        loglike = -0.5 * self.bic_n * (np.log(2.0 * math.pi * sigma2_mle) + 1.0)
        out: dict[str, np.ndarray] = {
            "beta": beta,
            "rss": rss,
            "sigma2": sigma2,
            "aic": 2.0 * self.k - 2.0 * loglike,
            "bic": math.log(max(self.bic_n, 2)) * self.k - 2.0 * loglike,
            "loglike": loglike,
        }
        if self.dip_slice is not None and self.dip_metric is not None:
            d = beta[self.dip_slice, :]
            amp = np.sqrt(np.sum(d * d, axis=0))
            snr2 = np.einsum("ib,ij,jb->b", d, self.dip_metric, d) / np.maximum(sigma2, 1e-300)
            out["amp"] = amp
            out["snr"] = np.sqrt(np.maximum(snr2, 0.0))
            out["dvec"] = d
        return out


def fit_single_dipole(z: np.ndarray, vec: np.ndarray) -> dict[str, float]:
    x, dip_slice = design_matrix(vec, include_dipole=True)
    fit = LinearCache(x, dip_slice=dip_slice).fit_batch(z)
    dvec = fit["dvec"][:, 0]
    l_deg, b_deg = unit_to_lb(dvec)
    return {"amp": float(fit["amp"][0]), "snr": float(fit["snr"][0]), "l_deg": l_deg, "b_deg": b_deg, "dvec": dvec}


def overlap_weights(features: np.ndarray, hemisphere: np.ndarray, seed: int) -> tuple[np.ndarray, np.ndarray, float]:
    model = make_pipeline(
        StandardScaler(),
        LogisticRegression(max_iter=600, C=1.0, solver="lbfgs", random_state=seed),
    )
    model.fit(features, hemisphere.astype(int))
    propensity = np.clip(model.predict_proba(features)[:, 1], 0.02, 0.98)
    weights = np.where(hemisphere, 1.0 - propensity, propensity)
    weights /= np.mean(weights)
    ess = float((np.sum(weights) ** 2) / np.sum(weights * weights))
    return weights, propensity, ess


def propensity_strata(propensity: np.ndarray, n_strata: int) -> np.ndarray:
    edges = np.unique(np.quantile(propensity[np.isfinite(propensity)], np.linspace(0.0, 1.0, n_strata + 1)))
    if len(edges) <= 2:
        return np.zeros(len(propensity), dtype=np.int64)
    return np.searchsorted(edges[1:-1], propensity, side="right").astype(np.int64)


def make_batch_z(z_base: np.ndarray, strata: np.ndarray, batch_size: int, rng: np.random.Generator) -> np.ndarray:
    z = np.empty((len(z_base), batch_size), dtype=float)
    unique_strata = np.unique(strata)
    for col in range(batch_size):
        for sid in unique_strata:
            sub = np.flatnonzero(strata == sid)
            if len(sub) > 1:
                z[sub, col] = rng.permutation(z_base[sub])
            else:
                z[sub, col] = z_base[sub]
    return z


def direction_catalogue(q: dict[str, np.ndarray], args: argparse.Namespace, rng: np.random.Generator) -> tuple[list[dict[str, Any]], float, dict[str, Any]]:
    raw_fits: dict[float, dict[str, float]] = {}
    for bcut in sorted(set([args.primary_bcut, 15.0, 25.0])):
        mask = (q["z"] >= args.z_min) & (q["z"] < args.z_max) & (np.abs(q["b"]) >= bcut)
        raw_fits[bcut] = fit_single_dipole(q["z"][mask], q["vec"][mask])
    ref_amp = raw_fits[args.primary_bcut]["amp"]
    directions: list[dict[str, Any]] = [
        {"direction_label": f"locked_b{int(args.primary_bcut)}", "unit": raw_fits[args.primary_bcut]["dvec"] / raw_fits[args.primary_bcut]["amp"]},
        {"direction_label": "locked_b25", "unit": raw_fits[25.0]["dvec"] / raw_fits[25.0]["amp"]},
        {"direction_label": "cmb", "unit": unit_from_lb(264.021, 48.253)},
        {"direction_label": "anti_cmb", "unit": -unit_from_lb(264.021, 48.253)},
    ]
    for i in range(args.n_random_directions):
        directions.append({"direction_label": f"random_{i + 1:02d}", "unit": random_unit(rng)})
    for item in directions:
        l_deg, b_deg = unit_to_lb(item["unit"])
        item["l_deg"] = l_deg
        item["b_deg"] = b_deg
    meta = {
        "reference_bcut": args.primary_bcut,
        "reference_amp": ref_amp,
        "raw_reference_fits": {
            str(k): {kk: vv for kk, vv in v.items() if kk != "dvec"} for k, v in raw_fits.items()
        },
    }
    return directions, ref_amp, meta


def recovery_pass(amp_ratio: float, sep: float, snr: float, amp_scale: float, snr_cut: float, sep_cut: float) -> bool:
    if amp_scale <= 0.0:
        return bool(snr < snr_cut)
    return bool(math.isfinite(amp_ratio) and 0.5 <= amp_ratio <= 1.5 and sep <= sep_cut and snr >= snr_cut)


def run(args: argparse.Namespace) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    rng = np.random.default_rng(args.seed)
    q, load_meta = load_data(args)
    directions, reference_amp, direction_meta = direction_catalogue(q, args, rng)
    amp_scales = parse_csv_floats(args.amplitude_scales)
    quality_names, quality_values = zip(*all_quality_variables(q))
    control_names, control_values = zip(*scanlaw_colour_controls(q))
    trial_rows: list[dict[str, Any]] = []
    setup_rows: list[dict[str, Any]] = []
    started = time.time()

    for bcut in args.b_cuts:
        base_mask = (
            (q["z"] >= args.z_min)
            & (q["z"] < args.z_max)
            & (np.abs(q["b"]) >= bcut)
            & finite_mask(q["z"], q["vec"][:, 0], q["vec"][:, 1], q["vec"][:, 2])
        )
        idx0 = np.flatnonzero(base_mask)
        finite = np.ones(len(idx0), dtype=bool)
        for values in list(quality_values) + list(control_values):
            finite &= np.isfinite(values[idx0])
        idx = idx0[finite]
        z_base = q["z"][idx].astype(float)
        vec = q["vec"][idx].astype(float)
        raw_x, raw_dip = design_matrix(vec, include_dipole=True)
        raw_cache = LinearCache(raw_x, dip_slice=raw_dip)
        controls = [values[idx] for values in control_values]
        restricted_x, _ = design_matrix(vec, controls=controls, include_dipole=False)
        full_x, full_dip = design_matrix(vec, controls=controls, include_dipole=True)
        restricted_cache = LinearCache(restricted_x, dip_slice=None, bic_n=len(idx))
        full_cache = LinearCache(full_x, dip_slice=full_dip, bic_n=len(idx))
        quality_features = np.column_stack([values[idx] for values in quality_values]).astype(float)
        print(f"bcut={bcut:g} N={len(idx)} controls={len(controls)}", flush=True)

        for direction in directions:
            dunit = direction["unit"]
            projection = vec @ dunit
            hemisphere = projection >= 0.0
            if np.sum(hemisphere) < args.min_hemisphere or np.sum(~hemisphere) < args.min_hemisphere:
                print(f"skip direction {direction['direction_label']} bcut={bcut:g}: small hemisphere", flush=True)
                continue
            weights, propensity, ess = overlap_weights(quality_features, hemisphere, args.seed + int(100 * bcut) + len(setup_rows))
            strata = propensity_strata(propensity, args.n_strata)
            matched_cache = LinearCache(raw_x, weights=weights, dip_slice=raw_dip)
            setup_rows.append(
                {
                    "direction_label": direction["direction_label"],
                    "direction_l_deg": direction["l_deg"],
                    "direction_b_deg": direction["b_deg"],
                    "bcut_deg": bcut,
                    "N": int(len(idx)),
                    "matched_ess": ess,
                    "propensity_min": float(np.min(propensity)),
                    "propensity_p50": float(np.quantile(propensity, 0.50)),
                    "propensity_max": float(np.max(propensity)),
                }
            )
            for amp_scale in amp_scales:
                injected_amp = reference_amp * amp_scale
                signal = injected_amp * projection
                n_done = 0
                while n_done < args.n_mocks:
                    batch = min(args.batch_size, args.n_mocks - n_done)
                    z_batch = make_batch_z(z_base, strata, batch, rng)
                    if injected_amp != 0.0:
                        z_batch += signal[:, None]
                    raw_fit = raw_cache.fit_batch(z_batch)
                    matched_fit = matched_cache.fit_batch(z_batch, use_ess_dof=True)
                    restricted_fit = restricted_cache.fit_batch(z_batch)
                    full_fit = full_cache.fit_batch(z_batch)
                    delta_bic = full_fit["bic"] - restricted_fit["bic"]
                    lrt = np.maximum(2.0 * (full_fit["loglike"] - restricted_fit["loglike"]), 0.0)
                    lrt_p = chi2.sf(lrt, 3)
                    partial_r2 = 1.0 - full_fit["rss"] / np.maximum(restricted_fit["rss"], 1e-300)
                    for j in range(batch):
                        mock_id = n_done + j
                        raw_d = raw_fit["dvec"][:, j]
                        matched_d = matched_fit["dvec"][:, j]
                        full_d = full_fit["dvec"][:, j]
                        raw_sep = angular_sep_deg(raw_d, dunit)
                        matched_sep = angular_sep_deg(matched_d, dunit)
                        full_sep = angular_sep_deg(full_d, dunit)
                        raw_amp_ratio = float(raw_fit["amp"][j] / injected_amp) if injected_amp > 0.0 else float("nan")
                        matched_amp_ratio = float(matched_fit["amp"][j] / injected_amp) if injected_amp > 0.0 else float("nan")
                        full_amp_ratio = float(full_fit["amp"][j] / injected_amp) if injected_amp > 0.0 else float("nan")
                        trial_rows.append(
                            {
                                "direction_label": direction["direction_label"],
                                "direction_l_deg": direction["l_deg"],
                                "direction_b_deg": direction["b_deg"],
                                "amp_scale": amp_scale,
                                "injected_amp": injected_amp,
                                "bcut_deg": bcut,
                                "mock_id": mock_id,
                                "raw_amp": float(raw_fit["amp"][j]),
                                "raw_snr": float(raw_fit["snr"][j]),
                                "raw_sep_from_injected_deg": raw_sep,
                                "raw_amp_ratio": raw_amp_ratio,
                                "raw_pass": recovery_pass(raw_amp_ratio, raw_sep, float(raw_fit["snr"][j]), amp_scale, args.snr_cut, args.sep_cut_deg),
                                "matched_amp": float(matched_fit["amp"][j]),
                                "matched_snr": float(matched_fit["snr"][j]),
                                "matched_sep_from_injected_deg": matched_sep,
                                "matched_amp_ratio": matched_amp_ratio,
                                "matched_pass": recovery_pass(matched_amp_ratio, matched_sep, float(matched_fit["snr"][j]), amp_scale, args.snr_cut, args.sep_cut_deg),
                                "scanlaw_colour_dipole_amp": float(full_fit["amp"][j]),
                                "scanlaw_colour_dipole_snr": float(full_fit["snr"][j]),
                                "scanlaw_colour_sep_from_injected_deg": full_sep,
                                "scanlaw_colour_amp_ratio": full_amp_ratio,
                                "scanlaw_colour_delta_bic": float(delta_bic[j]),
                                "scanlaw_colour_lrt_p": float(lrt_p[j]),
                                "scanlaw_colour_partial_r2": float(partial_r2[j]),
                                "scanlaw_colour_pass": bool(
                                    recovery_pass(full_amp_ratio, full_sep, float(full_fit["snr"][j]), amp_scale, args.snr_cut, args.sep_cut_deg)
                                    and float(delta_bic[j]) < args.bic_pass_delta
                                ),
                            }
                        )
                    n_done += batch
                print(
                    f"completed bcut={bcut:g} dir={direction['direction_label']} amp_scale={amp_scale:g} "
                    f"elapsed={(time.time() - started) / 60.0:.1f} min",
                    flush=True,
                )

    summary_rows = summarise_trials(trial_rows)
    config = {
        "date": "2026-07-15",
        "analysis": "su2_quaia_injection_recovery",
        "z_min": args.z_min,
        "z_max": args.z_max,
        "b_cuts": args.b_cuts,
        "n_mocks": args.n_mocks,
        "batch_size": args.batch_size,
        "amplitude_scales": amp_scales,
        "n_random_directions": args.n_random_directions,
        "seed": args.seed,
        "reference_amp": reference_amp,
        "quality_variables": list(quality_names),
        "scanlaw_colour_controls": list(control_names),
        "external_inputs": load_meta,
        "direction_meta": direction_meta,
        "gate": "Injected dipoles should survive matched-quality weighting and scan-law+colour likelihood if these gates are not over-correcting real angular modes.",
    }
    return trial_rows, summary_rows, setup_rows, config


def summarise_trials(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[tuple[Any, ...], list[dict[str, Any]]] = {}
    for row in rows:
        key = (row["direction_label"], row["amp_scale"], row["bcut_deg"])
        groups.setdefault(key, []).append(row)
    summary: list[dict[str, Any]] = []
    for (direction_label, amp_scale, bcut), sub in sorted(groups.items(), key=lambda item: (item[0][0], item[0][1], item[0][2])):
        def arr(name: str) -> np.ndarray:
            return np.array([r[name] for r in sub], dtype=float)

        def rate(name: str) -> float:
            return float(np.mean([bool(r[name]) for r in sub]))

        summary.append(
            {
                "direction_label": direction_label,
                "amp_scale": amp_scale,
                "bcut_deg": bcut,
                "n_mocks": len(sub),
                "raw_pass_rate": rate("raw_pass"),
                "matched_pass_rate": rate("matched_pass"),
                "scanlaw_colour_pass_rate": rate("scanlaw_colour_pass"),
                "raw_snr_p50": float(np.quantile(arr("raw_snr"), 0.50)),
                "raw_snr_p95": float(np.quantile(arr("raw_snr"), 0.95)),
                "matched_snr_p50": float(np.quantile(arr("matched_snr"), 0.50)),
                "matched_snr_p95": float(np.quantile(arr("matched_snr"), 0.95)),
                "scanlaw_colour_snr_p50": float(np.quantile(arr("scanlaw_colour_dipole_snr"), 0.50)),
                "scanlaw_colour_snr_p95": float(np.quantile(arr("scanlaw_colour_dipole_snr"), 0.95)),
                "raw_sep_p50_deg": float(np.nanquantile(arr("raw_sep_from_injected_deg"), 0.50)),
                "matched_sep_p50_deg": float(np.nanquantile(arr("matched_sep_from_injected_deg"), 0.50)),
                "scanlaw_colour_sep_p50_deg": float(np.nanquantile(arr("scanlaw_colour_sep_from_injected_deg"), 0.50)),
                "raw_amp_ratio_p50": float(np.nanquantile(arr("raw_amp_ratio"), 0.50)) if float(amp_scale) > 0.0 else float("nan"),
                "matched_amp_ratio_p50": float(np.nanquantile(arr("matched_amp_ratio"), 0.50)) if float(amp_scale) > 0.0 else float("nan"),
                "scanlaw_colour_amp_ratio_p50": float(np.nanquantile(arr("scanlaw_colour_amp_ratio"), 0.50)) if float(amp_scale) > 0.0 else float("nan"),
                "scanlaw_colour_delta_bic_p50": float(np.quantile(arr("scanlaw_colour_delta_bic"), 0.50)),
                "scanlaw_colour_delta_bic_p05": float(np.quantile(arr("scanlaw_colour_delta_bic"), 0.05)),
                "scanlaw_colour_delta_bic_p95": float(np.quantile(arr("scanlaw_colour_delta_bic"), 0.95)),
            }
        )
    return summary


def write_reports(out_dir: Path, summary: list[dict[str, Any]], setup: list[dict[str, Any]], config: dict[str, Any]) -> None:
    positive = [row for row in summary if float(row["amp_scale"]) > 0.0]
    one_x = [row for row in positive if abs(float(row["amp_scale"]) - 1.0) < 1e-9]
    if one_x:
        matched_mean = float(np.mean([row["matched_pass_rate"] for row in one_x]))
        scanlaw_mean = float(np.mean([row["scanlaw_colour_pass_rate"] for row in one_x]))
        bottom = (
            f"At 1x injected amplitude, mean matched-quality recovery is {matched_mean:.3g} "
            f"and mean scan-law+colour likelihood recovery is {scanlaw_mean:.3g}."
        )
    else:
        bottom = "No 1x injected-amplitude rows were available in this run."
    cols = [
        "direction_label",
        "amp_scale",
        "bcut_deg",
        "n_mocks",
        "raw_pass_rate",
        "matched_pass_rate",
        "scanlaw_colour_pass_rate",
        "matched_snr_p50",
        "scanlaw_colour_delta_bic_p50",
    ]
    one_x_rows = [row for row in summary if abs(float(row["amp_scale"]) - 1.0) < 1e-9]
    zero_rows = [row for row in summary if abs(float(row["amp_scale"])) < 1e-12]
    report = [
        "# SU2 / Quaia Injection-Recovery Stress Test",
        "",
        "Date: July 15, 2026",
        "",
        "## Purpose",
        "",
        "This overnight stress test injects known redshift dipoles into the locked Quaia `0.95 <= z < 1.45` sample and reruns the main gates. It asks whether matched catalogue-quality weighting and Gaia scan-law plus colour likelihood controls would preserve a real angular mode or erase it.",
        "",
        "## Model",
        "",
        "```text",
        "z_mock       = shuffle_strata(z_observed) + A_inj (n_i . d_hat)",
        "A_inj        = scale * A_locked_reference",
        "raw gate     = z_mock ~ 1 + n_x + n_y + n_z",
        "matched gate = weighted(z_mock ~ 1 + n_x + n_y + n_z)",
        "template gate= compare z_mock ~ S_scan + C_colour with z_mock ~ n + S_scan + C_colour",
        "```",
        "",
        "## Bottom Line",
        "",
        bottom,
        "",
        "## 1x Injection Summary",
        "",
        markdown_table(one_x_rows, cols),
        "",
        "## Zero-Injection False-Positive Summary",
        "",
        markdown_table(zero_rows, cols),
        "",
        "## Full Summary",
        "",
        markdown_table(summary, cols),
        "",
        "## Setup Diagnostics",
        "",
        markdown_table(setup[:20], ["direction_label", "bcut_deg", "N", "matched_ess", "propensity_min", "propensity_p50", "propensity_max"]),
        "",
        "## Configuration",
        "",
        "```json",
        json.dumps(config, indent=2),
        "```",
    ]
    (out_dir / "su2_quaia_injection_recovery_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    readout = [
        "# SU2 / Quaia Injection-Recovery Stress Test Readout",
        "",
        "Date: July 15, 2026",
        "",
        "## Bottom Line",
        "",
        bottom,
        "",
        "## 1x Injection Summary",
        "",
        markdown_table(one_x_rows, cols),
    ]
    (out_dir / "su2_quaia_injection_recovery_readout.md").write_text("\n".join(readout) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="SU2/Quaia injection-recovery stress test.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--quaia", type=Path, default=EXT.DEFAULT_QUAIA)
    parser.add_argument("--randoms", type=Path, default=EXT.DEFAULT_RANDOMS)
    parser.add_argument("--selection", type=Path, default=EXT.DEFAULT_SELECTION)
    parser.add_argument("--sfd-dir", type=Path, default=EXT.DEFAULT_SFD_DIR)
    parser.add_argument("--gaiascanlaw-data-dir", type=Path, default=None)
    parser.add_argument("--z-min", type=float, default=0.95)
    parser.add_argument("--z-max", type=float, default=1.45)
    parser.add_argument("--b-cuts", type=float, nargs="+", default=[10.0, 15.0, 20.0, 25.0, 30.0, 35.0])
    parser.add_argument("--primary-bcut", type=float, default=15.0)
    parser.add_argument("--amplitude-scales", default="0,0.25,0.5,1.0,1.5,2.0")
    parser.add_argument("--n-random-directions", type=int, default=6)
    parser.add_argument("--n-mocks", type=int, default=500)
    parser.add_argument("--batch-size", type=int, default=20)
    parser.add_argument("--n-strata", type=int, default=10)
    parser.add_argument("--min-hemisphere", type=int, default=100)
    parser.add_argument("--snr-cut", type=float, default=3.0)
    parser.add_argument("--sep-cut-deg", type=float, default=30.0)
    parser.add_argument("--bic-pass-delta", type=float, default=-10.0)
    parser.add_argument("--seed", type=int, default=190715)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    trial_rows, summary_rows, setup_rows, config = run(args)
    write_csv(args.out_dir / "su2_quaia_injection_recovery_trials.csv", trial_rows)
    write_csv(args.out_dir / "su2_quaia_injection_recovery_summary.csv", summary_rows)
    write_csv(args.out_dir / "su2_quaia_injection_recovery_setup.csv", setup_rows)
    (args.out_dir / "su2_quaia_injection_recovery_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")
    write_reports(args.out_dir, summary_rows, setup_rows, config)
    print(f"Saved trials: {args.out_dir / 'su2_quaia_injection_recovery_trials.csv'}", flush=True)
    print(f"Saved summary: {args.out_dir / 'su2_quaia_injection_recovery_summary.csv'}", flush=True)
    print(f"Saved report: {args.out_dir / 'su2_quaia_injection_recovery_report.md'}", flush=True)
    print(f"Saved readout: {args.out_dir / 'su2_quaia_injection_recovery_readout.md'}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
