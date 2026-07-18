from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import importlib.util
import json
import math
import shutil
import sys
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits


ROOT = Path(r"C:\Users\clive\Documents\Cosmology")
CODE_DIR = ROOT / "github_export" / "code" / "quaia_su2"
RESULT_DIR = ROOT / "github_export" / "results" / "2026-07-18" / "su2_quaia"
DEFAULT_OUT = ROOT / "plamb_runs" / "diagnostics" / "su2_quaia_injection_calibrated_angular_likelihood_20260718"
PREFIX = "su2_quaia_injection_calibrated_angular_likelihood"
DATE_SUFFIX = "2026-07-18"

EXT_PATH = CODE_DIR / "analyze_su2_quaia_external_dust_gate_2026-07-15.py"
WISE_PATH = CODE_DIR / "analyze_su2_quaia_wise_stellar_gate_2026-07-15.py"
PREREG_MD = RESULT_DIR / f"{PREFIX}_preregistration_{DATE_SUFFIX}.md"
PREREG_JSON = RESULT_DIR / f"{PREFIX}_preregistration_{DATE_SUFFIX}.json"
LCDM_BESTFIT = ROOT / "plamb_runs" / "updated_datasets_20260710" / "runs" / "desi_LCDM_Pantheon_SN_BAO_Planck" / "LCDM_emcee_bestfit.txt"
SU2R_BESTFIT = ROOT / "plamb_runs" / "updated_datasets_20260710" / "runs" / "desi_SU2R_Pantheon_SN_BAO_Planck" / "SU2R_emcee_bestfit.txt"


def load_module(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {name} from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


EXT = load_module("analyze_su2_quaia_external_dust_gate", EXT_PATH)
sys.modules["analyze_su2_quaia_external_dust_gate"] = EXT
WISE = load_module("analyze_su2_quaia_wise_stellar_gate", WISE_PATH)
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from PLamb_Test_10V6c_plus import mu_model  # noqa: E402


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def clean_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): clean_json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean_json(item) for item in value]
    if isinstance(value, np.ndarray):
        return clean_json(value.tolist())
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        number = float(value)
        return number if math.isfinite(number) else None
    if isinstance(value, Path):
        return str(value)
    return value


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(clean_json(value), indent=2) + "\n", encoding="utf-8")


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
    if isinstance(value, (float, np.floating)):
        return f"{float(value):.6g}" if math.isfinite(float(value)) else ""
    return str(value)


def markdown_table(rows: list[dict[str, Any]], columns: list[str]) -> str:
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(fmt(row.get(column, "")) for column in columns) + " |")
    return "\n".join(lines)


def parse_csv_floats(text: str) -> list[float]:
    return [float(item.strip()) for item in text.split(",") if item.strip()]


def read_bestfit(path: Path) -> dict[str, Any]:
    values: dict[str, Any] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
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


def model_mu(z: np.ndarray, values: dict[str, Any], nint: int) -> np.ndarray:
    return np.asarray(
        mu_model(
            z,
            float(values.get("H0", 70.0)),
            float(values.get("Om", 0.3)),
            float(values.get("Ol", 0.7)),
            0.0,
            0.0,
            0.0,
            0.0,
            nint_base=nint,
            Omega_chi0=float(values.get("Omega_chi0", 0.0)),
            w0_chi=float(values.get("w0_chi", -1.0)),
            wa_chi=float(values.get("wa_chi", 0.0)),
        ),
        dtype=float,
    )


def unit_from_lb(l_deg: float, b_deg: float) -> np.ndarray:
    l_rad = math.radians(l_deg)
    b_rad = math.radians(b_deg)
    return np.array(
        [math.cos(b_rad) * math.cos(l_rad), math.cos(b_rad) * math.sin(l_rad), math.sin(b_rad)],
        dtype=float,
    )


def unit_to_lb(vector: np.ndarray) -> tuple[float, float]:
    norm = float(np.linalg.norm(vector))
    if not math.isfinite(norm) or norm <= 0.0:
        return float("nan"), float("nan")
    unit = vector / norm
    return math.degrees(math.atan2(float(unit[1]), float(unit[0]))) % 360.0, math.degrees(
        math.asin(float(np.clip(unit[2], -1.0, 1.0)))
    )


def vector_diagnostics(vectors: np.ndarray, target: np.ndarray, target_amp: float) -> dict[str, np.ndarray]:
    amplitudes = np.linalg.norm(vectors, axis=0)
    dots = np.einsum("ib,i->b", vectors, target)
    cosines = dots / np.maximum(amplitudes, 1e-300)
    separations = np.degrees(np.arccos(np.clip(cosines, -1.0, 1.0)))
    ratios = amplitudes / target_amp if target_amp > 0.0 else np.full_like(amplitudes, np.nan)
    return {
        "amplitude": amplitudes,
        "amplitude_ratio": ratios,
        "direction_separation_deg": separations,
        "positive_orientation": dots > 0.0,
    }


def standardise(values: np.ndarray, mask: np.ndarray) -> tuple[np.ndarray, float, float]:
    subset = np.asarray(values[mask], dtype=float)
    centre = float(np.nanmedian(subset))
    scale = float(np.nanstd(subset))
    if not math.isfinite(scale) or scale <= 1e-12:
        scale = 1.0
    return (np.asarray(values, dtype=float) - centre) / scale, centre, scale


def scalar_curve(args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    grid = np.linspace(args.z_min, args.z_max, args.scalar_grid_points)
    lcdm = read_bestfit(args.lcdm_bestfit)
    su2r = read_bestfit(args.su2r_bestfit)
    delta = model_mu(grid, su2r, args.integration_nint) - model_mu(grid, lcdm, args.integration_nint)
    if not np.all(np.isfinite(delta)) or float(np.std(delta)) <= 0.0:
        raise RuntimeError("The locked SU2R-minus-LCDM scalar curve is non-finite or constant.")
    return grid, delta, {
        "lcdm_bestfit": str(args.lcdm_bestfit),
        "lcdm_sha256": sha256_file(args.lcdm_bestfit),
        "su2r_bestfit": str(args.su2r_bestfit),
        "su2r_sha256": sha256_file(args.su2r_bestfit),
        "grid_points": args.scalar_grid_points,
        "integration_nint": args.integration_nint,
        "delta_mu_min": float(np.min(delta)),
        "delta_mu_max": float(np.max(delta)),
    }


def load_inputs(args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[str], list[dict[str, Any]], dict[str, Any], list[dict[str, Any]]]:
    scanlaw_dir = args.gaiascanlaw_data_dir or EXT.resolve_gaiascanlaw_data_dir()
    q = EXT.load_data(args.quaia, args.selection, args.randoms, args.sfd_dir, scanlaw_dir)
    with fits.open(args.quaia, memmap=True) as hdul:
        data = hdul[1].data
        q["bp_rp"] = np.asarray(data["phot_bp_mean_mag"], dtype=float) - np.asarray(data["phot_rp_mean_mag"], dtype=float)
    q.update(WISE.load_wise_depth(args.wise_inventory, q))

    gaia_density_source = "omitted by command-line request"
    if not args.skip_gaia_density:
        density, gaia_density_source = WISE.load_gaia_density(
            args.gaia_density,
            q,
            args.gaia_density_hips_base_url,
            args.gaia_density_hips_cache_dir,
            args.gaia_density_hips_order,
        )
        q.update(density)

    selection = (
        (q["z"] >= args.z_min)
        & (q["z"] < args.z_max)
        & (np.abs(q["b"]) >= args.b_cut)
        & np.isfinite(q["z"])
        & np.all(np.isfinite(q["vec"]), axis=1)
    )
    q["w1_w2_z"], _, _ = standardise(q["w1_w2"], selection)
    q["bp_rp_z"], _, _ = standardise(q["bp_rp"], selection)
    colour_cross = q["w1_w2_z"] * q["bp_rp_z"]
    q["bp_rp_x_w1_w2_z"], _, _ = standardise(colour_cross, selection)

    groups: list[tuple[str, list[str]]] = [
        ("dust", ["ebv_sfd", "sfd_log1p", "sfd_sq"]),
        (
            "depth",
            [
                "selection_T",
                "random_density_log1p",
                "wise_w1_medcov_log1p",
                "wise_w2_medcov_log1p",
                "wise_w1_lowcov_log1p",
                "wise_w2_lowcov_log1p",
                "gaia_dr3_density_hpx4_log1p",
            ],
        ),
        (
            "gaia_scanning_geometry",
            [
                "gaia_scan_count_log1p_dr3",
                "gaia_scan_angle_cos2_mean_dr3",
                "gaia_scan_angle_sin2_mean_dr3",
                "gaia_scan_angle_resultant_dr3",
            ],
        ),
        ("colour", ["w1_w2_z", "bp_rp_z", "bp_rp_x_w1_w2_z"]),
        ("catalogue_quality", ["zerr", "g", "w1", "w2", "pmra_error", "pmdec_error"]),
    ]

    inventory: list[dict[str, Any]] = []
    retained: list[tuple[str, str, np.ndarray]] = []
    for group, names in groups:
        for name in names:
            if name not in q:
                inventory.append({"group": group, "template": name, "status": "omitted_unavailable"})
                continue
            values = np.asarray(q[name], dtype=float)
            finite_count = int(np.sum(selection & np.isfinite(values)))
            scale = float(np.nanstd(values[selection & np.isfinite(values)])) if finite_count else float("nan")
            status = "retained" if finite_count > 0 and math.isfinite(scale) and scale > 1e-12 else "omitted_zero_or_nonfinite"
            inventory.append(
                {
                    "group": group,
                    "template": name,
                    "status": status,
                    "finite_rows_before_joint_mask": finite_count,
                    "raw_std_before_joint_mask": scale,
                }
            )
            if status == "retained":
                retained.append((group, name, values))

    joint_mask = selection.copy()
    for _, _, values in retained:
        joint_mask &= np.isfinite(values)
    if int(np.sum(joint_mask)) < 10000:
        raise RuntimeError(f"Only {int(np.sum(joint_mask))} rows remain after the locked joint finite mask.")

    template_columns: list[np.ndarray] = []
    template_names: list[str] = []
    for group, name, values in retained:
        transformed, centre, scale = standardise(values, joint_mask)
        template_columns.append(transformed[joint_mask])
        template_names.append(name)
        row = next(item for item in inventory if item["group"] == group and item["template"] == name)
        row.update({"joint_rows": int(np.sum(joint_mask)), "centre": centre, "scale": scale})

    grid, delta, curve_meta = scalar_curve(args)
    scalar = np.interp(q["z"][joint_mask], grid, delta)
    scalar_mean = float(np.mean(scalar))
    scalar_std = float(np.std(scalar))
    scalar_weight = (scalar - scalar_mean) / scalar_std
    vector = np.asarray(q["vec"][joint_mask], dtype=float)
    su2_basis = scalar_weight[:, None] * vector
    templates = np.column_stack(template_columns).astype(float)
    response = np.asarray(q["z"][joint_mask], dtype=float)
    l_deg = np.asarray(q["l"][joint_mask], dtype=float)
    b_deg = np.asarray(q["b"][joint_mask], dtype=float)

    curve_rows = [
        {"z": float(z_value), "delta_mu_su2r_minus_lcdm": float(delta_value)}
        for z_value, delta_value in zip(grid, delta)
    ]
    meta = {
        "source_rows": int(len(q["z"])),
        "window_latitude_rows": int(np.sum(selection)),
        "joint_finite_rows": int(np.sum(joint_mask)),
        "template_count": int(templates.shape[1]),
        "template_names": template_names,
        "gaiascanlaw_data_dir": str(scanlaw_dir) if scanlaw_dir else None,
        "gaia_density_source": gaia_density_source,
        "scalar_mean": scalar_mean,
        "scalar_std": scalar_std,
        "curve": curve_meta,
    }
    return response, templates, su2_basis, l_deg, b_deg, template_names, inventory, meta, curve_rows


def make_cells_and_folds(l_deg: np.ndarray, b_deg: np.ndarray, cell_deg: float, n_folds: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    n_lon = max(int(round(360.0 / cell_deg)), 4)
    n_sin_b = max(int(round(180.0 / cell_deg)), 2)
    lon_bin = np.floor((np.mod(l_deg, 360.0) / 360.0) * n_lon).astype(int)
    sin_b = np.sin(np.radians(b_deg))
    sin_bin = np.floor(((sin_b + 1.0) / 2.0) * n_sin_b).astype(int)
    lon_bin = np.clip(lon_bin, 0, n_lon - 1)
    sin_bin = np.clip(sin_bin, 0, n_sin_b - 1)
    raw = sin_bin * n_lon + lon_bin
    unique_raw, cell = np.unique(raw, return_inverse=True)
    occupied_lon_bin = unique_raw % n_lon
    cell_centre_lon = (occupied_lon_bin + 0.5) * 360.0 / n_lon
    cell_fold = np.floor(cell_centre_lon / (360.0 / n_folds)).astype(int)
    cell_fold = np.clip(cell_fold, 0, n_folds - 1)
    object_fold = cell_fold[cell]
    return cell.astype(np.int64), cell_fold.astype(np.int64), object_fold.astype(np.int64), {
        "cell_degrees": cell_deg,
        "occupied_cells": int(len(unique_raw)),
        "n_longitude_bins": n_lon,
        "n_sine_latitude_bins": n_sin_b,
        "n_folds": n_folds,
        "fold_rule": "12-degree cell centre assigned to fixed Galactic-longitude sector",
    }


def aggregate_x_residual(x: np.ndarray, residual: np.ndarray, cell: np.ndarray, n_cells: int) -> np.ndarray:
    out = np.empty((x.shape[1], n_cells), dtype=float)
    for column in range(x.shape[1]):
        out[column] = np.bincount(cell, weights=x[:, column] * residual, minlength=n_cells)
    return out


def build_model_stats(x: np.ndarray, base: np.ndarray, signal: np.ndarray, residual: np.ndarray, cell: np.ndarray, object_fold: np.ndarray, n_folds: int) -> dict[str, Any]:
    xtx_global = x.T @ x
    xtx_test = np.stack([x[object_fold == fold].T @ x[object_fold == fold] for fold in range(n_folds)])
    return {
        "p": int(x.shape[1]),
        "xtx_global": xtx_global,
        "xtx_test": xtx_test,
        "pinv_train": np.stack([np.linalg.pinv(xtx_global - xtx_test[fold], rcond=1e-12) for fold in range(n_folds)]),
        "base_xty_global": x.T @ base,
        "base_xty_test": np.stack([x[object_fold == fold].T @ base[object_fold == fold] for fold in range(n_folds)]),
        "signal_xty_global": x.T @ signal,
        "signal_xty_test": np.stack([x[object_fold == fold].T @ signal[object_fold == fold] for fold in range(n_folds)]),
        "cell_x_residual": aggregate_x_residual(x, residual, cell, int(np.max(cell)) + 1),
    }


def sum_by_fold(values: np.ndarray, object_fold: np.ndarray, n_folds: int) -> np.ndarray:
    return np.array([float(np.sum(values[object_fold == fold])) for fold in range(n_folds)], dtype=float)


def build_response_stats(base: np.ndarray, signal: np.ndarray, residual: np.ndarray, cell: np.ndarray, object_fold: np.ndarray, n_folds: int, cell_fold: np.ndarray) -> dict[str, Any]:
    components = {
        "base2": base * base,
        "base_signal": base * signal,
        "signal2": signal * signal,
        "residual2": residual * residual,
    }
    out: dict[str, Any] = {
        "n_global": int(len(base)),
        "n_test": np.bincount(object_fold, minlength=n_folds).astype(float),
        "cell_fold": cell_fold,
    }
    for name, values in components.items():
        out[f"{name}_global"] = float(np.sum(values))
        out[f"{name}_test"] = sum_by_fold(values, object_fold, n_folds)
    n_cells = int(np.max(cell)) + 1
    out["cell_base_residual"] = np.bincount(cell, weights=base * residual, minlength=n_cells)
    out["cell_signal_residual"] = np.bincount(cell, weights=signal * residual, minlength=n_cells)
    return out


def precompute_sign_terms(models: list[dict[str, Any]], response_stats: dict[str, Any], signs: np.ndarray, n_folds: int) -> dict[str, Any]:
    model_terms: list[dict[str, Any]] = []
    for model in models:
        global_term = model["cell_x_residual"] @ signs
        test_terms = []
        for fold in range(n_folds):
            mask = response_stats["cell_fold"] == fold
            test_terms.append(model["cell_x_residual"][:, mask] @ signs[mask])
        model_terms.append({"global": global_term, "test": np.stack(test_terms)})
    global_base_residual = response_stats["cell_base_residual"] @ signs
    global_signal_residual = response_stats["cell_signal_residual"] @ signs
    test_base_residual = []
    test_signal_residual = []
    for fold in range(n_folds):
        mask = response_stats["cell_fold"] == fold
        test_base_residual.append(response_stats["cell_base_residual"][mask] @ signs[mask])
        test_signal_residual.append(response_stats["cell_signal_residual"][mask] @ signs[mask])
    return {
        "models": model_terms,
        "base_residual_global": global_base_residual,
        "signal_residual_global": global_signal_residual,
        "base_residual_test": np.stack(test_base_residual),
        "signal_residual_test": np.stack(test_signal_residual),
    }


def evaluate_model(model: dict[str, Any], response_stats: dict[str, Any], sign_terms: dict[str, Any], model_index: int, amplitude_scale: float, n_folds: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    term = sign_terms["models"][model_index]
    n_batch = term["global"].shape[1]
    fold_lpd = np.empty((n_folds, n_batch), dtype=float)
    fold_sse = np.empty((n_folds, n_batch), dtype=float)
    fold_beta = np.empty((n_folds, model["p"], n_batch), dtype=float)

    fixed_global = (
        response_stats["base2_global"]
        + 2.0 * amplitude_scale * response_stats["base_signal_global"]
        + amplitude_scale**2 * response_stats["signal2_global"]
        + response_stats["residual2_global"]
    )
    yty_global = fixed_global + 2.0 * (
        sign_terms["base_residual_global"] + amplitude_scale * sign_terms["signal_residual_global"]
    )
    xty_global = (
        model["base_xty_global"][:, None]
        + amplitude_scale * model["signal_xty_global"][:, None]
        + term["global"]
    )

    for fold in range(n_folds):
        fixed_test = (
            response_stats["base2_test"][fold]
            + 2.0 * amplitude_scale * response_stats["base_signal_test"][fold]
            + amplitude_scale**2 * response_stats["signal2_test"][fold]
            + response_stats["residual2_test"][fold]
        )
        yty_test = fixed_test + 2.0 * (
            sign_terms["base_residual_test"][fold]
            + amplitude_scale * sign_terms["signal_residual_test"][fold]
        )
        xty_test = (
            model["base_xty_test"][fold, :, None]
            + amplitude_scale * model["signal_xty_test"][fold, :, None]
            + term["test"][fold]
        )
        xty_train = xty_global - xty_test
        yty_train = yty_global - yty_test
        beta = model["pinv_train"][fold] @ xty_train
        fold_beta[fold] = beta
        train_sse = np.maximum(yty_train - np.sum(beta * xty_train, axis=0), 1e-20)
        train_dof = max(int(response_stats["n_global"] - response_stats["n_test"][fold]) - model["p"], 1)
        sigma2 = np.maximum(train_sse / train_dof, 1e-20)
        test_sse = np.maximum(
            yty_test
            - 2.0 * np.sum(beta * xty_test, axis=0)
            + np.einsum("ib,ij,jb->b", beta, model["xtx_test"][fold], beta),
            0.0,
        )
        fold_sse[fold] = test_sse
        n_test = response_stats["n_test"][fold]
        fold_lpd[fold] = -0.5 * (n_test * np.log(2.0 * math.pi * sigma2) + test_sse / sigma2)
    return fold_lpd, fold_sse, fold_beta


def evaluate_pair(models: list[dict[str, Any]], response_stats: dict[str, Any], sign_terms: dict[str, Any], amplitude_scale: float, n_folds: int) -> dict[str, np.ndarray]:
    lpd0, sse0, _ = evaluate_model(models[0], response_stats, sign_terms, 0, amplitude_scale, n_folds)
    lpd1, sse1, beta1 = evaluate_model(models[1], response_stats, sign_terms, 1, amplitude_scale, n_folds)
    weights = response_stats["n_test"] / np.sum(response_stats["n_test"])
    pooled_vector = np.einsum("f,fib->ib", weights, beta1[:, -3:, :])
    return {
        "delta_lpd": np.sum(lpd1 - lpd0, axis=0),
        "fold_delta_lpd": lpd1 - lpd0,
        "fold_delta_sse": sse1 - sse0,
        "fold_candidate_beta": beta1[:, -3:, :],
        "pooled_vector": pooled_vector,
    }


def in_sample_diagnostics(x0: np.ndarray, x1: np.ndarray, response: np.ndarray) -> dict[str, Any]:
    rows: dict[str, Any] = {}
    for label, design in [("template_only", x0), ("template_plus_su2", x1)]:
        beta = np.linalg.pinv(design.T @ design, rcond=1e-12) @ (design.T @ response)
        residual = response - design @ beta
        rss = float(np.dot(residual, residual))
        n = int(len(response))
        k = int(design.shape[1])
        bic = n * math.log(max(rss / n, 1e-300)) + k * math.log(n)
        rows[label] = {"rss": rss, "n": n, "k": k, "bic": bic, "beta": beta}
    return {
        "template_only_bic": rows["template_only"]["bic"],
        "template_plus_su2_bic": rows["template_plus_su2"]["bic"],
        "delta_bic_plus_minus_template": rows["template_plus_su2"]["bic"] - rows["template_only"]["bic"],
        "candidate_vector": rows["template_plus_su2"]["beta"][-3:],
        "note": "Diagnostic only; the preregistered decision rule forbids substitution for held-out predictive gain.",
    }


def build_reports(out_dir: Path, decision: dict[str, Any], power_rows: list[dict[str, Any]], fold_rows: list[dict[str, Any]], template_rows: list[dict[str, Any]], config: dict[str, Any]) -> tuple[Path, Path]:
    observed = decision["observed"]
    calibration = decision["calibration"]
    conditions = decision["conditions"]
    condition_rows = [{"condition": key, "pass": value} for key, value in conditions.items()]
    power_columns = [
        "amplitude_scale",
        "partition",
        "n_replicates",
        "predictive_detection_rate",
        "joint_recovery_rate",
        "delta_cvlpd_median",
        "amplitude_ratio_median",
        "direction_separation_median_deg",
    ]
    fold_columns = ["fold", "N_test", "delta_lpd", "delta_sse", "recovered_amp", "recovered_l_deg", "recovered_b_deg"]
    report = [
        "# SU2 / Quaia Injection-Calibrated Angular Likelihood",
        "",
        "Date: 18 July 2026",
        "",
        "## Bottom Line",
        "",
        f"The preregistered conjunctive outcome is **{decision['outcome']}**. {decision['interpretation']}",
        "",
        "## Locked Models",
        "",
        "```text",
        "M_T: z_i = alpha + T_i beta + epsilon_i",
        "M_S: z_i = alpha + T_i beta + d . [standardise(Delta mu_SU2R-LCDM(z_i)) n_i] + epsilon_i",
        "```",
        "",
        "`T_i` jointly contains SFD dust, Quaia and WISE depth, Gaia source-density and scanning-geometry fields, WISE/Gaia colour and catalogue-quality variables. Twelve fixed Galactic-longitude sky sectors are held out in turn.",
        "",
        "## Frozen Calibration",
        "",
        f"- threshold-calibration null replicates: `{calibration['threshold_replicates']}`",
        f"- untouched validation null replicates: `{calibration['validation_replicates']}`",
        f"- frozen 99th-percentile Delta CVLPD threshold: `{fmt(calibration['frozen_threshold'])}`",
        f"- validation-null false-positive rate: `{fmt(calibration['validation_false_positive_rate'])}`",
        f"- joint held-out/direction/amplitude recovery at 1x: `{fmt(calibration['joint_recovery_rate_1x'])}`",
        "",
        "## Observed Held-Out Result",
        "",
        f"- pooled Delta CVLPD: `{fmt(observed['delta_cvlpd'])}`",
        f"- empirical validation-null p-value: `{fmt(observed['empirical_validation_null_p'])}`",
        f"- positive held-out sectors: `{observed['positive_fold_count']}/{config['held_out']['n_folds']}`",
        f"- pooled recovered amplitude: `{fmt(observed['pooled_amplitude'])}`",
        f"- amplitude ratio to locked A_ref: `{fmt(observed['amplitude_ratio_to_reference'])}`",
        f"- direction separation from locked axis: `{fmt(observed['direction_separation_deg'])} deg`",
        "",
        "## Conjunctive Rule",
        "",
        markdown_table(condition_rows, ["condition", "pass"]),
        "",
        "No in-sample BIC result is used in this rule. The reported BIC difference is diagnostic only and cannot substitute for positive held-out gain.",
        "",
        "## Power Calibration",
        "",
        markdown_table(power_rows, power_columns),
        "",
        "A non-zero replicate counts as recovered only when its held-out statistic exceeds the frozen threshold, its recovered direction is within 30 degrees, its amplitude ratio is in [0.5, 1.5], and its orientation is positive.",
        "",
        "## Held-Out Sectors",
        "",
        markdown_table(fold_rows, fold_columns),
        "",
        "## Template Inventory",
        "",
        f"The joint fit retained `{config['inputs']['template_count']}` nuisance templates and `{config['inputs']['joint_finite_rows']}` catalogue rows. The machine-readable inventory records deterministic omissions and standardisation constants.",
        "",
        "## Scope and Limitations",
        "",
        "This is a targeted, conditional power audit on the already selected Quaia candidate. The scalar shape comes from existing late-time best fits, and the injection null conditions on the observed catalogue, mask and nuisance fields. A pass would not constitute an independent detection. A failure does not exclude SU2 generally. An underpowered result does not permit the threshold to be moved.",
        "",
        "The SU2 scalar curve is evaluated at catalogue redshift and enters only through a centred angular interaction. This avoids adding the scalar curve alone as a circular predictor of redshift, but the result remains a phenomenological association test rather than a first-principles quasar-selection likelihood.",
    ]
    readout = [
        "# SU2 / Quaia Injection-Calibrated Angular-Likelihood Readout",
        "",
        "Date: 18 July 2026",
        "",
        f"Outcome: **{decision['outcome']}**.",
        "",
        f"Validation false-positive rate was `{fmt(calibration['validation_false_positive_rate'])}` and joint 1x recovery was `{fmt(calibration['joint_recovery_rate_1x'])}`. Observed pooled held-out gain was `{fmt(observed['delta_cvlpd'])}` against a frozen threshold of `{fmt(calibration['frozen_threshold'])}`, with `{observed['positive_fold_count']}` positive sectors, amplitude ratio `{fmt(observed['amplitude_ratio_to_reference'])}` and direction separation `{fmt(observed['direction_separation_deg'])} deg`.",
        "",
        decision["interpretation"],
        "",
        "In-sample BIC is non-gating and cannot replace held-out gain.",
    ]
    report_path = out_dir / f"{PREFIX}_report.md"
    readout_path = out_dir / f"{PREFIX}_readout.md"
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")
    readout_path.write_text("\n".join(readout) + "\n", encoding="utf-8")
    return report_path, readout_path


def export_outputs(out_dir: Path, local_paths: list[Path], script_path: Path) -> Path:
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    exported: list[Path] = []
    for source in local_paths:
        target = RESULT_DIR / f"{source.stem}_{DATE_SUFFIX}{source.suffix}"
        shutil.copy2(source, target)
        exported.append(target)
    manifest_rows = []
    for path in [script_path, PREREG_MD, PREREG_JSON, *exported]:
        manifest_rows.append(
            {
                "path": str(path.relative_to(ROOT)),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    manifest_path = RESULT_DIR / f"{PREFIX}_manifest_{DATE_SUFFIX}.csv"
    write_csv(manifest_path, manifest_rows)
    return manifest_path


def run(args: argparse.Namespace) -> dict[str, Any]:
    if not PREREG_MD.exists() or not PREREG_JSON.exists():
        raise FileNotFoundError("The dated preregistration pair must exist before outcome analysis.")
    if args.threshold_replicates >= args.n_replicates:
        raise ValueError("threshold-replicates must leave an untouched validation partition.")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    response, templates, su2_basis, l_deg, b_deg, template_names, inventory, input_meta, curve_rows = load_inputs(args)
    locked_unit = unit_from_lb(args.locked_l_deg, args.locked_b_deg)
    signal = args.reference_amplitude * (su2_basis @ locked_unit)
    x0 = np.column_stack([np.ones(len(response)), templates])
    x1 = np.column_stack([x0, su2_basis])
    beta0 = np.linalg.pinv(x0.T @ x0, rcond=1e-12) @ (x0.T @ response)
    base = x0 @ beta0
    residual = response - base

    cell, cell_fold, object_fold, spatial_meta = make_cells_and_folds(
        l_deg, b_deg, args.cell_degrees, args.n_folds
    )

    model0 = build_model_stats(x0, base, signal, residual, cell, object_fold, args.n_folds)
    model1 = build_model_stats(x1, base, signal, residual, cell, object_fold, args.n_folds)
    models = [model0, model1]
    response_stats = build_response_stats(
        base, signal, residual, cell, object_fold, args.n_folds, cell_fold
    )

    observed_signs = np.ones((spatial_meta["occupied_cells"], 1), dtype=float)
    observed_terms = precompute_sign_terms(models, response_stats, observed_signs, args.n_folds)
    observed_eval = evaluate_pair(models, response_stats, observed_terms, 0.0, args.n_folds)
    observed_vector = observed_eval["pooled_vector"][:, 0]
    observed_vector_diag = vector_diagnostics(
        observed_vector[:, None], locked_unit, args.reference_amplitude
    )
    observed_amp = float(observed_vector_diag["amplitude"][0])
    observed_l, observed_b = unit_to_lb(observed_vector)

    rng = np.random.default_rng(args.seed)
    signs = rng.choice(
        np.array([-1.0, 1.0]),
        size=(spatial_meta["occupied_cells"], args.n_replicates),
    )
    calibration_terms = precompute_sign_terms(models, response_stats, signs, args.n_folds)
    amplitude_scales = parse_csv_floats(args.amplitude_scales)
    if 0.0 not in amplitude_scales or 1.0 not in amplitude_scales:
        raise ValueError("The locked amplitude scales must include both 0 and 1.")

    evaluations: dict[float, dict[str, np.ndarray]] = {}
    for amplitude_scale in amplitude_scales:
        evaluations[amplitude_scale] = evaluate_pair(
            models,
            response_stats,
            calibration_terms,
            amplitude_scale,
            args.n_folds,
        )
        print(f"evaluated amplitude scale {amplitude_scale:g}", flush=True)

    calibration_slice = slice(0, args.threshold_replicates)
    validation_slice = slice(args.threshold_replicates, args.n_replicates)
    null_delta = evaluations[0.0]["delta_lpd"]
    threshold = float(
        np.quantile(null_delta[calibration_slice], 1.0 - args.alpha, method="higher")
    )
    validation_null = null_delta[validation_slice]
    validation_fpr = float(np.mean((validation_null > 0.0) & (validation_null > threshold)))

    calibration_rows: list[dict[str, Any]] = []
    power_rows: list[dict[str, Any]] = []
    recovery_by_scale: dict[float, np.ndarray] = {}
    for amplitude_scale in amplitude_scales:
        evaluation = evaluations[amplitude_scale]
        target_amp = amplitude_scale * args.reference_amplitude
        diag = vector_diagnostics(evaluation["pooled_vector"], locked_unit, target_amp)
        predictive = (evaluation["delta_lpd"] > 0.0) & (evaluation["delta_lpd"] > threshold)
        if amplitude_scale > 0.0:
            recovered = (
                predictive
                & diag["positive_orientation"]
                & (diag["amplitude_ratio"] >= args.amplitude_ratio_min)
                & (diag["amplitude_ratio"] <= args.amplitude_ratio_max)
                & (diag["direction_separation_deg"] <= args.maximum_direction_separation_deg)
            )
        else:
            recovered = predictive.copy()
        recovery_by_scale[amplitude_scale] = recovered

        for replicate in range(args.n_replicates):
            partition = "threshold_calibration" if replicate < args.threshold_replicates else "untouched_validation"
            calibration_rows.append(
                {
                    "replicate": replicate,
                    "partition": partition,
                    "amplitude_scale": amplitude_scale,
                    "delta_cvlpd": float(evaluation["delta_lpd"][replicate]),
                    "pooled_amplitude": float(diag["amplitude"][replicate]),
                    "amplitude_ratio": float(diag["amplitude_ratio"][replicate]),
                    "direction_separation_deg": float(diag["direction_separation_deg"][replicate]),
                    "positive_orientation": bool(diag["positive_orientation"][replicate]),
                    "predictive_detection": bool(predictive[replicate]),
                    "joint_recovery": bool(recovered[replicate]),
                }
            )

        for partition, subset in [
            ("threshold_calibration", calibration_slice),
            ("untouched_validation", validation_slice),
        ]:
            subset_delta = evaluation["delta_lpd"][subset]
            subset_predictive = predictive[subset]
            subset_recovered = recovered[subset]
            subset_ratio = diag["amplitude_ratio"][subset]
            subset_sep = diag["direction_separation_deg"][subset]
            power_rows.append(
                {
                    "amplitude_scale": amplitude_scale,
                    "partition": partition,
                    "n_replicates": int(len(subset_delta)),
                    "predictive_detection_rate": float(np.mean(subset_predictive)),
                    "joint_recovery_rate": float(np.mean(subset_recovered)) if amplitude_scale > 0.0 else float("nan"),
                    "delta_cvlpd_median": float(np.median(subset_delta)),
                    "delta_cvlpd_q05": float(np.quantile(subset_delta, 0.05)),
                    "delta_cvlpd_q95": float(np.quantile(subset_delta, 0.95)),
                    "amplitude_ratio_median": float(np.nanmedian(subset_ratio)) if amplitude_scale > 0.0 else float("nan"),
                    "direction_separation_median_deg": float(np.nanmedian(subset_sep)) if amplitude_scale > 0.0 else float("nan"),
                }
            )

    recovery_1x = float(np.mean(recovery_by_scale[1.0][validation_slice]))
    observed_delta = float(observed_eval["delta_lpd"][0])
    observed_fold_delta = observed_eval["fold_delta_lpd"][:, 0]
    positive_fold_count = int(np.sum(observed_fold_delta > 0.0))
    observed_sep = float(observed_vector_diag["direction_separation_deg"][0])
    observed_ratio = float(observed_vector_diag["amplitude_ratio"][0])
    observed_orientation = bool(observed_vector_diag["positive_orientation"][0])
    empirical_p = float((1 + np.sum(validation_null >= observed_delta)) / (1 + len(validation_null)))

    conditions = {
        "validation_null_false_positive_rate_le_0p01": validation_fpr <= args.alpha,
        "joint_1x_recovery_ge_0p80": recovery_1x >= args.minimum_power,
        "observed_pooled_held_out_gain_positive_and_above_threshold": observed_delta > 0.0 and observed_delta > threshold,
        "observed_direction_within_30_deg": observed_sep <= args.maximum_direction_separation_deg,
        "observed_amplitude_ratio_in_0p5_to_1p5": args.amplitude_ratio_min <= observed_ratio <= args.amplitude_ratio_max,
        "observed_orientation_positive": observed_orientation,
        "at_least_8_of_12_positive_held_out_sectors": positive_fold_count >= args.minimum_positive_folds,
    }
    false_positive_control_ok = conditions["validation_null_false_positive_rate_le_0p01"]
    reference_power_ok = conditions["joint_1x_recovery_ge_0p80"]
    if not false_positive_control_ok and not reference_power_ok:
        outcome = "calibration failure and underpowered"
        interpretation = "The untouched null false-positive rate exceeded 1%, and the locked 1x effect was not jointly recovered in at least 80% of validation injections. Both are calibration outcomes: the threshold remains unchanged and the observed candidate is not promotable."
    elif not false_positive_control_ok:
        outcome = "calibration failure"
        interpretation = "The frozen rule did not control the false-positive rate on the untouched null partition; the observed candidate is not promotable."
    elif not reference_power_ok:
        outcome = "underpowered"
        interpretation = "The locked 1x effect was not jointly recovered in at least 80% of validation injections. The threshold remains unchanged, and the observed result is not interpreted as a physical rejection."
    elif all(conditions.values()):
        outcome = "pass for independent follow-up only"
        interpretation = "Every calibration and observed held-out condition passed. Independent-catalogue validation is still required before promotion."
    else:
        outcome = "gate failure"
        interpretation = "False-positive control and 1x power were adequate, but the observed candidate failed one or more locked held-out direction, amplitude or stability conditions."

    in_sample = in_sample_diagnostics(x0, x1, response)
    in_sample_vector = np.asarray(in_sample.pop("candidate_vector"), dtype=float)
    in_sample_l, in_sample_b = unit_to_lb(in_sample_vector)
    in_sample.update(
        {
            "candidate_amplitude": float(np.linalg.norm(in_sample_vector)),
            "candidate_l_deg": in_sample_l,
            "candidate_b_deg": in_sample_b,
        }
    )

    fold_rows: list[dict[str, Any]] = []
    for fold in range(args.n_folds):
        vector = observed_eval["fold_candidate_beta"][fold, :, 0]
        fold_l, fold_b = unit_to_lb(vector)
        fold_rows.append(
            {
                "fold": fold,
                "longitude_lo_deg": fold * 360.0 / args.n_folds,
                "longitude_hi_deg": (fold + 1) * 360.0 / args.n_folds,
                "N_test": int(response_stats["n_test"][fold]),
                "delta_lpd": float(observed_eval["fold_delta_lpd"][fold, 0]),
                "delta_sse": float(observed_eval["fold_delta_sse"][fold, 0]),
                "recovered_amp": float(np.linalg.norm(vector)),
                "recovered_l_deg": fold_l,
                "recovered_b_deg": fold_b,
            }
        )

    decision = {
        "analysis": PREFIX,
        "date": DATE_SUFFIX,
        "outcome": outcome,
        "interpretation": interpretation,
        "conjunctive_rule": True,
        "threshold_frozen_after_calibration_partition": True,
        "threshold_moved_for_power": False,
        "in_sample_bic_can_substitute": False,
        "calibration": {
            "threshold_replicates": args.threshold_replicates,
            "validation_replicates": args.n_replicates - args.threshold_replicates,
            "frozen_threshold": threshold,
            "validation_false_positive_rate": validation_fpr,
            "joint_recovery_rate_1x": recovery_1x,
        },
        "observed": {
            "delta_cvlpd": observed_delta,
            "empirical_validation_null_p": empirical_p,
            "positive_fold_count": positive_fold_count,
            "pooled_amplitude": observed_amp,
            "amplitude_ratio_to_reference": observed_ratio,
            "pooled_l_deg": observed_l,
            "pooled_b_deg": observed_b,
            "direction_separation_deg": observed_sep,
            "positive_orientation": observed_orientation,
        },
        "conditions": conditions,
        "in_sample_non_gating_diagnostics": in_sample,
    }

    config = {
        "analysis": PREFIX,
        "date": DATE_SUFFIX,
        "preregistration_markdown": str(PREREG_MD),
        "preregistration_markdown_sha256": sha256_file(PREREG_MD),
        "preregistration_json": str(PREREG_JSON),
        "preregistration_json_sha256": sha256_file(PREREG_JSON),
        "sample": {"z_min": args.z_min, "z_max": args.z_max, "b_cut": args.b_cut},
        "locked_axis": {
            "l_deg": args.locked_l_deg,
            "b_deg": args.locked_b_deg,
            "reference_amplitude": args.reference_amplitude,
        },
        "held_out": spatial_meta,
        "calibration": {
            "seed": args.seed,
            "n_replicates": args.n_replicates,
            "threshold_replicates": args.threshold_replicates,
            "amplitude_scales": amplitude_scales,
            "alpha": args.alpha,
            "minimum_power": args.minimum_power,
            "direction_separation_max_deg": args.maximum_direction_separation_deg,
            "amplitude_ratio_range": [args.amplitude_ratio_min, args.amplitude_ratio_max],
        },
        "templates": template_names,
        "inputs": input_meta,
    }

    paths = {
        "calibration_draws": args.out_dir / f"{PREFIX}_calibration_draws.csv",
        "power_summary": args.out_dir / f"{PREFIX}_power_summary.csv",
        "heldout_folds": args.out_dir / f"{PREFIX}_heldout_folds.csv",
        "template_inventory": args.out_dir / f"{PREFIX}_template_inventory.csv",
        "scalar_curve": args.out_dir / f"{PREFIX}_scalar_curve.csv",
        "decision": args.out_dir / f"{PREFIX}_decision.json",
        "config": args.out_dir / f"{PREFIX}_config.json",
    }
    write_csv(paths["calibration_draws"], calibration_rows)
    write_csv(paths["power_summary"], power_rows)
    write_csv(paths["heldout_folds"], fold_rows)
    write_csv(paths["template_inventory"], inventory)
    write_csv(paths["scalar_curve"], curve_rows)
    write_json(paths["decision"], decision)
    write_json(paths["config"], config)
    report_path, readout_path = build_reports(
        args.out_dir, decision, power_rows, fold_rows, inventory, config
    )

    compact_paths = [
        paths["power_summary"],
        paths["heldout_folds"],
        paths["template_inventory"],
        paths["scalar_curve"],
        paths["decision"],
        paths["config"],
        report_path,
        readout_path,
    ]
    manifest_path = None
    if not args.no_export:
        manifest_path = export_outputs(args.out_dir, compact_paths, Path(__file__))

    print(f"Outcome: {outcome}", flush=True)
    print(f"Saved report: {report_path}", flush=True)
    print(f"Saved decision: {paths['decision']}", flush=True)
    if manifest_path:
        print(f"Saved export manifest: {manifest_path}", flush=True)
    return {"decision": decision, "config": config, "paths": paths}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Injection- and power-calibrated held-out SU2/Quaia angular likelihood.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--quaia", type=Path, default=EXT.DEFAULT_QUAIA)
    parser.add_argument("--randoms", type=Path, default=EXT.DEFAULT_RANDOMS)
    parser.add_argument("--selection", type=Path, default=EXT.DEFAULT_SELECTION)
    parser.add_argument("--sfd-dir", type=Path, default=EXT.DEFAULT_SFD_DIR)
    parser.add_argument("--wise-inventory", type=Path, default=WISE.DEFAULT_WISE)
    parser.add_argument("--gaia-density", type=Path, default=WISE.DEFAULT_GAIA_DENSITY)
    parser.add_argument("--gaia-density-hips-base-url", default=WISE.DEFAULT_GAIA_HIPS_BASE)
    parser.add_argument("--gaia-density-hips-cache-dir", type=Path, default=WISE.DEFAULT_GAIA_HIPS_CACHE)
    parser.add_argument("--gaia-density-hips-order", type=int, default=2)
    parser.add_argument("--gaiascanlaw-data-dir", type=Path, default=None)
    parser.add_argument("--skip-gaia-density", action="store_true")
    parser.add_argument("--lcdm-bestfit", type=Path, default=LCDM_BESTFIT)
    parser.add_argument("--su2r-bestfit", type=Path, default=SU2R_BESTFIT)
    parser.add_argument("--z-min", type=float, default=1.0)
    parser.add_argument("--z-max", type=float, default=1.5)
    parser.add_argument("--b-cut", type=float, default=15.0)
    parser.add_argument("--locked-l-deg", type=float, default=138.99062177880762)
    parser.add_argument("--locked-b-deg", type=float, default=-36.60410263509791)
    parser.add_argument("--reference-amplitude", type=float, default=0.002548691334890594)
    parser.add_argument("--cell-degrees", type=float, default=12.0)
    parser.add_argument("--n-folds", type=int, default=12)
    parser.add_argument("--n-replicates", type=int, default=2000)
    parser.add_argument("--threshold-replicates", type=int, default=1000)
    parser.add_argument("--amplitude-scales", default="0,0.5,1,1.5,2")
    parser.add_argument("--seed", type=int, default=180718)
    parser.add_argument("--alpha", type=float, default=0.01)
    parser.add_argument("--minimum-power", type=float, default=0.80)
    parser.add_argument("--minimum-positive-folds", type=int, default=8)
    parser.add_argument("--maximum-direction-separation-deg", type=float, default=30.0)
    parser.add_argument("--amplitude-ratio-min", type=float, default=0.5)
    parser.add_argument("--amplitude-ratio-max", type=float, default=1.5)
    parser.add_argument("--scalar-grid-points", type=int, default=501)
    parser.add_argument("--integration-nint", type=int, default=256)
    parser.add_argument("--no-export", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
