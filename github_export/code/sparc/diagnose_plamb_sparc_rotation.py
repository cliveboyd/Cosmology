#!/usr/bin/env python3
r"""
Run a first PLAMB-style audit against SPARC galaxy rotation curves.

This is a phenomenological diagnostic, not a full galaxy-formation or
action-level theory fit. It compares observed SPARC rotation velocities with:

  1. Newtonian baryons with fixed stellar mass-to-light ratios.
  2. Newtonian baryons with global fitted disk/bulge mass-to-light ratios.
  3. A standard radial-acceleration-relation style acceleration boost.
  4. A PLAMB-local-inertia power boost tied to low baryonic acceleration.

Outputs:
    plamb_runs/diagnostics/sparc_rotation_audit/sparc_rotation_fit_summary.csv
    plamb_runs/diagnostics/sparc_rotation_audit/sparc_rotation_galaxy_residuals.csv
    plamb_runs/diagnostics/sparc_rotation_audit/sparc_rotation_acceleration_points.csv
    plamb_runs/diagnostics/sparc_rotation_audit/sparc_rotation_report.md
    plamb_runs/diagnostics/sparc_rotation_audit/sparc_rotation_acceleration_plane.png

Run:
    python diagnose_plamb_sparc_rotation.py
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

import numpy as np

try:
    from scipy.optimize import differential_evolution, minimize
except Exception as exc:  # pragma: no cover - reported at runtime
    differential_evolution = None
    minimize = None
    SCIPY_IMPORT_ERROR = exc
else:
    SCIPY_IMPORT_ERROR = None

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as exc:  # pragma: no cover - plotting is optional
    plt = None
    MATPLOTLIB_IMPORT_ERROR = exc
else:
    MATPLOTLIB_IMPORT_ERROR = None


ROOT = Path(__file__).resolve().parents[3]
SPARC_DIR = ROOT / "external_datasets" / "current_cosmology_sources" / "SPARC"
DEFAULT_SAMPLE = SPARC_DIR / "sparc_galaxy_sample.csv"
DEFAULT_POINTS = SPARC_DIR / "sparc_rotation_points.csv"
BASE_OUTDIR = ROOT / "plamb_runs" / "diagnostics" / "sparc_rotation_audit"
OUTDIR = BASE_OUTDIR

KPC_M = 3.0856775814913673e19
DEFAULT_G0 = 1.2e-10


@dataclass
class Dataset:
    galaxy: np.ndarray
    rad_kpc: np.ndarray
    vobs: np.ndarray
    sigma_v: np.ndarray
    vgas: np.ndarray
    vdisk: np.ndarray
    vbul: np.ndarray
    quality: np.ndarray
    distance_method: np.ndarray
    inc_deg: np.ndarray
    d_mpc: np.ndarray


@dataclass
class ModelSpec:
    name: str
    param_names: list[str]
    bounds: list[tuple[float, float]]
    initial: list[float]
    predictor: Callable[[np.ndarray, Dataset], np.ndarray]
    note: str


def parse_float(value: Any, default: float = float("nan")) -> float:
    if value in (None, ""):
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def parse_int(value: Any, default: int = -1) -> int:
    if value in (None, ""):
        return default
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def load_dataset(sample_path: Path, points_path: Path, quality_max: int, distance_method: str, err_floor: float) -> Dataset:
    if not sample_path.exists() or not points_path.exists():
        raise FileNotFoundError(
            "SPARC prepared CSVs were not found. Run: python fetch_sparc_rotation.py"
        )
    sample_rows = read_csv(sample_path)
    sample_by_galaxy = {row["Galaxy"]: row for row in sample_rows}

    rows: list[dict[str, Any]] = []
    for row in read_csv(points_path):
        galaxy = row.get("Galaxy", "")
        sample = sample_by_galaxy.get(galaxy, {})
        q = parse_int(sample.get("Q"), 99)
        f_d = parse_int(sample.get("f_D"), -1)
        if q > quality_max:
            continue
        if distance_method == "non_hubble" and f_d == 1:
            continue
        if distance_method == "hubble_flow" and f_d != 1:
            continue
        rad = parse_float(row.get("Rad_kpc"))
        vobs = parse_float(row.get("Vobs_kms"))
        err = parse_float(row.get("errV_kms"))
        if not (math.isfinite(rad) and math.isfinite(vobs) and math.isfinite(err)):
            continue
        if rad <= 0.0 or vobs <= 0.0 or err < 0.0:
            continue
        rows.append(
            {
                "galaxy": galaxy,
                "rad": rad,
                "vobs": vobs,
                "sigma_v": math.sqrt(err * err + err_floor * err_floor),
                "vgas": parse_float(row.get("Vgas_kms"), 0.0),
                "vdisk": parse_float(row.get("Vdisk_kms"), 0.0),
                "vbul": parse_float(row.get("Vbul_kms"), 0.0),
                "quality": q,
                "distance_method": f_d,
                "inc_deg": parse_float(sample.get("Inc_deg")),
                "d_mpc": parse_float(sample.get("D_Mpc")),
            }
        )

    if not rows:
        raise RuntimeError("No SPARC rows survived the requested filters.")

    return Dataset(
        galaxy=np.array([row["galaxy"] for row in rows], dtype=object),
        rad_kpc=np.array([row["rad"] for row in rows], dtype=float),
        vobs=np.array([row["vobs"] for row in rows], dtype=float),
        sigma_v=np.array([row["sigma_v"] for row in rows], dtype=float),
        vgas=np.array([row["vgas"] for row in rows], dtype=float),
        vdisk=np.array([row["vdisk"] for row in rows], dtype=float),
        vbul=np.array([row["vbul"] for row in rows], dtype=float),
        quality=np.array([row["quality"] for row in rows], dtype=int),
        distance_method=np.array([row["distance_method"] for row in rows], dtype=int),
        inc_deg=np.array([row["inc_deg"] for row in rows], dtype=float),
        d_mpc=np.array([row["d_mpc"] for row in rows], dtype=float),
    )


def baryon_v2(data: Dataset, ydisk: float, ybul: float) -> np.ndarray:
    # SPARC gas terms can be negative in some regions; preserve the sign of the
    # gas contribution while keeping stellar disk/bulge contributions positive.
    return data.vgas * np.abs(data.vgas) + ydisk * data.vdisk**2 + ybul * data.vbul**2


def acceleration_from_v2(v2_kms2: np.ndarray, radius_kpc: np.ndarray) -> np.ndarray:
    return np.maximum(v2_kms2, 1e-18) * 1.0e6 / (radius_kpc * KPC_M)


def velocity_from_acceleration(g_si: np.ndarray, radius_kpc: np.ndarray) -> np.ndarray:
    return np.sqrt(np.maximum(g_si, 0.0) * radius_kpc * KPC_M) / 1000.0


def fixed_baryon(theta: np.ndarray, data: Dataset) -> np.ndarray:
    del theta
    return np.sqrt(np.maximum(baryon_v2(data, 0.5, 0.7), 0.0))


def free_baryon(theta: np.ndarray, data: Dataset) -> np.ndarray:
    ydisk, ybul = theta
    return np.sqrt(np.maximum(baryon_v2(data, ydisk, ybul), 0.0))


def rar_predictor(theta: np.ndarray, data: Dataset) -> np.ndarray:
    ydisk, ybul, log10_gdagger = theta
    v2_bar = baryon_v2(data, ydisk, ybul)
    gbar = acceleration_from_v2(v2_bar, data.rad_kpc)
    gdagger = 10.0**log10_gdagger
    x = np.sqrt(np.maximum(gbar / gdagger, 1e-30))
    denominator = np.maximum(1.0 - np.exp(-x), 1e-12)
    gpred = gbar / denominator
    return velocity_from_acceleration(gpred, data.rad_kpc)


def plamb_inertia_power(theta: np.ndarray, data: Dataset) -> np.ndarray:
    ydisk, ybul, amp, index = theta
    v2_bar = baryon_v2(data, ydisk, ybul)
    gbar = acceleration_from_v2(v2_bar, data.rad_kpc)
    ratio = np.clip(DEFAULT_G0 / np.maximum(gbar, 1e-30), 1e-8, 1e8)
    boost = 1.0 + amp * ratio**index
    return velocity_from_acceleration(gbar * boost, data.rad_kpc)


def plamb_inertia_power_free_g0(theta: np.ndarray, data: Dataset) -> np.ndarray:
    ydisk, ybul, amp, index, log10_g0 = theta
    v2_bar = baryon_v2(data, ydisk, ybul)
    gbar = acceleration_from_v2(v2_bar, data.rad_kpc)
    g0 = 10.0**log10_g0
    ratio = np.clip(g0 / np.maximum(gbar, 1e-30), 1e-8, 1e8)
    boost = 1.0 + amp * ratio**index
    return velocity_from_acceleration(gbar * boost, data.rad_kpc)


def model_specs(include_free_g0: bool) -> list[ModelSpec]:
    specs = [
        ModelSpec(
            name="BARYON_FIXED_ML",
            param_names=[],
            bounds=[],
            initial=[],
            predictor=fixed_baryon,
            note="Newtonian baryons; fixed disk M/L=0.5, bulge M/L=0.7.",
        ),
        ModelSpec(
            name="BARYON_FREE_ML",
            param_names=["ydisk", "ybul"],
            bounds=[(0.05, 1.20), (0.05, 1.50)],
            initial=[0.5, 0.7],
            predictor=free_baryon,
            note="Newtonian baryons; global disk and bulge mass-to-light ratios fitted.",
        ),
        ModelSpec(
            name="RAR_FREE_GDAGGER",
            param_names=["ydisk", "ybul", "log10_gdagger"],
            bounds=[(0.05, 1.20), (0.05, 1.50), (-11.5, -9.3)],
            initial=[0.5, 0.7, math.log10(DEFAULT_G0)],
            predictor=rar_predictor,
            note="Radial-acceleration relation style interpolation with fitted g_dagger.",
        ),
        ModelSpec(
            name="PLAMB_INERTIA_POWER",
            param_names=["ydisk", "ybul", "inertia_amp", "inertia_index"],
            bounds=[(0.05, 1.20), (0.05, 1.50), (0.0, 3.0), (0.0, 1.2)],
            initial=[0.5, 0.7, 0.5, 0.5],
            predictor=plamb_inertia_power,
            note="PLAMB-style local inertia boost: g_pred = g_bar * (1 + A*(g0/g_bar)^n), fixed g0=1.2e-10 m/s^2.",
        ),
    ]
    if include_free_g0:
        specs.append(
            ModelSpec(
                name="PLAMB_INERTIA_POWER_FREE_G0",
                param_names=["ydisk", "ybul", "inertia_amp", "inertia_index", "log10_g0"],
                bounds=[(0.05, 1.20), (0.05, 1.50), (0.0, 3.0), (0.0, 1.2), (-11.7, -9.2)],
                initial=[0.5, 0.7, 0.5, 0.5, math.log10(DEFAULT_G0)],
                predictor=plamb_inertia_power_free_g0,
                note="PLAMB-style local inertia boost with fitted acceleration pivot g0.",
            )
        )
    return specs


def chi2_for(theta: np.ndarray, spec: ModelSpec, data: Dataset) -> float:
    try:
        pred = spec.predictor(theta, data)
    except FloatingPointError:
        return 1.0e100
    if pred.shape != data.vobs.shape or not np.all(np.isfinite(pred)):
        return 1.0e100
    residual = (data.vobs - pred) / data.sigma_v
    value = float(np.sum(residual**2))
    if not math.isfinite(value):
        return 1.0e100
    return value


def fit_model(spec: ModelSpec, data: Dataset, seed: int, maxiter: int) -> tuple[np.ndarray, float]:
    if not spec.bounds:
        theta = np.array([], dtype=float)
        return theta, chi2_for(theta, spec, data)
    if differential_evolution is None or minimize is None:
        raise RuntimeError(f"SciPy optimizer import failed: {SCIPY_IMPORT_ERROR}")

    objective = lambda theta: chi2_for(np.asarray(theta, dtype=float), spec, data)
    de = differential_evolution(
        objective,
        bounds=spec.bounds,
        seed=seed,
        tol=1e-6,
        polish=False,
        maxiter=maxiter,
        popsize=12,
        updating="immediate",
        workers=1,
    )
    local = minimize(
        objective,
        de.x,
        method="L-BFGS-B",
        bounds=spec.bounds,
        options={"maxiter": 2000, "ftol": 1e-9},
    )
    if local.fun <= de.fun and local.success:
        return np.asarray(local.x, dtype=float), float(local.fun)
    return np.asarray(de.x, dtype=float), float(de.fun)


def summarize_fit(spec: ModelSpec, theta: np.ndarray, chi2: float, n: int) -> dict[str, Any]:
    k = len(theta)
    aic = chi2 + 2.0 * k
    bic = chi2 + k * math.log(max(n, 1))
    row: dict[str, Any] = {
        "model": spec.name,
        "k": k,
        "N": n,
        "chi2": chi2,
        "dof": n - k,
        "chi2_dof": chi2 / max(n - k, 1),
        "AIC": aic,
        "BIC": bic,
        "note": spec.note,
    }
    for name, value in zip(spec.param_names, theta):
        row[name] = value
    if spec.name == "PLAMB_INERTIA_POWER":
        row["log10_g0"] = math.log10(DEFAULT_G0)
    return row


def per_galaxy_residuals(data: Dataset, fits: dict[str, tuple[ModelSpec, np.ndarray]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for galaxy in sorted(set(data.galaxy.tolist())):
        mask = data.galaxy == galaxy
        if not np.any(mask):
            continue
        sub = Dataset(
            galaxy=data.galaxy[mask],
            rad_kpc=data.rad_kpc[mask],
            vobs=data.vobs[mask],
            sigma_v=data.sigma_v[mask],
            vgas=data.vgas[mask],
            vdisk=data.vdisk[mask],
            vbul=data.vbul[mask],
            quality=data.quality[mask],
            distance_method=data.distance_method[mask],
            inc_deg=data.inc_deg[mask],
            d_mpc=data.d_mpc[mask],
        )
        row: dict[str, Any] = {
            "Galaxy": galaxy,
            "N_points": int(np.sum(mask)),
            "Q": int(sub.quality[0]) if len(sub.quality) else "",
            "f_D": int(sub.distance_method[0]) if len(sub.distance_method) else "",
            "Inc_deg": float(sub.inc_deg[0]) if len(sub.inc_deg) else "",
            "D_Mpc": float(sub.d_mpc[0]) if len(sub.d_mpc) else "",
            "Vobs_outer_kms": float(sub.vobs[np.argmax(sub.rad_kpc)]),
            "Rad_outer_kpc": float(np.max(sub.rad_kpc)),
        }
        for name, (spec, theta) in fits.items():
            pred = spec.predictor(theta, sub)
            residual = (sub.vobs - pred) / sub.sigma_v
            row[f"chi2_{name}"] = float(np.sum(residual**2))
            row[f"rms_residual_kms_{name}"] = float(np.sqrt(np.mean((sub.vobs - pred) ** 2)))
        rows.append(row)
    return rows


def acceleration_rows(data: Dataset, best_name: str, spec: ModelSpec, theta: np.ndarray) -> list[dict[str, Any]]:
    ydisk = float(theta[0]) if len(theta) >= 1 else 0.5
    ybul = float(theta[1]) if len(theta) >= 2 else 0.7
    v2_bar = baryon_v2(data, ydisk, ybul)
    gbar = acceleration_from_v2(v2_bar, data.rad_kpc)
    gobs = acceleration_from_v2(data.vobs**2, data.rad_kpc)
    pred = spec.predictor(theta, data)
    gpred = acceleration_from_v2(pred**2, data.rad_kpc)
    rows = []
    for i in range(len(data.vobs)):
        rows.append(
            {
                "Galaxy": str(data.galaxy[i]),
                "Rad_kpc": data.rad_kpc[i],
                "Vobs_kms": data.vobs[i],
                "Vpred_kms": pred[i],
                "residual_kms": data.vobs[i] - pred[i],
                "gbar_si": gbar[i],
                "gobs_si": gobs[i],
                "gpred_si": gpred[i],
                "best_model": best_name,
                "Q": int(data.quality[i]),
                "f_D": int(data.distance_method[i]),
            }
        )
    return rows


def make_plot(data: Dataset, best_name: str, spec: ModelSpec, theta: np.ndarray, path: Path) -> None:
    if plt is None:
        return
    ydisk = float(theta[0]) if len(theta) >= 1 else 0.5
    ybul = float(theta[1]) if len(theta) >= 2 else 0.7
    v2_bar = baryon_v2(data, ydisk, ybul)
    gbar = acceleration_from_v2(v2_bar, data.rad_kpc)
    gobs = acceleration_from_v2(data.vobs**2, data.rad_kpc)

    xgrid = np.logspace(np.log10(max(np.nanmin(gbar[gbar > 0]), 1e-14)), np.log10(np.nanmax(gbar) * 1.2), 240)
    fake = Dataset(
        galaxy=np.array(["curve"] * len(xgrid), dtype=object),
        rad_kpc=np.ones(len(xgrid), dtype=float),
        vobs=np.ones(len(xgrid), dtype=float),
        sigma_v=np.ones(len(xgrid), dtype=float),
        vgas=np.sqrt(xgrid * KPC_M) / 1000.0,
        vdisk=np.zeros(len(xgrid), dtype=float),
        vbul=np.zeros(len(xgrid), dtype=float),
        quality=np.ones(len(xgrid), dtype=int),
        distance_method=np.ones(len(xgrid), dtype=int),
        inc_deg=np.zeros(len(xgrid), dtype=float),
        d_mpc=np.zeros(len(xgrid), dtype=float),
    )
    if best_name == "BARYON_FIXED_ML":
        y_curve = xgrid
    else:
        pred = spec.predictor(theta, fake)
        y_curve = acceleration_from_v2(pred**2, fake.rad_kpc)

    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    ax.scatter(gbar, gobs, s=6, alpha=0.28, linewidths=0, color="#34495e", label="SPARC points")
    ax.plot(xgrid, xgrid, color="#7f8c8d", linewidth=1.3, label="g_obs = g_bar")
    ax.plot(xgrid, y_curve, color="#c0392b", linewidth=2.0, label=f"Best: {best_name}")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Baryonic acceleration g_bar [m/s^2]")
    ax.set_ylabel("Observed/predicted acceleration [m/s^2]")
    ax.set_title("SPARC Acceleration Plane")
    ax.grid(True, which="both", alpha=0.18)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def write_report(
    data: Dataset,
    summary_rows: list[dict[str, Any]],
    args: argparse.Namespace,
    best_aic: dict[str, Any],
    best_bic: dict[str, Any],
) -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    quality_counts: dict[str, int] = {}
    distance_counts: dict[str, int] = {}
    for value in data.quality:
        quality_counts[str(int(value))] = quality_counts.get(str(int(value)), 0) + 1
    for value in data.distance_method:
        distance_counts[str(int(value))] = distance_counts.get(str(int(value)), 0) + 1

    lines = [
        "# SPARC Rotation Audit for PLAMB/SU2",
        "",
        "## Executive Summary",
        "",
        f"- Points fitted: {len(data.vobs)}",
        f"- Galaxies fitted: {len(set(data.galaxy.tolist()))}",
        f"- Quality cut: Q <= {args.quality_max}",
        f"- Distance-method filter: {args.distance_method}",
        f"- Velocity error floor: {args.err_floor_kms:.3g} km/s",
        f"- Best AIC model: {best_aic['model']} (AIC={float(best_aic['AIC']):.3f}, chi2/dof={float(best_aic['chi2_dof']):.3f})",
        f"- Best BIC model: {best_bic['model']} (BIC={float(best_bic['BIC']):.3f}, chi2/dof={float(best_bic['chi2_dof']):.3f})",
        "",
        "## Model Branches",
        "",
        "- `BARYON_FIXED_ML`: Newtonian baryons with fixed disk/bulge mass-to-light ratios.",
        "- `BARYON_FREE_ML`: Newtonian baryons with global fitted disk/bulge mass-to-light ratios.",
        "- `RAR_FREE_GDAGGER`: radial-acceleration-relation form with fitted acceleration scale.",
        "- `PLAMB_INERTIA_POWER`: local-inertia boost tied to low baryonic acceleration.",
        "",
        "The PLAMB audit branch uses:",
        "",
        "`g_pred = g_bar * (1 + A * (g0/g_bar)^n)`",
        "",
        f"with fixed `g0 = {DEFAULT_G0:.3e} m/s^2` unless the optional free-g0 branch is enabled.",
        "",
        "## Fit Table",
        "",
        "| Model | k | chi2 | chi2/dof | AIC | BIC | Key parameters |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in summary_rows:
        params = []
        for key, value in row.items():
            if key in {"model", "k", "N", "chi2", "dof", "chi2_dof", "AIC", "BIC", "note"}:
                continue
            try:
                params.append(f"{key}={float(value):.5g}")
            except (TypeError, ValueError):
                params.append(f"{key}={value}")
        lines.append(
            f"| {row['model']} | {row['k']} | {float(row['chi2']):.3f} | "
            f"{float(row['chi2_dof']):.3f} | {float(row['AIC']):.3f} | "
            f"{float(row['BIC']):.3f} | {', '.join(params)} |"
        )

    lines.extend(
        [
            "",
            "## Dataset Assumption Notes",
            "",
            "- SPARC is a galaxy-rotation dataset, not a cosmological distance ladder dataset.",
            "- The sample table includes a distance-method code. Method `f_D=1` uses Hubble-flow distances with H0=73 and Virgo-centric infall correction in the original SPARC metadata.",
            "- A strong result must therefore be repeated with `--distance-method non_hubble` and by quality cuts before it is used in a paper-level PLAMB claim.",
            "- Global stellar mass-to-light ratios are intentionally simple here. A full paper-grade fit would allow controlled galaxy-level nuisance structure or hierarchical priors.",
            "",
            "## Outputs",
            "",
            f"- Fit summary: `{OUTDIR / 'sparc_rotation_fit_summary.csv'}`",
            f"- Per-galaxy residuals: `{OUTDIR / 'sparc_rotation_galaxy_residuals.csv'}`",
            f"- Acceleration points: `{OUTDIR / 'sparc_rotation_acceleration_points.csv'}`",
            f"- Plot: `{OUTDIR / 'sparc_rotation_acceleration_plane.png'}`",
            "",
            "## Suggested Next Checks",
            "",
            "1. Repeat with `--distance-method non_hubble` to remove Hubble-flow distances.",
            "2. Repeat with `--quality-max 1` for only highest-quality rotation curves.",
            "3. Add a hierarchical per-galaxy mass-to-light nuisance model.",
            "4. Compare residual directions with galaxy inclination, surface brightness, and environment.",
            "5. Add Galaxy Zoo/MaNGA spin-axis data as a separate angular/parity audit for SU2.",
        ]
    )
    (OUTDIR / "sparc_rotation_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    meta = {
        "quality_counts_points": quality_counts,
        "distance_method_counts_points": distance_counts,
        "best_aic_model": best_aic["model"],
        "best_bic_model": best_bic["model"],
        "args": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
    }
    (OUTDIR / "sparc_rotation_audit_metadata.json").write_text(json.dumps(meta, indent=2, sort_keys=True), encoding="utf-8")


def main() -> None:
    global OUTDIR

    parser = argparse.ArgumentParser(description="Run first PLAMB/SPARC rotation-curve diagnostic.")
    parser.add_argument("--sample", type=Path, default=DEFAULT_SAMPLE)
    parser.add_argument("--points", type=Path, default=DEFAULT_POINTS)
    parser.add_argument("--quality-max", type=int, default=2, help="Maximum SPARC quality flag Q to include.")
    parser.add_argument(
        "--distance-method",
        choices=["all", "non_hubble", "hubble_flow"],
        default="all",
        help="Filter by SPARC distance method. non_hubble excludes f_D=1.",
    )
    parser.add_argument("--err-floor-kms", type=float, default=3.0)
    parser.add_argument("--seed", type=int, default=260713)
    parser.add_argument("--maxiter", type=int, default=220)
    parser.add_argument("--include-free-g0", action="store_true")
    parser.add_argument(
        "--output-label",
        default="",
        help="Optional subdirectory label under sparc_rotation_audit for variant runs.",
    )
    args = parser.parse_args()
    if args.output_label:
        safe_label = "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in args.output_label)
        OUTDIR = BASE_OUTDIR / safe_label

    data = load_dataset(args.sample, args.points, args.quality_max, args.distance_method, args.err_floor_kms)
    specs = model_specs(include_free_g0=args.include_free_g0)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    summary_rows: list[dict[str, Any]] = []
    fits: dict[str, tuple[ModelSpec, np.ndarray]] = {}
    for index, spec in enumerate(specs, start=1):
        print(f"[{index}/{len(specs)}] fitting {spec.name}")
        theta, chi2 = fit_model(spec, data, seed=args.seed + index, maxiter=args.maxiter)
        summary = summarize_fit(spec, theta, chi2, len(data.vobs))
        summary_rows.append(summary)
        fits[spec.name] = (spec, theta)
        print(
            f"    chi2={summary['chi2']:.3f} chi2/dof={summary['chi2_dof']:.3f} "
            f"AIC={summary['AIC']:.3f} BIC={summary['BIC']:.3f}"
        )

    summary_rows.sort(key=lambda row: float(row["AIC"]))
    fieldnames = [
        "model",
        "k",
        "N",
        "chi2",
        "dof",
        "chi2_dof",
        "AIC",
        "BIC",
        "ydisk",
        "ybul",
        "log10_gdagger",
        "inertia_amp",
        "inertia_index",
        "log10_g0",
        "note",
    ]
    write_csv(OUTDIR / "sparc_rotation_fit_summary.csv", summary_rows, fieldnames)

    galaxy_rows = per_galaxy_residuals(data, fits)
    galaxy_fields = sorted({key for row in galaxy_rows for key in row.keys()})
    write_csv(OUTDIR / "sparc_rotation_galaxy_residuals.csv", galaxy_rows, galaxy_fields)

    best_aic = min(summary_rows, key=lambda row: float(row["AIC"]))
    best_bic = min(summary_rows, key=lambda row: float(row["BIC"]))
    best_spec, best_theta = fits[str(best_aic["model"])]
    acc_rows = acceleration_rows(data, str(best_aic["model"]), best_spec, best_theta)
    write_csv(
        OUTDIR / "sparc_rotation_acceleration_points.csv",
        acc_rows,
        [
            "Galaxy",
            "Rad_kpc",
            "Vobs_kms",
            "Vpred_kms",
            "residual_kms",
            "gbar_si",
            "gobs_si",
            "gpred_si",
            "best_model",
            "Q",
            "f_D",
        ],
    )
    make_plot(data, str(best_aic["model"]), best_spec, best_theta, OUTDIR / "sparc_rotation_acceleration_plane.png")
    write_report(data, summary_rows, args, best_aic, best_bic)

    print(f"Saved summary: {OUTDIR / 'sparc_rotation_fit_summary.csv'}")
    print(f"Saved report: {OUTDIR / 'sparc_rotation_report.md'}")
    print(
        f"Best AIC: {best_aic['model']} chi2/dof={float(best_aic['chi2_dof']):.4g} "
        f"AIC={float(best_aic['AIC']):.4g}"
    )


if __name__ == "__main__":
    main()
