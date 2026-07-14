#!/usr/bin/env python3
r"""
Diagnostic SU2/LCDM comparison against GWTC-5.0 / O4b distance-redshift rows.

This is not a full gravitational-wave likelihood. It uses the published GWOSC
summary columns for redshift and luminosity distance where both are present,
with asymmetric distance errors compressed into a conservative log-distance
sigma. The result is a cross-probe diagnostic only.

Outputs:
    plamb_runs/diagnostics/su2_gwtc5_o4b/su2_gwtc5_o4b_fit_summary.csv
    plamb_runs/diagnostics/su2_gwtc5_o4b/su2_gwtc5_o4b_scored_bestfits.csv
    plamb_runs/diagnostics/su2_gwtc5_o4b/su2_gwtc5_o4b_event_residuals.csv
    plamb_runs/diagnostics/su2_gwtc5_o4b/su2_gwtc5_o4b_report.md
    plamb_runs/diagnostics/su2_gwtc5_o4b/su2_gwtc5_o4b_hubble.png

Run:
    python diagnose_su2_gwtc5_o4b.py
"""

from __future__ import annotations

import argparse
import ast
import csv
import math
import sys
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
except Exception as exc:  # pragma: no cover - optional plotting
    plt = None
    MATPLOTLIB_IMPORT_ERROR = exc
else:
    MATPLOTLIB_IMPORT_ERROR = None

from PLamb_Test_10V6c_plus import luminosity_distance


ROOT = Path(__file__).resolve().parent
DEFAULT_EVENT_CSV = (
    ROOT
    / "external_datasets"
    / "current_cosmology_sources"
    / "GWTC5_O4b"
    / "gwtc5_o4b_event_summary.csv"
)
OUTDIR = ROOT / "plamb_runs" / "diagnostics" / "su2_gwtc5_o4b"


@dataclass
class EventRow:
    event_id: str
    z: float
    z_lower: float
    z_upper: float
    dl_mpc: float
    dl_lower: float
    dl_upper: float
    sigma_ln_dl: float
    p_astro: float
    network_snr: float
    utc_approx: str


@dataclass
class ModelSpec:
    name: str
    param_names: list[str]
    bounds: list[tuple[float, float]]
    initial: list[float]
    note: str
    to_cosmo: Callable[[np.ndarray], dict[str, float]] 


def parse_float(value: Any, default: float = float("nan")) -> float:
    if value in (None, ""):
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def log_distance_sigma(dl: float, lower: float, upper: float, floor: float) -> float:
    lower_abs = abs(lower) if math.isfinite(lower) else 0.0
    upper_abs = abs(upper) if math.isfinite(upper) else 0.0
    pieces: list[float] = []
    if upper_abs > 0.0:
        pieces.append(math.log(max(dl + upper_abs, 1e-12) / dl))
    if lower_abs > 0.0 and dl - lower_abs > 0.0:
        pieces.append(math.log(dl / max(dl - lower_abs, 1e-12)))
    elif lower_abs > 0.0:
        pieces.append(lower_abs / max(dl, 1e-12))
    if not pieces:
        return floor
    return max(float(np.mean(pieces)), floor)


def load_events(path: Path, min_pastro: float, sigma_floor: float) -> list[EventRow]:
    rows: list[EventRow] = []
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            z = parse_float(row.get("redshift"))
            dl = parse_float(row.get("luminosity_distance"))
            p_astro = parse_float(row.get("p_astro"), 0.0)
            if not (math.isfinite(z) and math.isfinite(dl) and z > 0.0 and dl > 0.0):
                continue
            if p_astro < min_pastro:
                continue
            dl_lower = parse_float(row.get("luminosity_distance_lower"), 0.0)
            dl_upper = parse_float(row.get("luminosity_distance_upper"), 0.0)
            rows.append(
                EventRow(
                    event_id=str(row.get("event_id", "")),
                    z=z,
                    z_lower=parse_float(row.get("redshift_lower"), 0.0),
                    z_upper=parse_float(row.get("redshift_upper"), 0.0),
                    dl_mpc=dl,
                    dl_lower=dl_lower,
                    dl_upper=dl_upper,
                    sigma_ln_dl=log_distance_sigma(dl, dl_lower, dl_upper, sigma_floor),
                    p_astro=p_astro,
                    network_snr=parse_float(row.get("network_snr"), float("nan")),
                    utc_approx=str(row.get("utc_approx", "")),
                )
            )
    rows.sort(key=lambda event: event.z)
    return rows


def model_specs() -> list[ModelSpec]:
    def lcdm(theta: np.ndarray) -> dict[str, float]:
        H0, Om = theta
        Ol = 1.0 - Om
        return {"H0": H0, "Om": Om, "Ol": Ol, "Omega_chi0": 0.0, "w0_chi": -1.0, "wa_chi": 0.0}

    def su2(theta: np.ndarray) -> dict[str, float]:
        H0, Om, Omega_chi0, w0_chi = theta
        Ol = 1.0 - Om - Omega_chi0
        return {"H0": H0, "Om": Om, "Ol": Ol, "Omega_chi0": Omega_chi0, "w0_chi": w0_chi, "wa_chi": 0.0}

    def su2r(theta: np.ndarray) -> dict[str, float]:
        H0, Om, su2_fraction, w0_chi = theta
        dark = 1.0 - Om
        Omega_chi0 = su2_fraction * dark
        Ol = (1.0 - su2_fraction) * dark
        return {
            "H0": H0,
            "Om": Om,
            "Ol": Ol,
            "Omega_chi0": Omega_chi0,
            "su2_fraction": su2_fraction,
            "Omega_dark0": dark,
            "w0_chi": w0_chi,
            "wa_chi": 0.0,
        }

    return [
        ModelSpec(
            name="LCDM",
            param_names=["H0", "Om"],
            bounds=[(45.0, 95.0), (0.05, 0.60)],
            initial=[70.0, 0.30],
            note="Flat LCDM, H0 and Om free.",
            to_cosmo=lcdm,
        ),
        ModelSpec(
            name="SU2",
            param_names=["H0", "Om", "Omega_chi0", "w0_chi"],
            bounds=[(45.0, 95.0), (0.05, 0.60), (0.0, 0.90), (-3.0, -0.20)],
            initial=[73.6, 0.25, 0.35, -1.2],
            note="Flat ordinary SU2 diagnostic; Ol = 1 - Om - Omega_chi0.",
            to_cosmo=su2,
        ),
        ModelSpec(
            name="SU2R",
            param_names=["H0", "Om", "su2_fraction", "w0_chi"],
            bounds=[(45.0, 95.0), (0.05, 0.60), (0.0, 1.0), (-3.0, -0.20)],
            initial=[73.6, 0.25, 0.75, -1.05],
            note="Flat reparameterized SU2R diagnostic; samples SU2 fraction of the dark sector.",
            to_cosmo=su2r,
        ),
    ]


def cosmo_is_valid(cosmo: dict[str, float]) -> bool:
    H0 = cosmo["H0"]
    Om = cosmo["Om"]
    Ol = cosmo["Ol"]
    Omega_chi0 = cosmo["Omega_chi0"]
    if H0 <= 0.0 or Om <= 0.0 or Ol < -1e-8 or Omega_chi0 < -1e-8:
        return False
    if abs(1.0 - Om - Ol - Omega_chi0) > 1e-5:
        return False
    return True


def dl_model(z: np.ndarray, cosmo: dict[str, float], nint: int) -> np.ndarray:
    return luminosity_distance(
        z,
        cosmo["H0"],
        cosmo["Om"],
        cosmo["Ol"],
        0.0,
        0.0,
        0.0,
        nint_base=nint,
        Omega_chi0=cosmo.get("Omega_chi0", 0.0),
        w0_chi=cosmo.get("w0_chi", -1.0),
        wa_chi=cosmo.get("wa_chi", 0.0),
    )


def chi2_for_cosmo(events: list[EventRow], cosmo: dict[str, float], nint: int) -> float:
    if not cosmo_is_valid(cosmo):
        return float("inf")
    z = np.array([event.z for event in events], dtype=float)
    dl_obs = np.array([event.dl_mpc for event in events], dtype=float)
    sigma = np.array([event.sigma_ln_dl for event in events], dtype=float)
    dl_th = dl_model(z, cosmo, nint)
    if not np.all(np.isfinite(dl_th)) or np.any(dl_th <= 0.0):
        return float("inf")
    resid = (np.log(dl_th) - np.log(dl_obs)) / sigma
    return float(np.sum(resid * resid))


def fit_model(
    events: list[EventRow],
    spec: ModelSpec,
    nint: int,
    seed: int,
    maxiter: int,
    popsize: int,
) -> dict[str, Any]:
    if differential_evolution is None or minimize is None:
        raise RuntimeError(f"scipy optimize is unavailable: {SCIPY_IMPORT_ERROR}")

    def objective(theta: np.ndarray) -> float:
        return chi2_for_cosmo(events, spec.to_cosmo(theta), nint)

    result_de = differential_evolution(
        objective,
        spec.bounds,
        seed=seed,
        maxiter=maxiter,
        popsize=popsize,
        polish=False,
        tol=5e-5,
        updating="immediate",
        workers=1,
    )
    result_min = minimize(objective, result_de.x, method="Nelder-Mead", options={"maxiter": 600})
    theta = result_min.x if result_min.success and result_min.fun <= result_de.fun else result_de.x
    chi2 = float(objective(theta))
    k = len(spec.param_names)
    n = len(events)
    cosmo = spec.to_cosmo(theta)
    row: dict[str, Any] = {
        "source": "GWTC5_O4b_summary_fit",
        "model": spec.name,
        "N": n,
        "k": k,
        "dof": n - k,
        "chi2": chi2,
        "chi2_dof": chi2 / max(n - k, 1),
        "AIC": chi2 + 2.0 * k,
        "BIC": chi2 + k * math.log(max(n, 1)),
        "note": spec.note,
    }
    for name, value in zip(spec.param_names, theta):
        row[name] = float(value)
    for key in ["H0", "Om", "Ol", "Omega_chi0", "Omega_dark0", "su2_fraction", "w0_chi", "wa_chi"]:
        if key in cosmo:
            row[key] = float(cosmo[key])
    return row


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


def score_bestfits(events: list[EventRow], nint: int) -> list[dict[str, Any]]:
    candidates = [
        (
            "LCDM_DESI_DR2_SN_BAO_Planck",
            "LCDM",
            ROOT / "plamb_runs" / "updated_datasets_20260710" / "runs" / "desi_LCDM_Pantheon_SN_BAO_Planck" / "LCDM_emcee_bestfit.txt",
            "Current DESI DR2 matched LCDM scalar best fit.",
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
            "Corrected ordinary SU2 posterior rerun; may be absent while background chain is running.",
        ),
        (
            "SU2_DESI_DR2_legacy",
            "SU2",
            ROOT / "plamb_runs" / "updated_datasets_20260710" / "runs" / "desi_SU2_Pantheon_SN_BAO_Planck" / "SU2_emcee_bestfit.txt",
            "Legacy ordinary SU2 scalar fit; use only as a reference if fixed-posterior fit is absent.",
        ),
        (
            "SU2R_DESI_DR2_SN_BAO_Planck",
            "SU2R",
            ROOT / "plamb_runs" / "updated_datasets_20260710" / "runs" / "desi_SU2R_Pantheon_SN_BAO_Planck" / "SU2R_emcee_bestfit.txt",
            "Current DESI DR2 matched SU2R scalar best fit.",
        ),
    ]
    rows: list[dict[str, Any]] = []
    for label, model, path, note in candidates:
        best = read_bestfit(path)
        if not best:
            rows.append({"source": label, "model": model, "available": "no", "path": str(path), "note": note})
            continue
        cosmo = {
            "H0": float(best.get("H0", float("nan"))),
            "Om": float(best.get("Om", float("nan"))),
            "Ol": float(best.get("Ol", float("nan"))),
            "Omega_chi0": float(best.get("Omega_chi0", 0.0)),
            "w0_chi": float(best.get("w0_chi", -1.0)),
            "wa_chi": float(best.get("wa_chi", 0.0)),
        }
        chi2 = chi2_for_cosmo(events, cosmo, nint)
        active = best.get("active", [])
        k = len(active) if isinstance(active, list) else 0
        n = len(events)
        row = {
            "source": label,
            "model": model,
            "available": "yes",
            "path": str(path),
            "N": n,
            "k_scalar_active": k,
            "chi2_gw": chi2,
            "chi2_dof_gw": chi2 / max(n, 1),
            "AIC_gw_scored": chi2 + 2.0 * k,
            "BIC_gw_scored": chi2 + k * math.log(max(n, 1)),
            "H0": cosmo["H0"],
            "Om": cosmo["Om"],
            "Ol": cosmo["Ol"],
            "Omega_chi0": cosmo["Omega_chi0"],
            "w0_chi": cosmo["w0_chi"],
            "wa_chi": cosmo["wa_chi"],
            "note": note,
        }
        if "su2_fraction" in best:
            row["su2_fraction"] = best.get("su2_fraction")
        rows.append(row)
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


def residual_rows(events: list[EventRow], fit_rows: list[dict[str, Any]], nint: int) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    specs_by_name = {spec.name: spec for spec in model_specs()}
    for fit in fit_rows:
        spec = specs_by_name[fit["model"]]
        theta = np.array([float(fit[name]) for name in spec.param_names], dtype=float)
        cosmo = spec.to_cosmo(theta)
        z = np.array([event.z for event in events], dtype=float)
        dl_th = dl_model(z, cosmo, nint)
        for event, model_dl in zip(events, dl_th):
            rows.append(
                {
                    "model": fit["model"],
                    "event_id": event.event_id,
                    "utc_approx": event.utc_approx,
                    "z": event.z,
                    "dl_obs_mpc": event.dl_mpc,
                    "dl_model_mpc": float(model_dl),
                    "sigma_ln_dl": event.sigma_ln_dl,
                    "residual_ln_dl": float(math.log(model_dl) - math.log(event.dl_mpc)),
                    "pull_ln_dl": float((math.log(model_dl) - math.log(event.dl_mpc)) / event.sigma_ln_dl),
                    "p_astro": event.p_astro,
                    "network_snr": event.network_snr,
                }
            )
    return rows


def write_plot(path: Path, events: list[EventRow], fit_rows: list[dict[str, Any]], nint: int) -> str:
    if plt is None:
        return f"Plot skipped because matplotlib is unavailable: {MATPLOTLIB_IMPORT_ERROR}"
    specs_by_name = {spec.name: spec for spec in model_specs()}
    z = np.array([event.z for event in events], dtype=float)
    dl = np.array([event.dl_mpc for event in events], dtype=float)
    yerr_low = np.array([abs(event.dl_lower) for event in events], dtype=float)
    yerr_high = np.array([abs(event.dl_upper) for event in events], dtype=float)
    z_grid = np.linspace(max(1e-4, float(np.min(z)) * 0.7), float(np.max(z)) * 1.08, 240)

    fig, ax = plt.subplots(figsize=(8.8, 5.6), dpi=150)
    ax.errorbar(z, dl, yerr=[yerr_low, yerr_high], fmt=".", ms=4, alpha=0.45, color="0.25", label="GWTC-5 O4b summary rows")
    colors = {"LCDM": "#1f77b4", "SU2": "#d62728", "SU2R": "#2ca02c"}
    for fit in fit_rows:
        spec = specs_by_name[fit["model"]]
        theta = np.array([float(fit[name]) for name in spec.param_names], dtype=float)
        cosmo = spec.to_cosmo(theta)
        ax.plot(z_grid, dl_model(z_grid, cosmo, nint), lw=2, label=f"{fit['model']} GW-only fit", color=colors.get(fit["model"]))
    ax.set_xlabel("GWTC redshift summary")
    ax.set_ylabel("Luminosity distance summary [Mpc]")
    ax.set_title("GWTC-5.0/O4b diagnostic distance-redshift comparison")
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)
    return f"Saved plot: {path}"


def markdown_table(rows: list[dict[str, Any]], columns: list[str]) -> list[str]:
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in rows:
        values: list[str] = []
        for col in columns:
            value = row.get(col, "")
            if isinstance(value, float):
                if math.isfinite(value):
                    values.append(f"{value:.6g}")
                else:
                    values.append("")
            else:
                values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    return lines


def write_report(
    path: Path,
    events: list[EventRow],
    fit_rows: list[dict[str, Any]],
    scored_rows: list[dict[str, Any]],
    args: argparse.Namespace,
    plot_status: str,
) -> None:
    best = min(fit_rows, key=lambda row: row["AIC"]) if fit_rows else {}
    scored_available = [row for row in scored_rows if row.get("available") == "yes"]
    lines = [
        "# SU2 / GWTC-5.0 O4b Diagnostic",
        "",
        "## Headline",
        "",
        f"- Usable GWTC-5.0/O4b rows: `{len(events)}` with redshift and luminosity-distance summary values.",
        f"- Minimum p_astro cut: `{args.min_pastro}`.",
        f"- Distance likelihood: Gaussian in `ln(D_L)` using compressed asymmetric GWOSC distance uncertainty, floor `{args.sigma_floor}`.",
        f"- Best GW-only AIC row: `{best.get('model', 'n/a')}`.",
        "",
        "This is a diagnostic cross-check, not a publication-grade GW likelihood. The GWTC redshift values are summary products tied to GW parameter inference and cosmological assumptions; they are not equivalent to independent host-galaxy redshifts for all events.",
        "",
        "## GW-Only Fits",
        "",
    ]
    lines.extend(markdown_table(fit_rows, ["model", "N", "k", "chi2", "chi2_dof", "AIC", "BIC", "H0", "Om", "Ol", "Omega_chi0", "su2_fraction", "w0_chi"]))
    lines.extend(
        [
            "",
            "## Existing Scalar Best Fits Scored On GWTC-5.0",
            "",
        ]
    )
    lines.extend(markdown_table(scored_available, ["source", "model", "chi2_gw", "chi2_dof_gw", "AIC_gw_scored", "BIC_gw_scored", "H0", "Om", "Omega_chi0", "su2_fraction", "w0_chi"]))
    missing = [row for row in scored_rows if row.get("available") != "yes"]
    if missing:
        lines.extend(["", "## Missing Scored Fits", ""])
        lines.extend(markdown_table(missing, ["source", "model", "available", "note", "path"]))
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- This run can be included as a GW diagnostic branch in the SU2 program.",
            "- It should not yet move the core SU2 posterior because it ignores the full GW selection function, event posterior correlations, merger-rate population, and sky-localization information.",
            "- If the GW-only fit strongly prefers parameter boundaries, that is more likely a sign of weak or summary-level data than decisive cosmology.",
            "- A proper next step is a hierarchical GW likelihood using posterior samples/skymaps for selected O4b events, not only the catalog summary columns.",
            "",
            "## Local Outputs",
            "",
            f"- Fit summary: `{OUTDIR / 'su2_gwtc5_o4b_fit_summary.csv'}`",
            f"- Scored best fits: `{OUTDIR / 'su2_gwtc5_o4b_scored_bestfits.csv'}`",
            f"- Event residuals: `{OUTDIR / 'su2_gwtc5_o4b_event_residuals.csv'}`",
            f"- Plot: `{OUTDIR / 'su2_gwtc5_o4b_hubble.png'}`",
            f"- Plot status: {plot_status}",
            "",
            "## Recommended Program Update",
            "",
            "1. Add this diagnostic to the SU2 report as an independent GW cross-probe.",
            "2. Rerun after the corrected ordinary SU2 DESI posterior completes, so the scored-bestfit table includes the fixed SU2 result.",
            "3. Build a `likelihood_gwtc5_o4b.py` layer only after selecting events with usable posterior samples and defining the population/selection model.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--event-csv", type=Path, default=DEFAULT_EVENT_CSV)
    parser.add_argument("--out", type=Path, default=OUTDIR)
    parser.add_argument("--min-pastro", type=float, default=0.5)
    parser.add_argument("--sigma-floor", type=float, default=0.05)
    parser.add_argument("--nint", type=int, default=140)
    parser.add_argument("--seed", type=int, default=25050)
    parser.add_argument("--maxiter", type=int, default=80, help="Differential evolution max iterations per model.")
    parser.add_argument("--popsize", type=int, default=7, help="Differential evolution population multiplier.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    if not args.event_csv.exists():
        print(f"Missing GWTC event CSV: {args.event_csv}", file=sys.stderr)
        return 2

    events = load_events(args.event_csv, args.min_pastro, args.sigma_floor)
    if not events:
        print("No usable GWTC events after filters.", file=sys.stderr)
        return 2

    print(f"Loaded {len(events)} GWTC-5.0/O4b events with redshift and D_L summaries.")
    fit_rows = []
    for index, spec in enumerate(model_specs()):
        print(f"[fit] {spec.name}", flush=True)
        fit_rows.append(fit_model(events, spec, args.nint, args.seed + index, args.maxiter, args.popsize))

    scored_rows = score_bestfits(events, args.nint)
    residuals = residual_rows(events, fit_rows, args.nint)

    fit_csv = args.out / "su2_gwtc5_o4b_fit_summary.csv"
    scored_csv = args.out / "su2_gwtc5_o4b_scored_bestfits.csv"
    residual_csv = args.out / "su2_gwtc5_o4b_event_residuals.csv"
    report_md = args.out / "su2_gwtc5_o4b_report.md"
    plot_png = args.out / "su2_gwtc5_o4b_hubble.png"

    write_csv(fit_csv, fit_rows)
    write_csv(scored_csv, scored_rows)
    write_csv(residual_csv, residuals)
    plot_status = write_plot(plot_png, events, fit_rows, args.nint)
    write_report(report_md, events, fit_rows, scored_rows, args, plot_status)

    best = min(fit_rows, key=lambda row: row["AIC"])
    print(f"Saved fit summary: {fit_csv}")
    print(f"Saved scored bestfits: {scored_csv}")
    print(f"Saved residuals: {residual_csv}")
    print(f"Saved report: {report_md}")
    print(plot_status)
    print(
        "Best GW-only diagnostic row: "
        f"{best['model']} chi2={best['chi2']:.3f} AIC={best['AIC']:.3f} "
        f"BIC={best['BIC']:.3f} H0={best.get('H0', float('nan')):.3f}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
