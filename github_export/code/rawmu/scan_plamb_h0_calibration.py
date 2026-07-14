r"""Calibration-constrained H0 sensitivity scan for PLAMB exact log kernel.

The local-flow diagnostics normally fit a free global magnitude offset.  That
offset absorbs nearly all H0 movement in SN-only raw-MU fits:

    MU_model = 5 log10[(c/H0) D_shape(z)] + 25 + offset

This diagnostic makes the calibration assumption explicit by scanning H0 while
placing a Gaussian prior on the offset.  It compares:

    * EXACT_LOG_Pfree: PLAMB exact logarithmic propagation kernel
    * LCDM_Omfree: flat LCDM with Om fitted

Offset prior modes:

    fixed       offset fixed to offset_mean
    numeric     offset Gaussian prior with sigma in magnitudes
    free        no calibration prior; H0 should become degenerate

Outputs identify how preferred H0 and PLAMB-vs-LCDM ranking change as the
absolute calibration is tightened or relaxed.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from scipy.optimize import minimize_scalar

    HAVE_SCIPY = True
except Exception:  # pragma: no cover
    minimize_scalar = None
    HAVE_SCIPY = False

from compare_plamb_rawmu_exact_kernels import full_params, mu_model
from diagnose_pantheon_rawmu_fr import read_table
from diagnose_plamb_log_kernel_residuals import (
    DEFAULT_DES,
    DEFAULT_DES_COV,
    DEFAULT_PANTHEON,
    DatasetConfig,
    dataset_configs,
    numeric,
    prepare_dataset,
    spec_by_name,
)


ROOT = Path(__file__).resolve().parent
DEFAULT_OUTDIR = ROOT / "plamb_runs" / "diagnostics" / "plamb_h0_calibration_scan"
DEFAULT_PANTHEON_COV = ROOT / "Pantheon+SH0ES_STAT+SYS.cov"


@dataclass(frozen=True)
class OffsetPrior:
    label: str
    sigma: float | None
    fixed: bool


def parse_h0_grid(raw: str) -> np.ndarray:
    text = str(raw).strip()
    if ":" in text:
        lo, hi, step = [float(part) for part in text.split(":", 2)]
        n = int(round((hi - lo) / step)) + 1
        return lo + step * np.arange(n, dtype=float)
    values = [float(item.strip()) for item in text.split(",") if item.strip()]
    return np.asarray(values, dtype=float)


def parse_float_list(raw: str) -> list[float]:
    out: list[float] = []
    for item in str(raw).split(","):
        item = item.strip()
        if item:
            out.append(float(item))
    return out


def parse_offset_priors(raw: str) -> list[OffsetPrior]:
    priors: list[OffsetPrior] = []
    for item in str(raw).split(","):
        token = item.strip().lower()
        if not token:
            continue
        if token in {"fixed", "fix", "zero"}:
            priors.append(OffsetPrior("fixed", 0.0, True))
        elif token in {"free", "none", "inf", "infinite"}:
            priors.append(OffsetPrior("free", None, False))
        else:
            sigma = float(token)
            if sigma <= 0.0:
                priors.append(OffsetPrior("fixed", 0.0, True))
            else:
                label = f"sigma_{sigma:g}".replace(".", "p")
                priors.append(OffsetPrior(label, sigma, False))
    return priors


def requested_z_columns(raw_df: pd.DataFrame, requested: str) -> list[str]:
    if requested.lower() == "all":
        return [col for col in ["zHD", "zCMB", "zHEL"] if col in raw_df.columns]
    return [col.strip() for col in requested.split(",") if col.strip() and col.strip() in raw_df.columns]


def make_args(base: argparse.Namespace, **updates: object) -> SimpleNamespace:
    values = vars(base).copy()
    values.update(updates)
    return SimpleNamespace(**values)


def solve_offset_with_prior(
    residual_base: np.ndarray,
    precision: np.ndarray,
    prior: OffsetPrior,
    offset_mean: float,
) -> tuple[float, float, float, float]:
    """Solve profiled scalar offset and return offset, data chi2, prior chi2, total chi2."""
    r = np.asarray(residual_base, dtype=float)
    one = np.ones_like(r)
    if prior.fixed:
        offset = float(offset_mean)
    else:
        denom = float(one @ precision @ one)
        numer = float(one @ precision @ r)
        if prior.sigma is not None:
            inv_var = 1.0 / (prior.sigma * prior.sigma)
            denom += inv_var
            numer += offset_mean * inv_var
        offset = numer / denom
    data_resid = r - offset
    data_chi2 = float(data_resid @ precision @ data_resid)
    if prior.sigma is None:
        prior_chi2 = 0.0
    elif prior.fixed:
        prior_chi2 = 0.0
    else:
        prior_chi2 = float(((offset - offset_mean) / prior.sigma) ** 2)
    return offset, data_chi2, prior_chi2, data_chi2 + prior_chi2


def score_params(
    z: np.ndarray,
    mu: np.ndarray,
    precision: np.ndarray,
    spec,
    h0: float,
    params: dict[str, float],
    prior: OffsetPrior,
    offset_mean: float,
) -> tuple[float, float, float, float]:
    model = mu_model(z, spec, h0, params)
    if not np.all(np.isfinite(model)):
        return float("nan"), float("inf"), float("inf"), float("inf")
    return solve_offset_with_prior(mu - model, precision, prior, offset_mean)


def objective_value(
    theta: float,
    z: np.ndarray,
    mu: np.ndarray,
    precision: np.ndarray,
    spec,
    h0: float,
    prior: OffsetPrior,
    offset_mean: float,
) -> float:
    params = full_params(spec, [float(theta)])
    _offset, _data, _prior, total = score_params(z, mu, precision, spec, h0, params, prior, offset_mean)
    return total if np.isfinite(total) else 1.0e100


def fit_at_h0(
    z: np.ndarray,
    mu: np.ndarray,
    precision: np.ndarray,
    spec,
    h0: float,
    prior: OffsetPrior,
    offset_mean: float,
    grid_steps: int,
    maxiter: int,
) -> dict[str, object]:
    if len(spec.free_params) != 1:
        raise ValueError(f"This scan expects one free shape parameter for {spec.name}.")
    name = spec.free_params[0]
    lo, hi = spec.bounds[name]
    if HAVE_SCIPY:
        opt = minimize_scalar(
            lambda value: objective_value(float(value), z, mu, precision, spec, h0, prior, offset_mean),
            bounds=(lo, hi),
            method="bounded",
            options={"maxiter": maxiter},
        )
        theta = float(opt.x)
        success = bool(opt.success)
        method = f"minimize_{name}_with_offset_prior"
    else:
        grid = np.linspace(lo, hi, int(grid_steps))
        vals = np.array([objective_value(v, z, mu, precision, spec, h0, prior, offset_mean) for v in grid])
        theta = float(grid[int(np.nanargmin(vals))])
        success = bool(np.isfinite(vals).any())
        method = f"grid_{name}_with_offset_prior"
    params = full_params(spec, [theta])
    offset, data_chi2, prior_chi2, total_chi2 = score_params(z, mu, precision, spec, h0, params, prior, offset_mean)
    return {
        "model": spec.name,
        "family": spec.family,
        "param_name": name,
        "param_value": theta,
        "H0": h0,
        "offset_mag": offset,
        "offset_data_chi2": data_chi2,
        "offset_prior_chi2": prior_chi2,
        "chi2_total": total_chi2,
        "method": method,
        "success": bool(success and np.isfinite(total_chi2)),
    }


def run_scan(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    h0_grid = parse_h0_grid(args.h0_grid)
    priors = parse_offset_priors(args.offset_priors)
    zmins = parse_float_list(args.z_min_values)
    specs = [spec_by_name("EXACT_LOG_Pfree"), spec_by_name("LCDM_Omfree")]
    rows: list[dict[str, object]] = []
    errors: list[dict[str, object]] = []
    fit_count = 0

    for cfg0 in dataset_configs(args):
        cfg = cfg0
        if cfg0.label == "PantheonPlusSH0ES":
            cfg = DatasetConfig(
                args.pantheon_label,
                cfg0.data_path,
                cfg0.cov_path,
                args.pantheon_mu_col,
                cfg0.err_col,
                calibrator_col=cfg0.calibrator_col,
                exclude_calibrators=cfg0.exclude_calibrators,
            )
        raw = read_table(cfg.data_path)
        for z_col in requested_z_columns(raw, args.z_col):
            for z_min in zmins:
                prep_args = make_args(args, z_min=z_min)
                try:
                    prepared = prepare_dataset(cfg, z_col, prep_args)
                except Exception as exc:
                    errors.append({"dataset": cfg.label, "z_col": z_col, "z_min": z_min, "note": str(exc)})
                    continue
                for prior in priors:
                    for h0 in h0_grid:
                        for spec in specs:
                            try:
                                fit = fit_at_h0(
                                    prepared.z,
                                    prepared.mu,
                                    prepared.precision,
                                    spec,
                                    float(h0),
                                    prior,
                                    args.offset_mean,
                                    args.grid_steps,
                                    args.maxiter,
                                )
                                k_offset = 0 if prior.fixed else 1
                                k_total = int(spec.k_shape) + k_offset
                                rows.append(
                                    {
                                        "dataset": cfg.label,
                                        "z_col": z_col,
                                        "z_min": z_min,
                                        "N": int(len(prepared.z)),
                                        "offset_prior": prior.label,
                                        "offset_sigma": np.nan if prior.sigma is None else prior.sigma,
                                        "offset_fixed": prior.fixed,
                                        "offset_mean": args.offset_mean,
                                        "H0": float(h0),
                                        "model": fit["model"],
                                        "family": fit["family"],
                                        "param_name": fit["param_name"],
                                        "param_value": fit["param_value"],
                                        "offset_mag": fit["offset_mag"],
                                        "offset_data_chi2": fit["offset_data_chi2"],
                                        "offset_prior_chi2": fit["offset_prior_chi2"],
                                        "chi2_total": fit["chi2_total"],
                                        "AIC_like": float(fit["chi2_total"]) + 2.0 * k_total,
                                        "BIC_like": float(fit["chi2_total"]) + math.log(max(len(prepared.z), 2)) * k_total,
                                        "k_total_fixed_H0": k_total,
                                        "success": fit["success"],
                                        "method": fit["method"],
                                        "cov_note": prepared.cov_note,
                                        "subset_note": prepared.subset_note,
                                    }
                                )
                                fit_count += 1
                                if args.progress_every and fit_count % int(args.progress_every) == 0:
                                    print(
                                        f"[fit {fit_count}] {cfg.label} {z_col} zmin={z_min:g} "
                                        f"offset={prior.label} H0={float(h0):.2f} model={spec.name} "
                                        f"chi2={float(fit['chi2_total']):.4g}",
                                        flush=True,
                                    )
                            except Exception as exc:
                                errors.append(
                                    {
                                        "dataset": cfg.label,
                                        "z_col": z_col,
                                        "z_min": z_min,
                                        "offset_prior": prior.label,
                                        "H0": h0,
                                        "model": spec.name,
                                        "note": str(exc),
                                    }
                                )
    grid = pd.DataFrame(rows)
    errors_df = pd.DataFrame(errors)
    compare = build_comparison(grid)
    return grid, compare, errors_df


def build_comparison(grid: pd.DataFrame) -> pd.DataFrame:
    if grid.empty:
        return pd.DataFrame()
    key_cols = ["dataset", "z_col", "z_min", "N", "offset_prior", "offset_sigma", "offset_fixed", "offset_mean", "H0"]
    log = grid[grid["model"] == "EXACT_LOG_Pfree"].copy()
    lcdm = grid[grid["model"] == "LCDM_Omfree"].copy()
    cols = [
        "param_value",
        "offset_mag",
        "offset_data_chi2",
        "offset_prior_chi2",
        "chi2_total",
        "AIC_like",
        "BIC_like",
        "success",
    ]
    log = log[key_cols + cols].rename(columns={c: f"{c}_log" for c in cols})
    lcdm = lcdm[key_cols + cols].rename(columns={c: f"{c}_lcdm" for c in cols})
    out = log.merge(lcdm, on=key_cols, how="inner")
    out["delta_chi2_log_minus_lcdm"] = out["chi2_total_log"] - out["chi2_total_lcdm"]
    out["delta_AIC_log_minus_lcdm"] = out["AIC_like_log"] - out["AIC_like_lcdm"]
    out["delta_BIC_log_minus_lcdm"] = out["BIC_like_log"] - out["BIC_like_lcdm"]
    out["preferred_AIC"] = np.where(out["delta_AIC_log_minus_lcdm"] < 0, "EXACT_LOG_Pfree", "LCDM_Omfree")
    return out


def build_profile_summary(compare: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if compare.empty:
        return pd.DataFrame()
    group_cols = ["dataset", "z_col", "z_min", "N", "offset_prior", "offset_sigma", "offset_fixed", "offset_mean"]
    for key, grp in compare.groupby(group_cols, dropna=False):
        row = {col: value for col, value in zip(group_cols, key)}
        log_best = grp.loc[grp["chi2_total_log"].idxmin()]
        lcdm_best = grp.loc[grp["chi2_total_lcdm"].idxmin()]
        delta_profile = float(log_best["chi2_total_log"] - lcdm_best["chi2_total_lcdm"])
        row.update(
            {
                "best_H0_log": float(log_best["H0"]),
                "best_p_log": float(log_best["param_value_log"]),
                "best_offset_log": float(log_best["offset_mag_log"]),
                "best_chi2_log": float(log_best["chi2_total_log"]),
                "best_H0_lcdm": float(lcdm_best["H0"]),
                "best_Om_lcdm": float(lcdm_best["param_value_lcdm"]),
                "best_offset_lcdm": float(lcdm_best["offset_mag_lcdm"]),
                "best_chi2_lcdm": float(lcdm_best["chi2_total_lcdm"]),
                "delta_profile_chi2_log_minus_lcdm": delta_profile,
                "delta_profile_AIC_log_minus_lcdm": delta_profile,
                "profile_preferred": "EXACT_LOG_Pfree" if delta_profile < 0.0 else "LCDM_Omfree",
            }
        )
        at_675 = grp.iloc[(grp["H0"].astype(float) - 67.5).abs().argsort()[:1]]
        if not at_675.empty:
            r675 = at_675.iloc[0]
            row.update(
                {
                    "nearest_675_H0": float(r675["H0"]),
                    "delta_AIC_at_675_log_minus_lcdm": float(r675["delta_AIC_log_minus_lcdm"]),
                    "preferred_at_675": str(r675["preferred_AIC"]),
                }
            )
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["dataset", "z_col", "z_min", "offset_prior"]).reset_index(drop=True)


def plot_curves(compare: pd.DataFrame, outdir: Path) -> list[Path]:
    paths: list[Path] = []
    if compare.empty:
        return paths
    plot_dir = outdir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    for (dataset, z_col, z_min, prior), grp in compare.groupby(["dataset", "z_col", "z_min", "offset_prior"], dropna=False):
        grp = grp.sort_values("H0")
        fig, ax = plt.subplots(figsize=(8.5, 5))
        ax.plot(grp["H0"], grp["delta_AIC_log_minus_lcdm"], marker="o", ms=2.5, lw=1.4)
        ax.axhline(0.0, color="black", lw=0.8)
        ax.axvline(67.5, color="#777777", lw=0.8, ls="--")
        ax.set_xlabel("H0")
        ax.set_ylabel("delta AIC: exact log - LCDM")
        ax.set_title(f"{dataset} {z_col} z_min={z_min:g} offset={prior}")
        ax.grid(alpha=0.25)
        safe = "".join(c if c.isalnum() or c in {"_", "-"} else "_" for c in f"{dataset}_{z_col}_{z_min}_{prior}")
        path = plot_dir / f"h0_delta_aic_{safe}.png"
        fig.tight_layout()
        fig.savefig(path, dpi=150)
        plt.close(fig)
        paths.append(path)
    return paths


def write_report(path: Path, profile: pd.DataFrame, compare: pd.DataFrame, plots: list[Path], args: argparse.Namespace) -> None:
    lines: list[str] = []
    lines.append("# PLAMB H0 / Calibration Sensitivity Scan")
    lines.append("")
    lines.append("This scan constrains or fixes the global magnitude offset so H0 can no longer be completely absorbed by calibration.")
    lines.append("")
    lines.append("## Configuration")
    lines.append("")
    lines.append(f"- datasets: `{args.dataset}`")
    lines.append(f"- z columns: `{args.z_col}`")
    lines.append(f"- z_min values: `{args.z_min_values}`")
    lines.append(f"- H0 grid: `{args.h0_grid}`")
    lines.append(f"- offset priors: `{args.offset_priors}`")
    lines.append(f"- offset prior mean: `{args.offset_mean}` mag")
    lines.append("")
    lines.append("## Profiled H0 Summary")
    lines.append("")
    if profile.empty:
        lines.append("_No profile rows._")
    else:
        cols = [
            "dataset",
            "z_col",
            "z_min",
            "offset_prior",
            "best_H0_log",
            "best_H0_lcdm",
            "best_offset_log",
            "best_offset_lcdm",
            "delta_profile_AIC_log_minus_lcdm",
            "profile_preferred",
            "delta_AIC_at_675_log_minus_lcdm",
            "preferred_at_675",
        ]
        lines.append(profile[cols].to_markdown(index=False, floatfmt=".6g"))
    lines.append("")
    lines.append("## Interpretation")
    lines.append("")
    lines.append("- In the `free` offset mode, H0 is expected to be nearly degenerate with the offset; flat chi2 curves are a calibration-degeneracy check.")
    lines.append("- Tight offset priors or fixed offset make H0 identifiable, but they also assume the published MU absolute calibration is directly usable with zero additional offset.")
    lines.append("- If PLAMB/log only prefers H0=67.5 under very tight offset assumptions, that is a calibration-dependent result, not a pure shape result.")
    lines.append("- If the same model remains preferred across loose priors and corrected redshift bases, that would be a stronger cosmological signal.")
    lines.append("")
    lines.append("## Plots")
    lines.append("")
    for plot in plots:
        lines.append(f"- `{plot}`")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pantheon-data", type=Path, default=DEFAULT_PANTHEON)
    parser.add_argument("--pantheon-cov", type=Path, default=DEFAULT_PANTHEON_COV)
    parser.add_argument("--pantheon-mu-col", default="MU_SH0ES")
    parser.add_argument("--pantheon-label", default="PantheonPlusSH0ES")
    parser.add_argument("--des-data", type=Path, default=DEFAULT_DES)
    parser.add_argument("--des-cov", type=Path, default=DEFAULT_DES_COV)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--dataset", default="all")
    parser.add_argument("--z-col", default="all")
    parser.add_argument("--z-min-values", default="0.01,0.025,0.05")
    parser.add_argument("--H0-value", type=float, default=67.5, help="Compatibility placeholder; scan uses --h0-grid.")
    parser.add_argument("--h0-grid", default="60:80:0.5")
    parser.add_argument("--offset-priors", default="fixed,0.02,0.05,0.10,0.20,free")
    parser.add_argument("--offset-mean", type=float, default=0.0)
    parser.add_argument("--max-z", type=float, default=None)
    parser.add_argument("--min-n", type=int, default=10)
    parser.add_argument("--keep-calibrators", action="store_true")
    parser.add_argument("--allow-precision-submatrix", action="store_true")
    parser.add_argument("--grid-steps", type=int, default=401)
    parser.add_argument("--maxiter", type=int, default=100)
    parser.add_argument("--progress-every", type=int, default=250)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    grid, compare, errors = run_scan(args)
    profile = build_profile_summary(compare)
    plots = plot_curves(compare, args.outdir)

    grid_path = args.outdir / "plamb_h0_calibration_grid.csv"
    compare_path = args.outdir / "plamb_h0_calibration_comparison.csv"
    profile_path = args.outdir / "plamb_h0_calibration_profile_summary.csv"
    errors_path = args.outdir / "plamb_h0_calibration_errors.csv"
    report_path = args.outdir / "plamb_h0_calibration_report.md"
    config_path = args.outdir / "plamb_h0_calibration_config.json"
    grid.to_csv(grid_path, index=False)
    compare.to_csv(compare_path, index=False)
    profile.to_csv(profile_path, index=False)
    errors.to_csv(errors_path, index=False)
    config_path.write_text(json.dumps(vars(args), indent=2, default=str), encoding="utf-8")
    write_report(report_path, profile, compare, plots, args)

    print(f"Saved grid: {grid_path}")
    print(f"Saved comparison: {compare_path}")
    print(f"Saved profile summary: {profile_path}")
    print(f"Saved report: {report_path}")
    print(f"Saved plots: {args.outdir / 'plots'}")
    if not profile.empty:
        show = profile.sort_values("delta_profile_AIC_log_minus_lcdm").head(10)
        for _, row in show.iterrows():
            print(
                f"{row['dataset']} {row['z_col']} zmin={float(row['z_min']):g} "
                f"offset={row['offset_prior']} bestH0(log/LCDM)="
                f"{float(row['best_H0_log']):.2f}/{float(row['best_H0_lcdm']):.2f} "
                f"dAIC_profile={float(row['delta_profile_AIC_log_minus_lcdm']):.3f} "
                f"winner={row['profile_preferred']}"
            )


if __name__ == "__main__":
    main()
