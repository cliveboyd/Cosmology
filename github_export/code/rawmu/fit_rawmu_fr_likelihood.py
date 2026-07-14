r"""
SN-only likelihood fits for Peter Lamb's raw-MU FR/no-expansion proposal.

This promotes the raw distance-modulus diagnostic into a simple model-comparison
likelihood without touching `PLamb_Test_10V6c_plus.py`.

The model is deliberately minimal:

    D_obs_model(z) = (c / H0) * z * (1 + beta z)
    MU_model      = 5 log10(D_obs_model/Mpc) + 25

Peter's current proposal corresponds approximately to `beta = 0.5`, because
the empirical distance from MU is divided by `(1 + z/2)` before testing the
linear Hubble-time relation.

The script fits SN-only chi2/AIC/BIC for:

    * RAW_LINEAR_H0free       beta=0, H0 free
    * FR_BETA05_H0free        beta=0.5, H0 free
    * FR_BETAfree_H0free      beta free, H0 free
    * RAW_LINEAR_H0fixed675   beta=0, H0 fixed at 67.5
    * FR_BETA05_H0fixed675    beta=0.5, H0 fixed at 67.5
    * FR_BETAfree_H0fixed675  beta free, H0 fixed at 67.5

Outputs:
    plamb_runs/diagnostics/rawmu_fr_likelihood/rawmu_fr_likelihood_summary.csv
    plamb_runs/diagnostics/rawmu_fr_likelihood/rawmu_fr_likelihood_key_results.csv
    plamb_runs/diagnostics/rawmu_fr_likelihood/rawmu_fr_likelihood_report.md
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from scipy.optimize import minimize, minimize_scalar

    HAVE_SCIPY = True
except Exception:  # pragma: no cover - fallback is only for lean installs.
    minimize = None
    minimize_scalar = None
    HAVE_SCIPY = False

from diagnose_pantheon_rawmu_fr import (
    CALIBRATOR_CANDIDATES,
    DEFAULT_PANTHEON,
    HF_CANDIDATES,
    MU_CANDIDATES,
    MU_ERR_CANDIDATES,
    PROB_CANDIDATES,
    C_KMS,
    DatasetSpec,
    auto_z_columns,
    find_column,
    numeric_series,
    read_table,
    sanitize_label,
    subset_masks,
)


ROOT = Path(__file__).resolve().parent
DEFAULT_OUTDIR = ROOT / "plamb_runs" / "diagnostics" / "rawmu_fr_likelihood"


DEFAULT_DATASETS = [
    DatasetSpec("PantheonPlusSH0ES", DEFAULT_PANTHEON),
    DatasetSpec(
        "DES_Dovekie_raw",
        ROOT
        / "external_datasets"
        / "current_cosmology_sources"
        / "DES_SN5YR_Dovekie"
        / "DES-Dovekie_HD.csv",
    ),
    DatasetSpec(
        "DES_Dovekie_PLamb",
        ROOT
        / "external_datasets"
        / "current_cosmology_sources"
        / "prepared_for_plamb"
        / "DES_SN5YR_Dovekie_PLamb.dat",
    ),
    DatasetSpec(
        "Union3_compressed",
        ROOT
        / "external_datasets"
        / "current_cosmology_sources"
        / "Union3_release"
        / "mu_mat_union3_cosmo=2_mu.fits",
    ),
    DatasetSpec(
        "Union3p1_UNITY1p8",
        ROOT
        / "external_datasets"
        / "current_cosmology_sources"
        / "Union3_release"
        / "mu_mat_union3.1_UNITY1.8_template_cosmo=2_0_mu.fits",
    ),
]


@dataclass(frozen=True)
class ModelSpec:
    name: str
    h0_fixed: float | None
    beta_fixed: float | None

    @property
    def k(self) -> int:
        return int(self.h0_fixed is None) + int(self.beta_fixed is None)


MODELS = [
    ModelSpec("RAW_LINEAR_H0free", None, 0.0),
    ModelSpec("FR_BETA05_H0free", None, 0.5),
    ModelSpec("FR_BETAfree_H0free", None, None),
    ModelSpec("RAW_LINEAR_H0fixed675", 67.5, 0.0),
    ModelSpec("FR_BETA05_H0fixed675", 67.5, 0.5),
    ModelSpec("FR_BETAfree_H0fixed675", 67.5, None),
]


def parse_dataset_arg(raw: str) -> DatasetSpec:
    if "=" in raw:
        label, path = raw.split("=", 1)
        return DatasetSpec(sanitize_label(label), Path(path))
    path = Path(raw)
    return DatasetSpec(sanitize_label(path.stem), path)


def mu_model_rawfr(z: np.ndarray, h0: float, beta: float) -> np.ndarray:
    d_model = (C_KMS / h0) * z * (1.0 + beta * z)
    if np.any(d_model <= 0.0) or not np.all(np.isfinite(d_model)):
        return np.full_like(z, np.inf, dtype=float)
    return 5.0 * np.log10(d_model) + 25.0


def chi2_mu(z: np.ndarray, mu: np.ndarray, sigma: np.ndarray, h0: float, beta: float) -> float:
    model = mu_model_rawfr(z, h0, beta)
    if not np.all(np.isfinite(model)):
        return float("inf")
    resid = mu - model
    return float(np.sum(np.square(resid / sigma)))


def fit_model(
    z: np.ndarray,
    mu: np.ndarray,
    sigma: np.ndarray,
    model: ModelSpec,
    h0_bounds: tuple[float, float],
    beta_bounds: tuple[float, float],
) -> dict[str, float | str | bool | int]:
    if not HAVE_SCIPY:
        return fit_model_grid(z, mu, sigma, model, h0_bounds, beta_bounds)

    if model.h0_fixed is not None and model.beta_fixed is not None:
        h0 = model.h0_fixed
        beta = model.beta_fixed
        chi2 = chi2_mu(z, mu, sigma, h0, beta)
        ok = np.isfinite(chi2)
        return {"H0": h0, "beta": beta, "chi2": chi2, "success": ok, "method": "fixed"}

    if model.h0_fixed is None and model.beta_fixed is not None:
        beta = model.beta_fixed

        def objective(h0: float) -> float:
            return chi2_mu(z, mu, sigma, h0, beta)

        opt = minimize_scalar(objective, bounds=h0_bounds, method="bounded", options={"xatol": 1e-6})
        return {
            "H0": float(opt.x),
            "beta": beta,
            "chi2": float(opt.fun),
            "success": bool(opt.success),
            "method": "minimize_scalar_H0",
        }

    if model.h0_fixed is not None and model.beta_fixed is None:
        h0 = model.h0_fixed

        def objective(beta: float) -> float:
            return chi2_mu(z, mu, sigma, h0, beta)

        opt = minimize_scalar(objective, bounds=beta_bounds, method="bounded", options={"xatol": 1e-6})
        return {
            "H0": h0,
            "beta": float(opt.x),
            "chi2": float(opt.fun),
            "success": bool(opt.success),
            "method": "minimize_scalar_beta",
        }

    def objective_vec(theta: np.ndarray) -> float:
        h0, beta = float(theta[0]), float(theta[1])
        return chi2_mu(z, mu, sigma, h0, beta)

    starts = [
        np.array([67.5, 0.5]),
        np.array([70.0, 0.5]),
        np.array([72.0, 0.5]),
        np.array([67.5, 0.0]),
    ]
    best = None
    for start in starts:
        opt = minimize(objective_vec, start, method="L-BFGS-B", bounds=[h0_bounds, beta_bounds])
        if best is None or float(opt.fun) < float(best.fun):
            best = opt
    assert best is not None
    return {
        "H0": float(best.x[0]),
        "beta": float(best.x[1]),
        "chi2": float(best.fun),
        "success": bool(best.success),
        "method": "L-BFGS-B",
    }


def fit_model_grid(
    z: np.ndarray,
    mu: np.ndarray,
    sigma: np.ndarray,
    model: ModelSpec,
    h0_bounds: tuple[float, float],
    beta_bounds: tuple[float, float],
) -> dict[str, float | str | bool | int]:
    h0_grid = [model.h0_fixed] if model.h0_fixed is not None else np.linspace(h0_bounds[0], h0_bounds[1], 401)
    beta_grid = [model.beta_fixed] if model.beta_fixed is not None else np.linspace(beta_bounds[0], beta_bounds[1], 401)
    best = (float("inf"), float("nan"), float("nan"))
    for h0 in h0_grid:
        for beta in beta_grid:
            val = chi2_mu(z, mu, sigma, float(h0), float(beta))
            if val < best[0]:
                best = (val, float(h0), float(beta))
    return {"H0": best[1], "beta": best[2], "chi2": best[0], "success": np.isfinite(best[0]), "method": "grid"}


def load_dataset_arrays(
    dataset: DatasetSpec,
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, list[str], str, str | None, str | None, str | None, str | None]:
    df = read_table(dataset.path)
    z_cols = auto_z_columns(df, args.z_cols)
    mu_col = find_column(df, args.mu_col, MU_CANDIDATES, required=True)
    mu_err_col = find_column(df, args.mu_err_col, MU_ERR_CANDIDATES, required=False)
    calibrator_col = find_column(df, args.calibrator_col, CALIBRATOR_CANDIDATES, required=False)
    hf_col = find_column(df, args.hf_col, HF_CANDIDATES, required=False)
    prob_col = find_column(df, args.prob_col, PROB_CANDIDATES, required=False)
    assert mu_col is not None
    return df, z_cols, mu_col, mu_err_col, calibrator_col, hf_col, prob_col


def sigma_from_muerr(muerr: np.ndarray, weight_mode: str) -> tuple[np.ndarray, str]:
    if weight_mode == "unweighted":
        return np.ones_like(muerr, dtype=float), "unweighted"
    ok = np.isfinite(muerr) & (muerr > 0.0)
    if np.all(ok):
        return muerr.astype(float), "weighted_diag_mu"
    return np.ones_like(muerr, dtype=float), "unweighted_fallback"


def analyse_dataset(dataset: DatasetSpec, args: argparse.Namespace) -> tuple[list[dict[str, object]], dict[str, object]]:
    if not dataset.path.exists():
        return [], {"dataset": dataset.label, "path": str(dataset.path), "status": "missing", "note": "File not found."}

    try:
        df, z_cols, mu_col, mu_err_col, calibrator_col, hf_col, prob_col = load_dataset_arrays(dataset, args)
    except Exception as exc:
        return [], {"dataset": dataset.label, "path": str(dataset.path), "status": "error", "note": str(exc)}

    mu = numeric_series(df, mu_col)
    muerr = numeric_series(df, mu_err_col)
    calibrator = numeric_series(df, calibrator_col, default=0.0)
    hf = numeric_series(df, hf_col)
    prob = numeric_series(df, prob_col)
    weight_modes = ["weighted", "unweighted"] if args.weight_mode == "both" else [args.weight_mode]
    rows: list[dict[str, object]] = []

    for z_col in z_cols:
        z = numeric_series(df, z_col)
        masks = subset_masks(
            z=z,
            mu=mu,
            calibrator=calibrator,
            hf=hf,
            prob=prob,
            zmin=args.zmin,
            lowz_max=args.lowz_max,
            midz_max=args.midz_max,
            prob_min=args.prob_min,
        )
        for subset_name, mask in masks.items():
            finite = mask & np.isfinite(z) & np.isfinite(mu) & (z > 0.0)
            if np.count_nonzero(finite) < args.min_n:
                continue
            zz = z[finite].astype(float)
            mm = mu[finite].astype(float)
            ee = muerr[finite].astype(float)

            for wm in weight_modes:
                sigma, weighting_label = sigma_from_muerr(ee, wm)
                if not np.all(np.isfinite(sigma)) or np.any(sigma <= 0.0):
                    sigma = np.ones_like(mm)
                    weighting_label = "unweighted_fallback"

                baseline_chi2: dict[str, float] = {}
                model_rows: list[dict[str, object]] = []
                for model in MODELS:
                    fit = fit_model(
                        zz,
                        mm,
                        sigma,
                        model,
                        h0_bounds=(args.h0_min, args.h0_max),
                        beta_bounds=(args.beta_min, args.beta_max),
                    )
                    chi2 = float(fit["chi2"])
                    k = model.k
                    n = int(len(zz))
                    dof = max(n - k, 0)
                    aic = chi2 + 2.0 * k
                    bic = chi2 + k * math.log(n) if n > 0 else float("nan")
                    row = {
                        "dataset": dataset.label,
                        "source_path": str(dataset.path),
                        "z_col": z_col,
                        "mu_col": mu_col,
                        "mu_err_col": mu_err_col or "",
                        "subset": subset_name,
                        "weighting": weighting_label,
                        "model": model.name,
                        "N": n,
                        "k": k,
                        "dof": dof,
                        "H0": float(fit["H0"]),
                        "beta": float(fit["beta"]),
                        "chi2": chi2,
                        "chi2_dof": chi2 / dof if dof > 0 else float("nan"),
                        "AIC": aic,
                        "BIC": bic,
                        "success": bool(fit["success"]),
                        "method": str(fit["method"]),
                        "H0_minus_67p5": float(fit["H0"]) - 67.5,
                        "beta_minus_0p5": float(fit["beta"]) - 0.5,
                    }
                    model_rows.append(row)
                    baseline_chi2[model.name] = chi2

                raw = next((r for r in model_rows if r["model"] == "RAW_LINEAR_H0free"), None)
                beta05 = next((r for r in model_rows if r["model"] == "FR_BETA05_H0free"), None)
                best_chi2 = min(float(r["chi2"]) for r in model_rows)
                for row in model_rows:
                    row["delta_chi2_vs_best"] = float(row["chi2"]) - best_chi2
                    row["delta_chi2_vs_raw_H0free"] = float(row["chi2"]) - float(raw["chi2"]) if raw else float("nan")
                    row["delta_chi2_vs_beta05_H0free"] = float(row["chi2"]) - float(beta05["chi2"]) if beta05 else float("nan")
                    rows.append(row)

    return rows, {
        "dataset": dataset.label,
        "path": str(dataset.path),
        "status": "ok",
        "note": f"{len(rows)} model-fit rows.",
    }


def key_results(summary: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    if summary.empty:
        return summary
    subsets = [
        "SH0ES_HF_flag",
        f"lowz_{args.zmin:g}_{args.lowz_max:g}_noncal",
        f"full_noncal_zgt{args.zmin:g}",
        "all_finite_zpos",
    ]
    models = [
        "RAW_LINEAR_H0free",
        "FR_BETA05_H0free",
        "FR_BETAfree_H0free",
        "FR_BETA05_H0fixed675",
        "FR_BETAfree_H0fixed675",
    ]
    key = summary[
        (summary["subset"].isin(subsets))
        & (summary["model"].isin(models))
        & (summary["weighting"].isin(["weighted_diag_mu", "unweighted"]))
    ].copy()
    return key.sort_values(["dataset", "z_col", "subset", "weighting", "AIC", "model"])


def write_report(path: Path, summary: pd.DataFrame, key: pd.DataFrame, readiness: pd.DataFrame) -> None:
    lines: list[str] = []
    lines.append("# Raw-MU FR Likelihood Report")
    lines.append("")
    lines.append("This standalone likelihood converts Peter Lamb's updated raw-MU correction into SN-only chi2/AIC/BIC model comparisons.")
    lines.append("")
    lines.append("Model distance:")
    lines.append("")
    lines.append("`D_obs_model = (c/H0) * z * (1 + beta z)`")
    lines.append("")
    lines.append("Peter's current proposal is the fixed case `beta = 0.5`.")
    lines.append("")

    if not readiness.empty:
        lines.append("## Dataset Readiness")
        lines.append("")
        lines.append(readiness.to_markdown(index=False))
        lines.append("")

    if not key.empty:
        show_cols = [
            "dataset",
            "z_col",
            "subset",
            "weighting",
            "model",
            "N",
            "k",
            "H0",
            "beta",
            "chi2_dof",
            "AIC",
            "BIC",
            "delta_chi2_vs_raw_H0free",
        ]
        lines.append("## Key AIC/BIC Rows")
        lines.append("")
        lines.append(key[show_cols].to_markdown(index=False, floatfmt=".5g"))
        lines.append("")

        best = key[key["weighting"] == "weighted_diag_mu"].copy()
        if not best.empty:
            idx = best.groupby(["dataset", "z_col", "subset"])["AIC"].idxmin()
            best = best.loc[idx].sort_values(["dataset", "z_col", "subset"])
            lines.append("## Best Weighted Model By Dataset/Subsample")
            lines.append("")
            lines.append(best[show_cols].to_markdown(index=False, floatfmt=".5g"))
            lines.append("")

    lines.append("## Caveats")
    lines.append("")
    lines.append("- This is a diagonal-MU or unweighted SN-only likelihood, not a full covariance analysis.")
    lines.append("- DES distance moduli carry the DES release zero-point/H0 convention, so fitted H0 should be read as an effective scale until calibration is audited.")
    lines.append("- Union3 rows are compressed distance nodes, not individual supernovae.")
    lines.append("- This script is intentionally standalone; promote to `PLamb_Test_10V6c_plus.py` only after the SN-only behaviour is understood.")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--dataset", action="append", default=None, help="Extra/replacement dataset as LABEL=PATH.")
    parser.add_argument("--replace-defaults", action="store_true")
    parser.add_argument("--z-cols", default="auto")
    parser.add_argument("--mu-col", default=None)
    parser.add_argument("--mu-err-col", default=None)
    parser.add_argument("--calibrator-col", default=None)
    parser.add_argument("--hf-col", default=None)
    parser.add_argument("--prob-col", default=None)
    parser.add_argument("--prob-min", type=float, default=None)
    parser.add_argument("--zmin", type=float, default=0.01)
    parser.add_argument("--lowz-max", type=float, default=0.15)
    parser.add_argument("--midz-max", type=float, default=0.7)
    parser.add_argument("--min-n", type=int, default=5)
    parser.add_argument("--h0-min", type=float, default=40.0)
    parser.add_argument("--h0-max", type=float, default=100.0)
    parser.add_argument("--beta-min", type=float, default=-0.25)
    parser.add_argument("--beta-max", type=float, default=1.25)
    parser.add_argument("--weight-mode", choices=["weighted", "unweighted", "both"], default="weighted")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    datasets = [] if args.replace_defaults else list(DEFAULT_DATASETS)
    if args.dataset:
        datasets.extend(parse_dataset_arg(raw) for raw in args.dataset)

    all_rows: list[dict[str, object]] = []
    readiness_rows: list[dict[str, object]] = []
    for dataset in datasets:
        rows, ready = analyse_dataset(dataset, args)
        all_rows.extend(rows)
        readiness_rows.append(ready)

    summary = pd.DataFrame(all_rows)
    readiness = pd.DataFrame(readiness_rows)
    key = key_results(summary, args)

    summary_path = args.outdir / "rawmu_fr_likelihood_summary.csv"
    key_path = args.outdir / "rawmu_fr_likelihood_key_results.csv"
    readiness_path = args.outdir / "rawmu_fr_likelihood_readiness.csv"
    report_path = args.outdir / "rawmu_fr_likelihood_report.md"

    summary.to_csv(summary_path, index=False)
    key.to_csv(key_path, index=False)
    readiness.to_csv(readiness_path, index=False)
    write_report(report_path, summary, key, readiness)

    print(f"Saved summary: {summary_path}")
    print(f"Saved key results: {key_path}")
    print(f"Saved readiness: {readiness_path}")
    print(f"Saved report: {report_path}")

    if not key.empty:
        preferred = key[
            (key["weighting"] == "weighted_diag_mu")
            & (key["model"].isin(["FR_BETA05_H0free", "FR_BETAfree_H0free", "RAW_LINEAR_H0free"]))
            & (key["subset"].isin(["all_finite_zpos", f"full_noncal_zgt{args.zmin:g}"]))
        ].copy()
        preferred = preferred.sort_values(["dataset", "z_col", "AIC"]).groupby(["dataset", "z_col"]).head(3)
        print("\nWeighted SN-only likelihood headline:")
        for _, row in preferred.iterrows():
            print(
                f"{row['dataset']} {row['z_col']} {row['subset']} {row['model']}: "
                f"H0={row['H0']:.3f} beta={row['beta']:.3f} "
                f"chi2/dof={row['chi2_dof']:.3f} AIC={row['AIC']:.3f} "
                f"dchi2_raw={row['delta_chi2_vs_raw_H0free']:+.3f}"
            )


if __name__ == "__main__":
    main()

