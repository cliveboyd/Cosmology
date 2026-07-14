r"""
Covariance-aware SN-only likelihood for the raw-MU FR/no-expansion proposal.

This is the covariance stress test for Peter Lamb's updated distance-modulus
rule. It uses the same simple model as `fit_rawmu_fr_likelihood.py`, but uses
the full covariance or inverse covariance where available:

    D_obs_model = (c / H0) * z * (1 + beta z)
    MU_model    = 5 log10(D_obs_model/Mpc) + 25
    chi2        = residual.T @ C_inv @ residual

Peter's current rule corresponds to `beta = 0.5`.

Supported covariance formats:
    * Pantheon+ `.cov`: first scalar is N, followed by flattened NxN covariance.
    * Full ASCII NxN covariance.
    * DES SN5YR `.npz`: packed upper-triangular inverse covariance.
    * Union3 compressed MU FITS matrix: first row z, first column MU, rest inverse covariance.

Outputs:
    plamb_runs/diagnostics/rawmu_fr_cov_likelihood/rawmu_fr_cov_likelihood_summary.csv
    plamb_runs/diagnostics/rawmu_fr_cov_likelihood/rawmu_fr_cov_likelihood_key_results.csv
    plamb_runs/diagnostics/rawmu_fr_cov_likelihood/rawmu_fr_cov_likelihood_report.md
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
except Exception:  # pragma: no cover
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
DEFAULT_OUTDIR = ROOT / "plamb_runs" / "diagnostics" / "rawmu_fr_cov_likelihood"


@dataclass(frozen=True)
class CovDatasetSpec:
    label: str
    data_path: Path
    cov_path: Path | None


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


DEFAULT_DATASETS = [
    CovDatasetSpec(
        "PantheonPlusSH0ES",
        DEFAULT_PANTHEON,
        ROOT / "Pantheon+SH0ES_STAT+SYS.cov",
    ),
    CovDatasetSpec(
        "DES_Dovekie_raw",
        ROOT
        / "external_datasets"
        / "current_cosmology_sources"
        / "DES_SN5YR_Dovekie"
        / "DES-Dovekie_HD.csv",
        ROOT
        / "external_datasets"
        / "current_cosmology_sources"
        / "DES_SN5YR_Dovekie"
        / "STAT+SYS.npz",
    ),
    CovDatasetSpec(
        "DES_Dovekie_PLamb",
        ROOT
        / "external_datasets"
        / "current_cosmology_sources"
        / "prepared_for_plamb"
        / "DES_SN5YR_Dovekie_PLamb.dat",
        ROOT
        / "external_datasets"
        / "current_cosmology_sources"
        / "prepared_for_plamb"
        / "DES_SN5YR_Dovekie_STAT+SYS.cov",
    ),
    CovDatasetSpec(
        "Union3_compressed",
        ROOT
        / "external_datasets"
        / "current_cosmology_sources"
        / "Union3_release"
        / "mu_mat_union3_cosmo=2_mu.fits",
        ROOT
        / "external_datasets"
        / "current_cosmology_sources"
        / "Union3_release"
        / "mu_mat_union3_cosmo=2_mu.fits",
    ),
    CovDatasetSpec(
        "Union3p1_UNITY1p8",
        ROOT
        / "external_datasets"
        / "current_cosmology_sources"
        / "Union3_release"
        / "mu_mat_union3.1_UNITY1.8_template_cosmo=2_0_mu.fits",
        ROOT
        / "external_datasets"
        / "current_cosmology_sources"
        / "Union3_release"
        / "mu_mat_union3.1_UNITY1.8_template_cosmo=2_0_mu.fits",
    ),
]


def parse_dataset_arg(raw: str) -> CovDatasetSpec:
    parts = raw.split("=")
    if len(parts) == 3:
        return CovDatasetSpec(sanitize_label(parts[0]), Path(parts[1]), Path(parts[2]))
    if len(parts) == 2:
        path = Path(parts[1])
        return CovDatasetSpec(sanitize_label(parts[0]), path, None)
    path = Path(raw)
    return CovDatasetSpec(sanitize_label(path.stem), path, None)


def load_covariance_or_precision(path: Path, n: int) -> tuple[np.ndarray | None, np.ndarray | None, str]:
    """Return `(covariance, precision, status_note)` for supported formats."""
    if path is None or not path.exists():
        return None, None, "missing covariance"

    suffix = path.suffix.lower()
    if suffix == ".npz":
        npz = np.load(path)
        n_file = int(np.asarray(npz[npz.files[0]]).ravel()[0])
        packed = np.asarray(npz[npz.files[1]], dtype=float)
        expected = n_file * (n_file + 1) // 2
        if packed.size != expected:
            raise ValueError(f"Packed DES covariance has {packed.size} entries, expected {expected}.")
        precision = np.zeros((n_file, n_file), dtype=float)
        precision[np.triu_indices(n_file)] = packed
        lower = np.tril_indices(n_file, -1)
        precision[lower] = precision.T[lower]
        if n_file != n:
            raise ValueError(f"DES precision size {n_file} does not match data N={n}.")
        return None, precision, "DES packed inverse covariance"

    if suffix in {".fits", ".fit", ".fz"}:
        from astropy.io import fits

        data = np.asarray(fits.getdata(path), dtype=float)
        if data.ndim == 2 and data.shape[0] == data.shape[1] and data.shape[0] == n + 1:
            precision = np.asarray(data[1:, 1:], dtype=float)
            return None, precision, "Union compressed inverse covariance"
        raise ValueError(f"Unsupported FITS covariance format: {path}")

    arr = np.loadtxt(path)
    if arr.ndim == 2 and arr.shape == (n, n):
        return np.asarray(arr, dtype=float), None, "full covariance matrix"

    flat = np.asarray(arr, dtype=float).ravel()
    nn = n * n
    if flat.size == nn:
        return flat.reshape(n, n), None, "flat covariance matrix"
    if flat.size == nn + 1 and int(round(flat[0])) == n:
        return flat[1:].reshape(n, n), None, "Pantheon-style covariance with leading N"
    if flat.size > nn:
        return flat[-nn:].reshape(n, n), None, "covariance using trailing NxN block"
    raise ValueError(f"Could not reshape covariance size {flat.size} to ({n},{n}).")


def precision_for_subset(
    cov: np.ndarray | None,
    precision: np.ndarray | None,
    mask: np.ndarray,
    allow_precision_submatrix: bool,
) -> tuple[np.ndarray, str]:
    idx = np.flatnonzero(mask)
    n_full = len(mask)
    if len(idx) == 0:
        raise ValueError("Empty subset.")

    full_subset = len(idx) == n_full and np.all(idx == np.arange(n_full))
    if cov is not None:
        c_sub = cov if full_subset else cov[np.ix_(idx, idx)]
        try:
            return np.linalg.inv(c_sub), "covariance subset inverted"
        except np.linalg.LinAlgError:
            return np.linalg.pinv(c_sub), "covariance subset pseudo-inverted"

    if precision is None:
        raise ValueError("No covariance or precision available.")

    if full_subset or allow_precision_submatrix:
        p_sub = precision if full_subset else precision[np.ix_(idx, idx)]
        note = "precision matrix" if full_subset else "precision submatrix approximation"
        return p_sub, note

    # For a true marginal subset, invert full precision to covariance first.
    try:
        c_full = np.linalg.inv(precision)
    except np.linalg.LinAlgError:
        c_full = np.linalg.pinv(precision)
    c_sub = c_full[np.ix_(idx, idx)]
    try:
        return np.linalg.inv(c_sub), "precision inverted to covariance, subset inverted"
    except np.linalg.LinAlgError:
        return np.linalg.pinv(c_sub), "precision inverted to covariance, subset pseudo-inverted"


def mu_model_rawfr(z: np.ndarray, h0: float, beta: float) -> np.ndarray:
    d_model = (C_KMS / h0) * z * (1.0 + beta * z)
    if np.any(d_model <= 0.0) or not np.all(np.isfinite(d_model)):
        return np.full_like(z, np.inf, dtype=float)
    return 5.0 * np.log10(d_model) + 25.0


def chi2_cov(
    z: np.ndarray,
    mu: np.ndarray,
    precision: np.ndarray,
    h0: float,
    beta: float,
    offset_mode: str,
) -> tuple[float, float]:
    model = mu_model_rawfr(z, h0, beta)
    if not np.all(np.isfinite(model)):
        return float("inf"), 0.0
    resid = mu - model
    offset = 0.0
    if offset_mode == "profile":
        ones = np.ones_like(resid)
        denom = float(ones @ precision @ ones)
        if denom > 0.0:
            offset = float((ones @ precision @ resid) / denom)
            resid = resid - offset
    chi2 = float(resid @ precision @ resid)
    return chi2, offset


def fit_model(
    z: np.ndarray,
    mu: np.ndarray,
    precision: np.ndarray,
    model: ModelSpec,
    h0_bounds: tuple[float, float],
    beta_bounds: tuple[float, float],
    offset_mode: str,
) -> dict[str, float | str | bool]:
    def objective_pair(h0: float, beta: float) -> float:
        return chi2_cov(z, mu, precision, h0, beta, offset_mode)[0]

    if not HAVE_SCIPY:
        return fit_model_grid(z, mu, precision, model, h0_bounds, beta_bounds, offset_mode)

    if model.h0_fixed is not None and model.beta_fixed is not None:
        chi2, offset = chi2_cov(z, mu, precision, model.h0_fixed, model.beta_fixed, offset_mode)
        return {
            "H0": model.h0_fixed,
            "beta": model.beta_fixed,
            "offset": offset,
            "chi2": chi2,
            "success": np.isfinite(chi2),
            "method": "fixed",
        }

    if model.h0_fixed is None and model.beta_fixed is not None:
        beta = model.beta_fixed
        opt = minimize_scalar(lambda h0: objective_pair(h0, beta), bounds=h0_bounds, method="bounded")
        chi2, offset = chi2_cov(z, mu, precision, float(opt.x), beta, offset_mode)
        return {
            "H0": float(opt.x),
            "beta": beta,
            "offset": offset,
            "chi2": chi2,
            "success": bool(opt.success),
            "method": "minimize_scalar_H0",
        }

    if model.h0_fixed is not None and model.beta_fixed is None:
        h0 = model.h0_fixed
        opt = minimize_scalar(lambda beta: objective_pair(h0, beta), bounds=beta_bounds, method="bounded")
        chi2, offset = chi2_cov(z, mu, precision, h0, float(opt.x), offset_mode)
        return {
            "H0": h0,
            "beta": float(opt.x),
            "offset": offset,
            "chi2": chi2,
            "success": bool(opt.success),
            "method": "minimize_scalar_beta",
        }

    starts = [
        np.array([67.5, 0.5]),
        np.array([68.0, 0.53]),
        np.array([71.5, 0.53]),
        np.array([67.5, 0.0]),
    ]
    best = None
    for start in starts:
        opt = minimize(
            lambda theta: objective_pair(float(theta[0]), float(theta[1])),
            start,
            method="L-BFGS-B",
            bounds=[h0_bounds, beta_bounds],
        )
        if best is None or float(opt.fun) < float(best.fun):
            best = opt
    assert best is not None
    chi2, offset = chi2_cov(z, mu, precision, float(best.x[0]), float(best.x[1]), offset_mode)
    return {
        "H0": float(best.x[0]),
        "beta": float(best.x[1]),
        "offset": offset,
        "chi2": chi2,
        "success": bool(best.success),
        "method": "L-BFGS-B",
    }


def fit_model_grid(
    z: np.ndarray,
    mu: np.ndarray,
    precision: np.ndarray,
    model: ModelSpec,
    h0_bounds: tuple[float, float],
    beta_bounds: tuple[float, float],
    offset_mode: str,
) -> dict[str, float | str | bool]:
    h0_grid = [model.h0_fixed] if model.h0_fixed is not None else np.linspace(h0_bounds[0], h0_bounds[1], 301)
    beta_grid = [model.beta_fixed] if model.beta_fixed is not None else np.linspace(beta_bounds[0], beta_bounds[1], 301)
    best = (float("inf"), float("nan"), float("nan"), 0.0)
    for h0 in h0_grid:
        for beta in beta_grid:
            chi2, offset = chi2_cov(z, mu, precision, float(h0), float(beta), offset_mode)
            if chi2 < best[0]:
                best = (chi2, float(h0), float(beta), offset)
    return {"H0": best[1], "beta": best[2], "offset": best[3], "chi2": best[0], "success": np.isfinite(best[0]), "method": "grid"}


def analyse_dataset(spec: CovDatasetSpec, args: argparse.Namespace) -> tuple[list[dict[str, object]], dict[str, object]]:
    if not spec.data_path.exists():
        return [], {"dataset": spec.label, "status": "missing", "note": f"Data file not found: {spec.data_path}"}
    if spec.cov_path is None or not spec.cov_path.exists():
        return [], {"dataset": spec.label, "status": "missing_cov", "note": f"Covariance file not found: {spec.cov_path}"}

    try:
        df = read_table(spec.data_path)
        z_cols = auto_z_columns(df, args.z_cols)
        mu_col = find_column(df, args.mu_col, MU_CANDIDATES, required=True)
        mu_err_col = find_column(df, args.mu_err_col, MU_ERR_CANDIDATES, required=False)
        calibrator_col = find_column(df, args.calibrator_col, CALIBRATOR_CANDIDATES, required=False)
        hf_col = find_column(df, args.hf_col, HF_CANDIDATES, required=False)
        prob_col = find_column(df, args.prob_col, PROB_CANDIDATES, required=False)
        assert mu_col is not None
        mu = numeric_series(df, mu_col)
        calibrator = numeric_series(df, calibrator_col, default=0.0)
        hf = numeric_series(df, hf_col)
        prob = numeric_series(df, prob_col)
        cov, precision, cov_note = load_covariance_or_precision(spec.cov_path, len(df))
    except Exception as exc:
        return [], {"dataset": spec.label, "status": "error", "note": str(exc)}

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

            try:
                p_sub, subset_cov_note = precision_for_subset(
                    cov,
                    precision,
                    finite,
                    allow_precision_submatrix=args.allow_precision_submatrix,
                )
            except Exception as exc:
                rows.append(
                    {
                        "dataset": spec.label,
                        "source_path": str(spec.data_path),
                        "cov_path": str(spec.cov_path),
                        "z_col": z_col,
                        "subset": subset_name,
                        "model": "COVARIANCE_SUBSET_ERROR",
                        "status": "error",
                        "note": str(exc),
                    }
                )
                continue

            zz = z[finite].astype(float)
            mm = mu[finite].astype(float)
            baseline = None
            model_rows: list[dict[str, object]] = []
            for model in MODELS:
                fit = fit_model(
                    zz,
                    mm,
                    p_sub,
                    model,
                    h0_bounds=(args.h0_min, args.h0_max),
                    beta_bounds=(args.beta_min, args.beta_max),
                    offset_mode=args.offset_mode,
                )
                k = model.k + int(args.offset_mode == "profile")
                n = int(len(zz))
                chi2 = float(fit["chi2"])
                dof = max(n - k, 0)
                aic = chi2 + 2.0 * k
                bic = chi2 + k * math.log(n) if n > 0 else float("nan")
                row = {
                    "dataset": spec.label,
                    "source_path": str(spec.data_path),
                    "cov_path": str(spec.cov_path),
                    "cov_note": cov_note,
                    "subset_cov_note": subset_cov_note,
                    "z_col": z_col,
                    "mu_col": mu_col,
                    "mu_err_col": mu_err_col or "",
                    "subset": subset_name,
                    "offset_mode": args.offset_mode,
                    "model": model.name,
                    "N": n,
                    "k": k,
                    "dof": dof,
                    "H0": float(fit["H0"]),
                    "beta": float(fit["beta"]),
                    "offset_mag": float(fit["offset"]),
                    "chi2": chi2,
                    "chi2_dof": chi2 / dof if dof > 0 else float("nan"),
                    "AIC": aic,
                    "BIC": bic,
                    "success": bool(fit["success"]),
                    "method": str(fit["method"]),
                    "status": "ok",
                    "note": "",
                    "H0_minus_67p5": float(fit["H0"]) - 67.5,
                    "beta_minus_0p5": float(fit["beta"]) - 0.5,
                }
                model_rows.append(row)
                if model.name == "RAW_LINEAR_H0free":
                    baseline = row
            best_chi2 = min(float(r["chi2"]) for r in model_rows)
            beta05 = next((r for r in model_rows if r["model"] == "FR_BETA05_H0free"), None)
            for row in model_rows:
                row["delta_chi2_vs_best"] = float(row["chi2"]) - best_chi2
                row["delta_chi2_vs_raw_H0free"] = float(row["chi2"]) - float(baseline["chi2"]) if baseline else float("nan")
                row["delta_chi2_vs_beta05_H0free"] = float(row["chi2"]) - float(beta05["chi2"]) if beta05 else float("nan")
                rows.append(row)

    return rows, {
        "dataset": spec.label,
        "status": "ok",
        "note": f"{len(rows)} rows; covariance={cov_note}",
        "data_path": str(spec.data_path),
        "cov_path": str(spec.cov_path),
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
    ok = summary[summary.get("status", "ok") == "ok"].copy()
    key = ok[(ok["subset"].isin(subsets)) & (ok["model"].isin(models))].copy()
    return key.sort_values(["dataset", "z_col", "subset", "AIC", "model"])


def write_report(path: Path, summary: pd.DataFrame, key: pd.DataFrame, readiness: pd.DataFrame, args: argparse.Namespace) -> None:
    lines: list[str] = []
    lines.append("# Raw-MU FR Covariance Likelihood Report")
    lines.append("")
    lines.append("This report reruns Peter Lamb's raw-MU FR/no-expansion likelihood using full covariance or precision matrices.")
    lines.append("")
    lines.append("Model distance:")
    lines.append("")
    lines.append("`D_obs_model = (c/H0) * z * (1 + beta z)`")
    lines.append("")
    lines.append("Peter's fixed correction is `beta = 0.5`.")
    lines.append("")
    lines.append(f"Offset mode: `{args.offset_mode}`.")
    lines.append("")

    if not readiness.empty:
        lines.append("## Dataset Readiness")
        lines.append("")
        lines.append(readiness.to_markdown(index=False))
        lines.append("")

    if not key.empty:
        cols = [
            "dataset",
            "z_col",
            "subset",
            "model",
            "N",
            "k",
            "H0",
            "beta",
            "offset_mag",
            "chi2_dof",
            "AIC",
            "BIC",
            "delta_chi2_vs_raw_H0free",
        ]
        lines.append("## Key AIC/BIC Rows")
        lines.append("")
        lines.append(key[cols].to_markdown(index=False, floatfmt=".5g"))
        lines.append("")

        idx = key.groupby(["dataset", "z_col", "subset"])["AIC"].idxmin()
        best = key.loc[idx].sort_values(["dataset", "z_col", "subset"])
        lines.append("## Best Model By Dataset/Subsample")
        lines.append("")
        lines.append(best[cols].to_markdown(index=False, floatfmt=".5g"))
        lines.append("")

    lines.append("## Notes")
    lines.append("")
    lines.append("- `offset_mode=fixed` is the direct test of the calibrated MU scale.")
    lines.append("- `offset_mode=profile` profiles a constant magnitude offset, similar in spirit to standard SN nuisance-magnitude marginalization; H0 is then not independently meaningful.")
    lines.append("- DES `.npz` files store the inverse covariance in packed upper-triangular form; this script follows the DES likelihood unpacking convention.")
    lines.append("- For subsets of an inverse covariance, the script inverts to covariance, subsets, and reinverts unless `--allow-precision-submatrix` is set.")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--dataset", action="append", default=None, help="Dataset as LABEL=DATA_PATH=COV_PATH.")
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
    parser.add_argument("--offset-mode", choices=["fixed", "profile"], default="fixed")
    parser.add_argument("--allow-precision-submatrix", action="store_true")
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

    suffix = f"_{args.offset_mode}"
    summary_path = args.outdir / f"rawmu_fr_cov_likelihood_summary{suffix}.csv"
    key_path = args.outdir / f"rawmu_fr_cov_likelihood_key_results{suffix}.csv"
    readiness_path = args.outdir / f"rawmu_fr_cov_likelihood_readiness{suffix}.csv"
    report_path = args.outdir / f"rawmu_fr_cov_likelihood_report{suffix}.md"

    summary.to_csv(summary_path, index=False)
    key.to_csv(key_path, index=False)
    readiness.to_csv(readiness_path, index=False)
    write_report(report_path, summary, key, readiness, args)

    print(f"Saved summary: {summary_path}")
    print(f"Saved key results: {key_path}")
    print(f"Saved readiness: {readiness_path}")
    print(f"Saved report: {report_path}")

    if not key.empty:
        view = key[
            (key["model"].isin(["FR_BETA05_H0free", "FR_BETAfree_H0free", "FR_BETA05_H0fixed675"]))
            & (key["subset"].isin(["all_finite_zpos", f"full_noncal_zgt{args.zmin:g}"]))
        ].copy()
        view = view.sort_values(["dataset", "z_col", "AIC"]).groupby(["dataset", "z_col"]).head(3)
        print("\nCovariance-aware headline:")
        for _, row in view.iterrows():
            print(
                f"{row['dataset']} {row['z_col']} {row['subset']} {row['model']}: "
                f"H0={row['H0']:.3f} beta={row['beta']:.3f} offset={row['offset_mag']:.4f} "
                f"chi2/dof={row['chi2_dof']:.3f} AIC={row['AIC']:.3f} "
                f"dchi2_raw={row['delta_chi2_vs_raw_H0free']:+.3f}"
            )


if __name__ == "__main__":
    main()

