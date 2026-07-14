r"""
Raw distance-modulus diagnostic for Peter Lamb's FR/no-expansion proposal.

This script deliberately does not alter the main PLAMB likelihood. It tests
the specific email-level prescription:

    1. Read observed supernova distance modulus MU.
    2. Convert directly to empirical luminosity distance:

           DL_obs_Mpc = 10 ** ((MU - 25) / 5)

    3. Apply no cosmological redshift correction at the conversion stage.
    4. Optionally convert the empirical distance to an FR travel-time proxy
       by dividing by a faster-light factor, especially (1 + z/2).
    5. Fit distance proxy versus redshift through the origin, and also with
       a free intercept as a diagnostic.

The output is intended to answer: does the raw Pantheon+/DES distance modulus
data support a straight-line no-expansion relation with H0 near 67.5, and how
does that depend on the chosen redshift correction?

Default input:
    C:\Users\clive\Documents\Cosmology\Pantheon+SH0ES.dat

Main outputs:
    plamb_runs/diagnostics/pantheon_rawmu_fr/rawmu_fr_fit_summary.csv
    plamb_runs/diagnostics/pantheon_rawmu_fr/rawmu_fr_beta_scan.csv
    plamb_runs/diagnostics/pantheon_rawmu_fr/rawmu_fr_point_table.csv
    plamb_runs/diagnostics/pantheon_rawmu_fr/rawmu_fr_summary.md
    plamb_runs/diagnostics/pantheon_rawmu_fr/rawmu_fr_<dataset>_<zcol>.png
"""

from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

try:
    import matplotlib.pyplot as plt

    HAVE_MATPLOTLIB = True
except Exception:  # pragma: no cover - plotting is optional.
    plt = None
    HAVE_MATPLOTLIB = False


ROOT = Path(__file__).resolve().parent
DEFAULT_OUTDIR = ROOT / "plamb_runs" / "diagnostics" / "pantheon_rawmu_fr"
DEFAULT_PANTHEON = ROOT / "Pantheon+SH0ES.dat"
C_KMS = 299792.458


Z_CANDIDATES = [
    "zCMB",
    "zHD",
    "zHEL",
    "z",
    "zcmb",
    "zhd",
    "z_cmb",
    "redshift",
    "Z",
]

MU_CANDIDATES = [
    "MU_SH0ES",
    "MU",
    "mu",
    "distance_modulus",
    "DISTMOD",
    "DM",
]

MU_ERR_CANDIDATES = [
    "MU_SH0ES_ERR_DIAG",
    "MUERR",
    "MU_ERR",
    "mu_err",
    "distance_modulus_err",
    "DISTMODERR",
    "DMERR",
]

CALIBRATOR_CANDIDATES = ["IS_CALIBRATOR", "is_calibrator", "calibrator"]
HF_CANDIDATES = ["USED_IN_SH0ES_HF", "used_in_shoes_hf", "HF", "hubble_flow"]
PROB_CANDIDATES = ["PROB", "PROBIa", "PROB_Ia", "prob", "prob_ia"]
PROB_CANDIDATES += ["PROBIA_BEAMS", "PROB1A_BEAMS", "PROBIA", "PROB1A"]


@dataclass(frozen=True)
class DatasetSpec:
    label: str
    path: Path


@dataclass(frozen=True)
class FitResult:
    dataset: str
    source_path: str
    z_col: str
    mu_col: str
    mu_err_col: str
    subset: str
    correction: str
    fit_type: str
    weighting: str
    weighted: bool
    n: int
    beta: float
    slope_mpc_per_z: float
    slope_err: float
    h0: float
    h0_err: float
    intercept_mpc: float
    intercept_err: float
    chi2: float
    dof: int
    chi2_dof: float
    rms_mpc: float
    frac_rms: float
    median_z: float
    min_z: float
    max_z: float


def sanitize_label(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
    return value.strip("_") or "dataset"


def parse_dataset_arg(raw: str) -> DatasetSpec:
    if "=" in raw:
        label, path = raw.split("=", 1)
        return DatasetSpec(sanitize_label(label), Path(path).expanduser())
    path = Path(raw).expanduser()
    return DatasetSpec(sanitize_label(path.stem), path)


def read_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        lowered = str(path).lower()
        if "path\\to" in lowered or "your_" in lowered or "example" in lowered:
            raise FileNotFoundError(
                f"Dataset not found: {path}\n"
                "That looks like an example placeholder. Replace it with a real data file path, for example:\n"
                r"  --dataset DES_Dovekie=C:\Users\clive\Documents\Cosmology\external_datasets\current_cosmology_sources\DES_SN5YR_Dovekie\DES-Dovekie_HD.csv"
            )
        raise FileNotFoundError(f"Dataset not found: {path}")

    suffix = path.suffix.lower()
    if suffix in {".fits", ".fit", ".fz"}:
        try:
            from astropy.table import Table
        except Exception as exc:  # pragma: no cover - depends on optional astropy.
            raise RuntimeError(
                f"Reading FITS requires astropy, but it is not available: {exc}"
            ) from exc
        try:
            table = Table.read(path)
            df = table.to_pandas()
        except Exception:
            df = read_union3_mu_matrix(path)
    elif suffix == ".csv":
        snana_df = read_snana_ascii(path)
        if snana_df is not None:
            df = snana_df
        else:
            df = pd.read_csv(path, comment="#")
        if len(df.columns) == 1 and str(df.columns[0]).startswith("#"):
            df = pd.read_csv(path, comment="#", sep=r"\s+", engine="python")
    else:
        df = pd.read_csv(path, sep=r"\s+", comment="#", engine="python")

    df.columns = [str(c).strip() for c in df.columns]
    return df


def read_snana_ascii(path: Path) -> pd.DataFrame | None:
    names: list[str] | None = None
    rows: list[list[str]] = []
    with path.open("r", encoding="utf-8-sig", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if line.upper().startswith("VARNAMES:"):
                names = line.split(":", 1)[1].split()
                continue
            if line.upper().startswith("SN:") and names:
                values = line.split()[1:]
                if len(values) >= len(names):
                    rows.append(values[: len(names)])

    if not names or not rows:
        return None

    df = pd.DataFrame(rows, columns=names)
    for col in df.columns:
        if col.upper() not in {"CID", "NAME", "FIELD"}:
            converted = pd.to_numeric(df[col], errors="coerce")
            if converted.notna().any():
                df[col] = converted
    return df


def read_union3_mu_matrix(path: Path) -> pd.DataFrame:
    """Read Union3 compressed MU matrix FITS files.

    Union3 documents these files as a square matrix where the first row stores
    redshift nodes, the first column stores compressed distance moduli, and the
    remaining block is the inverse covariance matrix.
    """
    from astropy.io import fits

    data = np.asarray(fits.getdata(path), dtype=float)
    if data.ndim != 2 or data.shape[0] != data.shape[1] or data.shape[0] < 3:
        raise ValueError(f"FITS file is not a Union3 compressed MU matrix: {path}")

    z = np.asarray(data[0, 1:], dtype=float)
    mu = np.asarray(data[1:, 0], dtype=float)
    invcov = np.asarray(data[1:, 1:], dtype=float)
    muerr = np.full_like(mu, np.nan, dtype=float)
    try:
        cov = np.linalg.inv(invcov)
        diag = np.diag(cov)
        muerr = np.sqrt(np.where(diag > 0.0, diag, np.nan))
    except np.linalg.LinAlgError:
        diag_inv = np.diag(invcov)
        muerr = np.sqrt(np.where(diag_inv > 0.0, 1.0 / diag_inv, np.nan))

    return pd.DataFrame(
        {
            "z": z,
            "MU": mu,
            "MUERR": muerr,
            "UNION3_NODE": np.arange(len(z), dtype=int),
        }
    )


def find_column(
    df: pd.DataFrame,
    explicit: str | None,
    candidates: Iterable[str],
    required: bool = False,
) -> str | None:
    if explicit:
        if explicit in df.columns:
            return explicit
        lower = {c.lower(): c for c in df.columns}
        if explicit.lower() in lower:
            return lower[explicit.lower()]
        if required:
            raise KeyError(f"Column {explicit!r} not found. Available columns: {list(df.columns)}")
        return None

    lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand in df.columns:
            return cand
        if cand.lower() in lower:
            return lower[cand.lower()]

    if required:
        raise KeyError(f"None of {list(candidates)} found. Available columns: {list(df.columns)}")
    return None


def auto_z_columns(df: pd.DataFrame, explicit: str | None) -> list[str]:
    if explicit and explicit.lower() != "auto":
        cols = [c.strip() for c in explicit.split(",") if c.strip()]
        missing = [c for c in cols if c not in df.columns]
        if missing:
            raise KeyError(f"Requested z column(s) missing: {missing}")
        return cols

    found: list[str] = []
    lower = {c.lower(): c for c in df.columns}
    for cand in Z_CANDIDATES:
        col = cand if cand in df.columns else lower.get(cand.lower())
        if col and col not in found:
            found.append(col)
    if not found:
        raise KeyError(f"No redshift columns found. Available columns: {list(df.columns)}")
    return found


def numeric_series(df: pd.DataFrame, col: str | None, default: float = np.nan) -> np.ndarray:
    if col is None:
        return np.full(len(df), default, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)


def distance_from_mu(mu: np.ndarray) -> np.ndarray:
    return np.power(10.0, (mu - 25.0) / 5.0)


def sigma_distance_from_mu(mu: np.ndarray, mu_err: np.ndarray) -> np.ndarray:
    d = distance_from_mu(mu)
    return (math.log(10.0) / 5.0) * d * mu_err


def correction_denominator(z: np.ndarray, name: str, beta: float = 0.5) -> np.ndarray:
    if name == "raw_DL_from_MU":
        return np.ones_like(z, dtype=float)
    if name == "DL_div_1plus_z_over_2":
        return 1.0 + 0.5 * z
    if name == "DL_div_1plus_z":
        return 1.0 + z
    if name == "DL_div_1plus_beta_z":
        return 1.0 + beta * z
    raise ValueError(f"Unknown correction: {name}")


def build_weights(sigma: np.ndarray, weighting: str) -> tuple[np.ndarray, bool]:
    if weighting == "unweighted":
        return np.ones_like(sigma, dtype=float), False
    ok = np.isfinite(sigma) & (sigma > 0.0)
    if np.count_nonzero(ok) == len(sigma):
        return 1.0 / np.square(sigma), True
    return np.ones_like(sigma, dtype=float), False


def fit_through_origin(
    z: np.ndarray,
    y: np.ndarray,
    sigma: np.ndarray,
    weighting: str,
) -> dict[str, float | int | bool]:
    w, weighted = build_weights(sigma, weighting)
    den = float(np.sum(w * z * z))
    if den <= 0.0:
        raise ValueError("Cannot fit through origin with zero redshift leverage.")

    slope = float(np.sum(w * z * y) / den)
    slope_err = float(math.sqrt(1.0 / den)) if weighted else float("nan")
    intercept = 0.0
    intercept_err = 0.0
    resid = y - slope * z
    chi2 = float(np.sum(w * resid * resid))
    dof = max(int(len(z) - 1), 0)
    return fit_stats(z, y, resid, slope, slope_err, intercept, intercept_err, chi2, dof, weighted)


def fit_with_intercept(
    z: np.ndarray,
    y: np.ndarray,
    sigma: np.ndarray,
    weighting: str,
) -> dict[str, float | int | bool]:
    w, weighted = build_weights(sigma, weighting)
    x = np.column_stack([np.ones_like(z), z])
    sw = np.sqrt(w)
    xw = x * sw[:, None]
    yw = y * sw
    coeff, *_ = np.linalg.lstsq(xw, yw, rcond=None)
    intercept = float(coeff[0])
    slope = float(coeff[1])

    if weighted:
        xtwx = x.T @ (w[:, None] * x)
        try:
            cov = np.linalg.inv(xtwx)
            intercept_err = float(math.sqrt(max(cov[0, 0], 0.0)))
            slope_err = float(math.sqrt(max(cov[1, 1], 0.0)))
        except np.linalg.LinAlgError:
            intercept_err = float("nan")
            slope_err = float("nan")
    else:
        intercept_err = float("nan")
        slope_err = float("nan")

    resid = y - (intercept + slope * z)
    chi2 = float(np.sum(w * resid * resid))
    dof = max(int(len(z) - 2), 0)
    return fit_stats(z, y, resid, slope, slope_err, intercept, intercept_err, chi2, dof, weighted)


def fit_stats(
    z: np.ndarray,
    y: np.ndarray,
    resid: np.ndarray,
    slope: float,
    slope_err: float,
    intercept: float,
    intercept_err: float,
    chi2: float,
    dof: int,
    weighted: bool,
) -> dict[str, float | int | bool]:
    h0 = C_KMS / slope if slope > 0.0 else float("nan")
    h0_err = abs(C_KMS * slope_err / (slope * slope)) if slope > 0.0 and np.isfinite(slope_err) else float("nan")
    frac = resid / np.maximum(np.abs(y), 1e-12)
    return {
        "weighted": weighted,
        "slope_mpc_per_z": slope,
        "slope_err": slope_err,
        "h0": h0,
        "h0_err": h0_err,
        "intercept_mpc": intercept,
        "intercept_err": intercept_err,
        "chi2": chi2,
        "dof": dof,
        "chi2_dof": chi2 / dof if dof > 0 else float("nan"),
        "rms_mpc": float(np.sqrt(np.mean(resid * resid))),
        "frac_rms": float(np.sqrt(np.mean(frac * frac))),
        "median_z": float(np.median(z)),
        "min_z": float(np.min(z)),
        "max_z": float(np.max(z)),
    }


def make_fit_result(
    dataset: DatasetSpec,
    z_col: str,
    mu_col: str,
    mu_err_col: str | None,
    subset: str,
    correction: str,
    fit_type: str,
    weighting: str,
    beta: float,
    z: np.ndarray,
    y: np.ndarray,
    sigma: np.ndarray,
) -> FitResult:
    if fit_type == "through_origin":
        stats = fit_through_origin(z, y, sigma, weighting)
    elif fit_type == "intercept":
        stats = fit_with_intercept(z, y, sigma, weighting)
    else:
        raise ValueError(f"Unknown fit type: {fit_type}")

    return FitResult(
        dataset=dataset.label,
        source_path=str(dataset.path),
        z_col=z_col,
        mu_col=mu_col,
        mu_err_col=mu_err_col or "",
        subset=subset,
        correction=correction,
        fit_type=fit_type,
        weighting=weighting,
        weighted=bool(stats["weighted"]),
        n=int(len(z)),
        beta=float(beta),
        slope_mpc_per_z=float(stats["slope_mpc_per_z"]),
        slope_err=float(stats["slope_err"]),
        h0=float(stats["h0"]),
        h0_err=float(stats["h0_err"]),
        intercept_mpc=float(stats["intercept_mpc"]),
        intercept_err=float(stats["intercept_err"]),
        chi2=float(stats["chi2"]),
        dof=int(stats["dof"]),
        chi2_dof=float(stats["chi2_dof"]),
        rms_mpc=float(stats["rms_mpc"]),
        frac_rms=float(stats["frac_rms"]),
        median_z=float(stats["median_z"]),
        min_z=float(stats["min_z"]),
        max_z=float(stats["max_z"]),
    )


def result_to_dict(result: FitResult) -> dict[str, object]:
    return {
        "dataset": result.dataset,
        "source_path": result.source_path,
        "z_col": result.z_col,
        "mu_col": result.mu_col,
        "mu_err_col": result.mu_err_col,
        "subset": result.subset,
        "correction": result.correction,
        "fit_type": result.fit_type,
        "weighting": result.weighting,
        "weighted": result.weighted,
        "N": result.n,
        "beta": result.beta,
        "slope_Mpc_per_z": result.slope_mpc_per_z,
        "slope_err": result.slope_err,
        "H0_km_s_Mpc": result.h0,
        "H0_err": result.h0_err,
        "intercept_Mpc": result.intercept_mpc,
        "intercept_err": result.intercept_err,
        "chi2": result.chi2,
        "dof": result.dof,
        "chi2_dof": result.chi2_dof,
        "rms_Mpc": result.rms_mpc,
        "frac_rms": result.frac_rms,
        "median_z": result.median_z,
        "min_z": result.min_z,
        "max_z": result.max_z,
    }


def subset_masks(
    z: np.ndarray,
    mu: np.ndarray,
    calibrator: np.ndarray,
    hf: np.ndarray,
    prob: np.ndarray,
    zmin: float,
    lowz_max: float,
    midz_max: float,
    prob_min: float | None,
) -> dict[str, np.ndarray]:
    finite = np.isfinite(z) & np.isfinite(mu) & (z > 0.0)
    if prob_min is not None:
        finite &= np.isfinite(prob) & (prob >= prob_min)

    is_cal = np.isfinite(calibrator) & (calibrator != 0.0)
    noncal = ~is_cal
    masks = {
        "all_finite_zpos": finite,
        f"noncal_zgt{zmin:g}": finite & noncal & (z > zmin),
        f"lowz_{zmin:g}_{lowz_max:g}_noncal": finite & noncal & (z > zmin) & (z < lowz_max),
        f"midz_{lowz_max:g}_{midz_max:g}_noncal": finite & noncal & (z >= lowz_max) & (z < midz_max),
        f"highz_ge{midz_max:g}_noncal": finite & noncal & (z >= midz_max),
        f"full_noncal_zgt{zmin:g}": finite & noncal & (z > zmin),
    }
    if np.any(np.isfinite(hf)):
        masks["SH0ES_HF_flag"] = finite & (hf != 0.0)
    return masks


def analyse_dataset(
    dataset: DatasetSpec,
    args: argparse.Namespace,
    outdir: Path,
) -> tuple[list[FitResult], list[dict[str, object]], list[dict[str, object]]]:
    df = read_table(dataset.path)
    z_cols = auto_z_columns(df, args.z_cols)
    mu_col = find_column(df, args.mu_col, MU_CANDIDATES, required=True)
    mu_err_col = find_column(df, args.mu_err_col, MU_ERR_CANDIDATES, required=False)
    calibrator_col = find_column(df, args.calibrator_col, CALIBRATOR_CANDIDATES, required=False)
    hf_col = find_column(df, args.hf_col, HF_CANDIDATES, required=False)
    prob_col = find_column(df, args.prob_col, PROB_CANDIDATES, required=False)

    assert mu_col is not None
    mu = numeric_series(df, mu_col)
    mu_err = numeric_series(df, mu_err_col)
    calibrator = numeric_series(df, calibrator_col, default=0.0)
    hf = numeric_series(df, hf_col)
    prob = numeric_series(df, prob_col)
    dl_obs = distance_from_mu(mu)
    sigma_dl = sigma_distance_from_mu(mu, mu_err)

    fit_results: list[FitResult] = []
    beta_rows: list[dict[str, object]] = []
    point_rows: list[dict[str, object]] = []
    if args.weight_mode == "both":
        weightings = ["weighted_diag_mu", "unweighted"]
    elif args.weight_mode == "weighted":
        weightings = ["weighted_diag_mu"]
    else:
        weightings = ["unweighted"]

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

        finite = np.isfinite(z) & np.isfinite(mu) & (z > 0.0)
        denom_peter = correction_denominator(z, "DL_div_1plus_z_over_2")
        denom_oneplus = correction_denominator(z, "DL_div_1plus_z")
        for idx in np.flatnonzero(finite):
            point_rows.append(
                {
                    "dataset": dataset.label,
                    "source_path": str(dataset.path),
                    "row_index": int(idx),
                    "z_col": z_col,
                    "z": float(z[idx]),
                    "mu_col": mu_col,
                    "MU": float(mu[idx]),
                    "mu_err_col": mu_err_col or "",
                    "MU_err": float(mu_err[idx]) if np.isfinite(mu_err[idx]) else float("nan"),
                    "DL_obs_Mpc": float(dl_obs[idx]),
                    "sigma_DL_Mpc": float(sigma_dl[idx]) if np.isfinite(sigma_dl[idx]) else float("nan"),
                    "DL_div_1plus_z_over_2_Mpc": float(dl_obs[idx] / denom_peter[idx]),
                    "DL_div_1plus_z_Mpc": float(dl_obs[idx] / denom_oneplus[idx]),
                    "is_calibrator": float(calibrator[idx]) if np.isfinite(calibrator[idx]) else float("nan"),
                    "used_in_shoes_hf": float(hf[idx]) if np.isfinite(hf[idx]) else float("nan"),
                    "prob": float(prob[idx]) if np.isfinite(prob[idx]) else float("nan"),
                }
            )

        for subset_name, mask in masks.items():
            mask = mask & np.isfinite(dl_obs)
            if np.count_nonzero(mask) < args.min_n:
                continue

            zz = z[mask]
            dd = dl_obs[mask]
            sd = sigma_dl[mask]

            for weighting in weightings:
                for correction, beta in [
                    ("raw_DL_from_MU", 0.0),
                    ("DL_div_1plus_z_over_2", 0.5),
                    ("DL_div_1plus_z", 1.0),
                ]:
                    denom = correction_denominator(zz, correction, beta)
                    ok = np.isfinite(denom) & (np.abs(denom) > 1e-12)
                    if np.count_nonzero(ok) < args.min_n:
                        continue
                    y = dd[ok] / denom[ok]
                    sig = sd[ok] / np.abs(denom[ok])
                    for fit_type in ["through_origin", "intercept"]:
                        fit_results.append(
                            make_fit_result(
                                dataset=dataset,
                                z_col=z_col,
                                mu_col=mu_col,
                                mu_err_col=mu_err_col,
                                subset=subset_name,
                                correction=correction,
                                fit_type=fit_type,
                                weighting=weighting,
                                beta=beta,
                                z=zz[ok],
                                y=y,
                                sigma=sig,
                            )
                        )

                best_result: FitResult | None = None
                best_score = float("inf")
                for beta in np.linspace(args.beta_min, args.beta_max, args.beta_steps):
                    denom = correction_denominator(zz, "DL_div_1plus_beta_z", float(beta))
                    ok = np.isfinite(denom) & (denom > 1e-6)
                    if np.count_nonzero(ok) < args.min_n:
                        continue
                    y = dd[ok] / denom[ok]
                    sig = sd[ok] / np.abs(denom[ok])
                    result = make_fit_result(
                        dataset=dataset,
                        z_col=z_col,
                        mu_col=mu_col,
                        mu_err_col=mu_err_col,
                        subset=subset_name,
                        correction="DL_div_1plus_beta_z",
                        fit_type="through_origin",
                        weighting=weighting,
                        beta=float(beta),
                        z=zz[ok],
                        y=y,
                        sigma=sig,
                    )
                    row = result_to_dict(result)
                    beta_rows.append(row)
                    score = result.chi2_dof if np.isfinite(result.chi2_dof) else result.frac_rms
                    if score < best_score:
                        best_score = score
                        best_result = result

                if best_result is not None:
                    fit_results.append(
                        FitResult(
                            dataset=best_result.dataset,
                            source_path=best_result.source_path,
                            z_col=best_result.z_col,
                            mu_col=best_result.mu_col,
                            mu_err_col=best_result.mu_err_col,
                            subset=best_result.subset,
                            correction="DL_div_1plus_beta_z_best",
                            fit_type=best_result.fit_type,
                            weighting=best_result.weighting,
                            weighted=best_result.weighted,
                            n=best_result.n,
                            beta=best_result.beta,
                            slope_mpc_per_z=best_result.slope_mpc_per_z,
                            slope_err=best_result.slope_err,
                            h0=best_result.h0,
                            h0_err=best_result.h0_err,
                            intercept_mpc=best_result.intercept_mpc,
                            intercept_err=best_result.intercept_err,
                            chi2=best_result.chi2,
                            dof=best_result.dof,
                            chi2_dof=best_result.chi2_dof,
                            rms_mpc=best_result.rms_mpc,
                            frac_rms=best_result.frac_rms,
                            median_z=best_result.median_z,
                            min_z=best_result.min_z,
                            max_z=best_result.max_z,
                        )
                    )

    return fit_results, beta_rows, point_rows


def write_summary(
    outpath: Path,
    fit_df: pd.DataFrame,
    beta_df: pd.DataFrame,
    point_df: pd.DataFrame,
    args: argparse.Namespace,
) -> None:
    lines: list[str] = []
    lines.append("# Pantheon raw-MU FR diagnostic")
    lines.append("")
    lines.append("This diagnostic tests Peter Lamb's latest no-expansion distance-modulus prescription without changing the main PLAMB/SU2 likelihood code.")
    lines.append("")
    lines.append("Method:")
    lines.append("")
    lines.append("- Convert observed supernova distance modulus directly to empirical distance: `DL_obs_Mpc = 10 ** ((MU - 25) / 5)`.")
    lines.append("- Fit distance proxy versus redshift through the origin and with a free intercept.")
    lines.append("- Compare no correction, Peter's `/ (1 + z/2)` travel-time proxy, `/ (1 + z)`, and a fitted `/ (1 + beta z)` scan.")
    lines.append("- Report diagonal-MU weighted fits and/or unweighted fits according to `--weight-mode`. This is not a full covariance likelihood.")
    lines.append("")
    lines.append("Important interpretation caveats:")
    lines.append("")
    lines.append("- Pantheon+ already contains `MU_SH0ES`; the issue is not missing MU, but the correction convention and subset choice.")
    lines.append("- This raw diagnostic does not refit light-curve standardization, selection corrections, or the full Pantheon+ covariance.")
    lines.append("- A stable raw-MU signal should be reproduced in DES or other SN compilations before modifying the main PLAMB branch.")
    lines.append("")

    headline_subsets = [
        "SH0ES_HF_flag",
        f"lowz_{args.zmin:g}_{args.lowz_max:g}_noncal",
        f"full_noncal_zgt{args.zmin:g}",
    ]
    headline_corrections = [
        "raw_DL_from_MU",
        "DL_div_1plus_z_over_2",
        "DL_div_1plus_z",
        "DL_div_1plus_beta_z_best",
    ]
    headline = fit_df[
        (fit_df["fit_type"] == "through_origin")
        & (fit_df["subset"].isin(headline_subsets))
        & (fit_df["correction"].isin(headline_corrections))
    ].copy()
    if not headline.empty:
        headline = headline.sort_values(["dataset", "z_col", "subset", "correction"])
        cols = ["dataset", "z_col", "subset", "weighting", "correction", "N", "beta", "H0_km_s_Mpc", "chi2_dof", "frac_rms"]
        lines.append("## Headline through-origin fits")
        lines.append("")
        lines.append(headline[cols].to_markdown(index=False, floatfmt=".5g"))
        lines.append("")

    if not beta_df.empty:
        best_beta = fit_df[
            (fit_df["fit_type"] == "through_origin")
            & (fit_df["correction"] == "DL_div_1plus_beta_z_best")
        ].copy()
        if not best_beta.empty:
            best_beta = best_beta.sort_values(["dataset", "z_col", "subset"])
            cols = ["dataset", "z_col", "subset", "weighting", "N", "beta", "H0_km_s_Mpc", "chi2_dof", "frac_rms"]
            lines.append("## Best fitted beta values")
            lines.append("")
            lines.append("The fitted rule is `DL_obs / (1 + beta z)`. Peter's email corresponds approximately to `beta = 0.5`.")
            lines.append("")
            lines.append(best_beta[cols].to_markdown(index=False, floatfmt=".5g"))
            lines.append("")

    lines.append("## Output files")
    lines.append("")
    lines.append(f"- Fit summary rows: `{len(fit_df)}`")
    lines.append(f"- Beta scan rows: `{len(beta_df)}`")
    lines.append(f"- Point table rows: `{len(point_df)}`")
    lines.append("")
    lines.append("Recommended next check: run the same diagnostic on the DES SN distance-modulus table, then compare the preferred correction and H0 by subset.")
    lines.append("")

    outpath.write_text("\n".join(lines), encoding="utf-8")


def make_plots(outdir: Path, fit_df: pd.DataFrame, point_df: pd.DataFrame, args: argparse.Namespace) -> list[Path]:
    if not HAVE_MATPLOTLIB or args.no_plots:
        return []

    paths: list[Path] = []
    for (dataset, z_col), pts in point_df.groupby(["dataset", "z_col"]):
        pts = pts.replace([np.inf, -np.inf], np.nan).dropna(subset=["z", "DL_obs_Mpc", "DL_div_1plus_z_over_2_Mpc"])
        if pts.empty:
            continue
        pts = pts.sort_values("z")
        max_points = min(len(pts), args.plot_max_points)
        plot_pts = pts.sample(max_points, random_state=42) if len(pts) > max_points else pts

        fig, ax = plt.subplots(figsize=(9.5, 6.0))
        ax.scatter(plot_pts["z"], plot_pts["DL_obs_Mpc"], s=10, alpha=0.24, label="raw DL from MU")
        ax.scatter(
            plot_pts["z"],
            plot_pts["DL_div_1plus_z_over_2_Mpc"],
            s=10,
            alpha=0.24,
            label="DL / (1 + z/2)",
        )

        zmax = float(np.nanmax(pts["z"]))
        xline = np.linspace(0.0, zmax, 300)
        subset_order = [
            "SH0ES_HF_flag",
            f"lowz_{args.zmin:g}_{args.lowz_max:g}_noncal",
            f"full_noncal_zgt{args.zmin:g}",
        ]
        chosen_subset = None
        for subset in subset_order:
            if np.any(
                (fit_df["dataset"] == dataset)
                & (fit_df["z_col"] == z_col)
                & (fit_df["subset"] == subset)
                & (fit_df["fit_type"] == "through_origin")
            ):
                chosen_subset = subset
                break
        if chosen_subset:
            preferred_weighting = "weighted_diag_mu"
            if not np.any(
                (fit_df["dataset"] == dataset)
                & (fit_df["z_col"] == z_col)
                & (fit_df["subset"] == chosen_subset)
                & (fit_df["fit_type"] == "through_origin")
                & (fit_df["weighting"] == preferred_weighting)
            ):
                preferred_weighting = "unweighted"
            for correction, style in [
                ("raw_DL_from_MU", "-"),
                ("DL_div_1plus_z_over_2", "--"),
                ("DL_div_1plus_beta_z_best", ":"),
            ]:
                row = fit_df[
                    (fit_df["dataset"] == dataset)
                    & (fit_df["z_col"] == z_col)
                    & (fit_df["subset"] == chosen_subset)
                    & (fit_df["fit_type"] == "through_origin")
                    & (fit_df["weighting"] == preferred_weighting)
                    & (fit_df["correction"] == correction)
                ]
                if row.empty:
                    continue
                row0 = row.iloc[0]
                label = f"{correction}, H0={row0['H0_km_s_Mpc']:.2f}"
                if correction.endswith("best"):
                    label += f", beta={row0['beta']:.3f}"
                ax.plot(xline, row0["slope_Mpc_per_z"] * xline, style, linewidth=2.0, label=label)

        ax.set_xlabel(z_col)
        ax.set_ylabel("Distance proxy (Mpc)")
        ax.set_title(f"{dataset}: raw MU FR diagnostic ({z_col})")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
        fig.tight_layout()
        path = outdir / f"rawmu_fr_{sanitize_label(dataset)}_{sanitize_label(z_col)}.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        paths.append(path)
    return paths


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset",
        action="append",
        default=None,
        help=(
            "Dataset as LABEL=PATH. May be repeated. "
            "If omitted, uses PantheonPlusSH0ES=Pantheon+SH0ES.dat."
        ),
    )
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--z-cols", default="auto", help="Comma-separated z columns, or auto.")
    parser.add_argument("--mu-col", default=None, help="Override distance-modulus column.")
    parser.add_argument("--mu-err-col", default=None, help="Override distance-modulus error column.")
    parser.add_argument("--calibrator-col", default=None)
    parser.add_argument("--hf-col", default=None)
    parser.add_argument("--prob-col", default=None)
    parser.add_argument("--prob-min", type=float, default=None, help="Optional probability cut, e.g. 0.67 for DES Ia probability.")
    parser.add_argument("--zmin", type=float, default=0.01)
    parser.add_argument("--lowz-max", type=float, default=0.15)
    parser.add_argument("--midz-max", type=float, default=0.7)
    parser.add_argument("--min-n", type=int, default=5)
    parser.add_argument("--beta-min", type=float, default=-0.25)
    parser.add_argument("--beta-max", type=float, default=1.25)
    parser.add_argument("--beta-steps", type=int, default=301)
    parser.add_argument(
        "--weight-mode",
        choices=["weighted", "unweighted", "both"],
        default="both",
        help="Use diagonal MU errors, ignore them, or report both.",
    )
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--plot-max-points", type=int, default=5000)
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    if args.dataset:
        datasets = [parse_dataset_arg(raw) for raw in args.dataset]
    else:
        datasets = [DatasetSpec("PantheonPlusSH0ES", DEFAULT_PANTHEON)]

    all_results: list[FitResult] = []
    all_beta_rows: list[dict[str, object]] = []
    all_point_rows: list[dict[str, object]] = []
    for dataset in datasets:
        results, beta_rows, point_rows = analyse_dataset(dataset, args, outdir)
        all_results.extend(results)
        all_beta_rows.extend(beta_rows)
        all_point_rows.extend(point_rows)

    fit_df = pd.DataFrame([result_to_dict(result) for result in all_results])
    beta_df = pd.DataFrame(all_beta_rows)
    point_df = pd.DataFrame(all_point_rows)

    fit_path = outdir / "rawmu_fr_fit_summary.csv"
    beta_path = outdir / "rawmu_fr_beta_scan.csv"
    point_path = outdir / "rawmu_fr_point_table.csv"
    summary_path = outdir / "rawmu_fr_summary.md"

    fit_df.to_csv(fit_path, index=False)
    beta_df.to_csv(beta_path, index=False)
    point_df.to_csv(point_path, index=False)
    write_summary(summary_path, fit_df, beta_df, point_df, args)
    plot_paths = make_plots(outdir, fit_df, point_df, args)

    print(f"Saved fit summary: {fit_path}")
    print(f"Saved beta scan: {beta_path}")
    print(f"Saved point table: {point_path}")
    print(f"Saved summary: {summary_path}")
    for path in plot_paths:
        print(f"Saved plot: {path}")

    headline = fit_df[
        (fit_df["fit_type"] == "through_origin")
        & (fit_df["correction"].isin(["raw_DL_from_MU", "DL_div_1plus_z_over_2", "DL_div_1plus_beta_z_best"]))
    ].copy()
    priority = [
        "SH0ES_HF_flag",
        f"lowz_{args.zmin:g}_{args.lowz_max:g}_noncal",
        f"full_noncal_zgt{args.zmin:g}",
    ]
    shown = 0
    for subset in priority:
        rows = headline[headline["subset"] == subset]
        if rows.empty:
            continue
        print(f"\n=== {subset} ===")
        rows = rows.sort_values(["dataset", "z_col", "weighting", "correction"])
        for _, row in rows.iterrows():
            print(
                f"{row['dataset']} {row['z_col']} {row['weighting']} {row['correction']}: "
                f"N={int(row['N'])} beta={row['beta']:.3f} "
                f"H0={row['H0_km_s_Mpc']:.3f} chi2/dof={row['chi2_dof']:.3f} "
                f"frac_rms={row['frac_rms']:.4f}"
            )
            shown += 1
        if shown >= 24:
            break


if __name__ == "__main__":
    main()
