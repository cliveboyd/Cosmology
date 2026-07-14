r"""Dedicated local-flow/redshift-correction audit for the PLAMB log kernel.

This diagnostic follows the exact logarithmic PLAMB branch versus flat LCDM
question into the redshift machinery itself:

* compare zHD, zCMB, and zHEL where available;
* repeat full-covariance fits under low-z cuts;
* quantify object-level redshift corrections such as zHD-zHEL;
* test whether keeping only the lowest correction-magnitude objects changes
  the PLAMB/LCDM ranking;
* summarize survey blocks, velocity-scale corrections, and metadata links.

Sign conventions:

    delta_AIC_log_minus_lcdm < 0  => exact log kernel preferred
    delta_AIC_log_minus_lcdm > 0  => LCDM preferred

    diag_chi2_delta_lcdm_minus_log > 0 => exact log kernel improves that row

The diagonal row/block summaries are localization diagnostics only.  Formal
model ranking comes from the full covariance refits.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from scipy.stats import spearmanr

    HAVE_SCIPY = True
except Exception:  # pragma: no cover
    spearmanr = None
    HAVE_SCIPY = False

from diagnose_pantheon_rawmu_fr import C_KMS, read_table
from diagnose_plamb_log_kernel_residuals import (
    DEFAULT_DES,
    DEFAULT_DES_COV,
    DEFAULT_PANTHEON,
    DatasetConfig,
    PreparedDataset,
    cid_values,
    comparison_from_residuals,
    dataset_configs,
    fit_two_models,
    numeric,
    survey_values,
)
from fit_rawmu_fr_cov_likelihood import load_covariance_or_precision


ROOT = Path(__file__).resolve().parent
DEFAULT_OUTDIR = ROOT / "plamb_runs" / "diagnostics" / "plamb_local_flow_redshift_audit"
DEFAULT_PANTHEON_COV = ROOT / "Pantheon+SH0ES_STAT+SYS.cov"
Z_PRIORITY = ["zHD", "zCMB", "zHEL"]
FLOW_FEATURES = [
    "dz_hd_hel",
    "abs_dz_hd_hel",
    "v_hd_hel_kms",
    "abs_v_hd_hel_kms",
    "dz_hd_cmb",
    "abs_dz_hd_cmb",
    "v_hd_cmb_kms",
    "abs_v_hd_cmb_kms",
    "dz_cmb_hel",
    "abs_dz_cmb_hel",
    "v_cmb_hel_kms",
    "abs_v_cmb_hel_kms",
    "zHDERR",
    "zCMBERR",
    "zHELERR",
    "VPEC",
    "VPECERR",
    "MUERR_VPEC",
    "MUERR_SYS",
]

_COV_CACHE: dict[tuple[str, int], tuple[np.ndarray | None, np.ndarray | None, str]] = {}
_PREC_TO_COV_CACHE: dict[tuple[str, int], np.ndarray] = {}


def parse_float_list(raw: str) -> list[float]:
    out: list[float] = []
    raw_text = str(raw).strip()
    if raw_text.lower() in {"", "none", "off", "false", "no"}:
        return out
    for item in raw_text.split(","):
        item = item.strip()
        if item:
            out.append(float(item))
    return out


def make_args(base: argparse.Namespace, **updates: object) -> SimpleNamespace:
    values = vars(base).copy()
    values.update(updates)
    return SimpleNamespace(**values)


def redshift_blend(z_start: np.ndarray, z_end: np.ndarray, alpha: float, mode: str) -> np.ndarray:
    if mode == "multiplicative":
        start = 1.0 + np.asarray(z_start, dtype=float)
        end = 1.0 + np.asarray(z_end, dtype=float)
        ok = np.isfinite(start) & np.isfinite(end) & (start > 0.0) & (end > 0.0)
        out = np.full_like(start, np.nan, dtype=float)
        out[ok] = start[ok] * np.power(end[ok] / start[ok], alpha) - 1.0
        return out
    return np.asarray(z_start, dtype=float) + alpha * (np.asarray(z_end, dtype=float) - np.asarray(z_start, dtype=float))


def requested_z_columns(raw_df: pd.DataFrame, requested: str) -> list[str]:
    if requested.lower() == "all":
        return [col for col in Z_PRIORITY if col in raw_df.columns]
    return [col.strip() for col in requested.split(",") if col.strip() and col.strip() in raw_df.columns]


def base_mask(raw: pd.DataFrame, cfg: DatasetConfig, args: argparse.Namespace) -> np.ndarray:
    mu = numeric(raw, cfg.mu_col)
    sigma = numeric(raw, cfg.err_col)
    mask = np.isfinite(mu) & np.isfinite(sigma) & (sigma > 0.0)
    if cfg.exclude_calibrators and cfg.calibrator_col and cfg.calibrator_col in raw.columns and not args.keep_calibrators:
        mask &= numeric(raw, cfg.calibrator_col, default=0.0) == 0.0
    return mask


def cached_covariance(path: Path, n: int) -> tuple[np.ndarray | None, np.ndarray | None, str]:
    key = (str(path), int(n))
    if key not in _COV_CACHE:
        _COV_CACHE[key] = load_covariance_or_precision(path, n)
    return _COV_CACHE[key]


def invert_matrix(mat: np.ndarray) -> np.ndarray:
    try:
        return np.linalg.inv(mat)
    except np.linalg.LinAlgError:
        return np.linalg.pinv(mat)


def precision_for_mask(
    cov: np.ndarray | None,
    precision: np.ndarray | None,
    cov_path: Path,
    n_full: int,
    mask: np.ndarray,
    allow_precision_submatrix: bool,
) -> tuple[np.ndarray, str]:
    idx = np.flatnonzero(mask)
    if len(idx) == 0:
        raise ValueError("Empty subset.")
    full_subset = len(idx) == n_full and np.all(idx == np.arange(n_full))

    if cov is not None:
        c_sub = cov if full_subset else cov[np.ix_(idx, idx)]
        return invert_matrix(c_sub), "covariance subset inverted"

    if precision is None:
        raise ValueError("No covariance or precision available.")

    if full_subset:
        return precision, "precision matrix"

    if allow_precision_submatrix:
        return precision[np.ix_(idx, idx)], "precision submatrix approximation"

    key = (str(cov_path), int(n_full))
    if key not in _PREC_TO_COV_CACHE:
        _PREC_TO_COV_CACHE[key] = invert_matrix(precision)
    c_sub = _PREC_TO_COV_CACHE[key][np.ix_(idx, idx)]
    return invert_matrix(c_sub), "precision inverted to covariance, subset inverted"


def add_redshift_corrections(raw: pd.DataFrame, cfg: DatasetConfig, args: argparse.Namespace) -> pd.DataFrame:
    mask = base_mask(raw, cfg, args)
    out = pd.DataFrame(
        {
            "dataset": cfg.label,
            "source_index": np.arange(len(raw), dtype=int),
            "CID": raw["CID"].astype(str) if "CID" in raw.columns else np.arange(len(raw)).astype(str),
            "IDSURVEY": raw["IDSURVEY"].astype(str) if "IDSURVEY" in raw.columns else cfg.label,
            "base_mask": mask,
        }
    )
    for col in ["zHD", "zCMB", "zHEL", cfg.mu_col, cfg.err_col, "zHDERR", "zCMBERR", "zHELERR", "VPEC", "VPECERR", "MUERR_VPEC", "MUERR_SYS"]:
        if col in raw.columns and col not in out.columns:
            out[col] = numeric(raw, col)

    zhd = out["zHD"].to_numpy(dtype=float) if "zHD" in out.columns else np.full(len(out), np.nan)
    zcmb = out["zCMB"].to_numpy(dtype=float) if "zCMB" in out.columns else np.full(len(out), np.nan)
    zhel = out["zHEL"].to_numpy(dtype=float) if "zHEL" in out.columns else np.full(len(out), np.nan)
    denom_hel = np.where(np.isfinite(zhel), 1.0 + zhel, np.where(np.isfinite(zhd), 1.0 + zhd, np.nan))
    denom_cmb = np.where(np.isfinite(zcmb), 1.0 + zcmb, np.where(np.isfinite(zhd), 1.0 + zhd, np.nan))

    out["dz_hd_hel"] = zhd - zhel
    out["abs_dz_hd_hel"] = np.abs(out["dz_hd_hel"])
    out["v_hd_hel_kms"] = C_KMS * out["dz_hd_hel"] / denom_hel
    out["abs_v_hd_hel_kms"] = np.abs(out["v_hd_hel_kms"])

    out["dz_hd_cmb"] = zhd - zcmb
    out["abs_dz_hd_cmb"] = np.abs(out["dz_hd_cmb"])
    out["v_hd_cmb_kms"] = C_KMS * out["dz_hd_cmb"] / denom_cmb
    out["abs_v_hd_cmb_kms"] = np.abs(out["v_hd_cmb_kms"])

    out["dz_cmb_hel"] = zcmb - zhel
    out["abs_dz_cmb_hel"] = np.abs(out["dz_cmb_hel"])
    out["v_cmb_hel_kms"] = C_KMS * out["dz_cmb_hel"] / denom_hel
    out["abs_v_cmb_hel_kms"] = np.abs(out["v_cmb_hel_kms"])
    return out


def prepare_masked_dataset(
    cfg: DatasetConfig,
    z_col: str,
    args: argparse.Namespace,
    extra_mask: np.ndarray | None = None,
    subset_note_extra: str = "",
) -> PreparedDataset:
    raw = read_table(cfg.data_path)
    if z_col not in raw.columns:
        raise KeyError(f"{cfg.label} has no z column '{z_col}'.")
    if cfg.mu_col not in raw.columns:
        raise KeyError(f"{cfg.label} has no MU column '{cfg.mu_col}'.")

    z_all = numeric(raw, z_col)
    mu_all = numeric(raw, cfg.mu_col)
    sigma_all = numeric(raw, cfg.err_col)
    mask = base_mask(raw, cfg, args)
    mask &= np.isfinite(z_all) & np.isfinite(mu_all) & (z_all > args.z_min)
    if args.max_z is not None:
        mask &= z_all <= args.max_z
    if extra_mask is not None:
        mask &= extra_mask
    if int(np.count_nonzero(mask)) < args.min_n:
        raise ValueError(f"{cfg.label} {z_col} selected fewer than {args.min_n} rows.")

    cov, precision, cov_note = cached_covariance(cfg.cov_path, len(raw))
    precision_sub, subset_note = precision_for_mask(
        cov,
        precision,
        cfg.cov_path,
        len(raw),
        mask,
        allow_precision_submatrix=args.allow_precision_submatrix,
    )
    if subset_note_extra:
        subset_note = f"{subset_note}; {subset_note_extra}"

    subset = raw.loc[mask].copy().reset_index(drop=True)
    source_indices = np.flatnonzero(mask)
    sigma = sigma_all[mask].astype(float)
    sigma[~np.isfinite(sigma) | (sigma <= 0.0)] = np.nan
    return PreparedDataset(
        label=cfg.label,
        z_col=z_col,
        mu_col=cfg.mu_col,
        err_col=cfg.err_col,
        df=subset,
        source_indices=source_indices,
        z=z_all[mask].astype(float),
        mu=mu_all[mask].astype(float),
        sigma_diag=sigma,
        precision=precision_sub,
        cov_note=cov_note,
        subset_note=subset_note,
    )


def prepare_custom_z_dataset(
    cfg: DatasetConfig,
    z_label: str,
    z_values: np.ndarray,
    args: argparse.Namespace,
    extra_mask: np.ndarray | None = None,
    subset_note_extra: str = "",
) -> PreparedDataset:
    raw = read_table(cfg.data_path)
    if cfg.mu_col not in raw.columns:
        raise KeyError(f"{cfg.label} has no MU column '{cfg.mu_col}'.")

    z_all = np.asarray(z_values, dtype=float)
    if len(z_all) != len(raw):
        raise ValueError(f"Custom z array length {len(z_all)} does not match {cfg.label} N={len(raw)}.")
    mu_all = numeric(raw, cfg.mu_col)
    sigma_all = numeric(raw, cfg.err_col)
    mask = base_mask(raw, cfg, args)
    mask &= np.isfinite(z_all) & np.isfinite(mu_all) & (z_all > args.z_min)
    if args.max_z is not None:
        mask &= z_all <= args.max_z
    if extra_mask is not None:
        mask &= extra_mask
    if int(np.count_nonzero(mask)) < args.min_n:
        raise ValueError(f"{cfg.label} {z_label} selected fewer than {args.min_n} rows.")

    cov, precision, cov_note = cached_covariance(cfg.cov_path, len(raw))
    precision_sub, subset_note = precision_for_mask(
        cov,
        precision,
        cfg.cov_path,
        len(raw),
        mask,
        allow_precision_submatrix=args.allow_precision_submatrix,
    )
    if subset_note_extra:
        subset_note = f"{subset_note}; {subset_note_extra}"

    subset = raw.loc[mask].copy().reset_index(drop=True)
    source_indices = np.flatnonzero(mask)
    sigma = sigma_all[mask].astype(float)
    sigma[~np.isfinite(sigma) | (sigma <= 0.0)] = np.nan
    return PreparedDataset(
        label=cfg.label,
        z_col=z_label,
        mu_col=cfg.mu_col,
        err_col=cfg.err_col,
        df=subset,
        source_indices=source_indices,
        z=z_all[mask].astype(float),
        mu=mu_all[mask].astype(float),
        sigma_diag=sigma,
        precision=precision_sub,
        cov_note=cov_note,
        subset_note=subset_note,
    )


def summarize_fit_pair(fits: pd.DataFrame, sensitivity: str, subset_label: str, extra: dict[str, object]) -> dict[str, object]:
    log = fits[fits["model"] == "EXACT_LOG_Pfree"].iloc[0]
    lcdm = fits[fits["model"] == "LCDM_Omfree"].iloc[0]
    row: dict[str, object] = {
        "sensitivity": sensitivity,
        "subset": subset_label,
        "dataset": log["dataset"],
        "z_col": log["z_col"],
        "N": int(log["N"]),
        "H0": float(log["H0"]),
        "log_p": float(log["p"]),
        "log_beta_equiv": float(log["beta_equiv"]),
        "log_gamma_equiv": float(log["gamma_equiv"]),
        "lcdm_Om": float(lcdm["Om"]),
        "chi2_log": float(log["chi2"]),
        "chi2_lcdm": float(lcdm["chi2"]),
        "delta_chi2_log_minus_lcdm": float(log["chi2"] - lcdm["chi2"]),
        "chi2_dof_log": float(log["chi2_dof"]),
        "chi2_dof_lcdm": float(lcdm["chi2_dof"]),
        "AIC_log": float(log["AIC_like"]),
        "AIC_lcdm": float(lcdm["AIC_like"]),
        "delta_AIC_log_minus_lcdm": float(log["AIC_like"] - lcdm["AIC_like"]),
        "BIC_log": float(log["BIC_like"]),
        "BIC_lcdm": float(lcdm["BIC_like"]),
        "delta_BIC_log_minus_lcdm": float(log["BIC_like"] - lcdm["BIC_like"]),
        "preferred_AIC": "EXACT_LOG_Pfree" if float(log["AIC_like"]) < float(lcdm["AIC_like"]) else "LCDM_Omfree",
        "cov_note": str(log["cov_note"]),
        "subset_note": str(log["subset_note"]),
    }
    row.update(extra)
    return row


def fit_pair_for_subset(
    cfg: DatasetConfig,
    z_col: str,
    args: argparse.Namespace,
    sensitivity: str,
    subset_label: str,
    extra: dict[str, object],
    extra_mask: np.ndarray | None = None,
) -> tuple[dict[str, object], pd.DataFrame]:
    prepared = prepare_masked_dataset(cfg, z_col, args, extra_mask=extra_mask, subset_note_extra=subset_label)
    fits, residuals = fit_two_models(prepared, args)
    comp = comparison_from_residuals(residuals)
    comp["sensitivity"] = sensitivity
    comp["subset"] = subset_label
    return summarize_fit_pair(fits, sensitivity, subset_label, extra), comp


def fit_pair_for_custom_z(
    cfg: DatasetConfig,
    z_label: str,
    z_values: np.ndarray,
    args: argparse.Namespace,
    sensitivity: str,
    subset_label: str,
    extra: dict[str, object],
    extra_mask: np.ndarray | None = None,
) -> tuple[dict[str, object], pd.DataFrame]:
    prepared = prepare_custom_z_dataset(
        cfg,
        z_label,
        z_values,
        args,
        extra_mask=extra_mask,
        subset_note_extra=subset_label,
    )
    fits, residuals = fit_two_models(prepared, args)
    comp = comparison_from_residuals(residuals)
    comp["sensitivity"] = sensitivity
    comp["subset"] = subset_label
    return summarize_fit_pair(fits, sensitivity, subset_label, extra), comp


def note_fit_progress(args: argparse.Namespace, fit_rows: list[dict[str, object]], row: dict[str, object]) -> None:
    every = int(getattr(args, "progress_every", 0) or 0)
    if every <= 0 or len(fit_rows) % every:
        return
    print(
        f"[fit {len(fit_rows)}] {row.get('sensitivity')} {row.get('dataset')} {row.get('z_col')} "
        f"{row.get('subset')} dAIC={float(row.get('delta_AIC_log_minus_lcdm', float('nan'))):.4g} "
        f"preferred={row.get('preferred_AIC')}",
        flush=True,
    )


def run_formal_refits(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    fit_rows: list[dict[str, object]] = []
    comp_frames: list[pd.DataFrame] = []
    errors: list[dict[str, object]] = []

    for cfg in dataset_configs(args):
        raw = read_table(cfg.data_path)
        correction = add_redshift_corrections(raw, cfg, args)
        z_cols = requested_z_columns(raw, args.z_col)
        if not z_cols:
            errors.append({"dataset": cfg.label, "z_col": "", "status": "error", "note": "No requested z columns present."})
            continue

        for z_min in parse_float_list(args.z_min_values):
            run_args = make_args(args, z_min=z_min)
            for z_col in z_cols:
                try:
                    row, comp = fit_pair_for_subset(
                        cfg,
                        z_col,
                        run_args,
                        "z_min",
                        f"z_min_{z_min:g}",
                        {"z_min": z_min, "flow_keep_fraction": np.nan, "flow_abs_v_cut_kms": np.nan},
                    )
                    fit_rows.append(row)
                    note_fit_progress(args, fit_rows, row)
                    if abs(z_min - args.object_z_min) < 1e-12:
                        comp_frames.append(comp)
                except Exception as exc:
                    errors.append({"dataset": cfg.label, "z_col": z_col, "sensitivity": "z_min", "z_min": z_min, "status": "error", "note": str(exc)})

        flow = correction["abs_v_hd_hel_kms"].to_numpy(dtype=float)
        flow_ok = np.isfinite(flow)
        for keep_fraction in parse_float_list(args.flow_keep_fractions):
            if keep_fraction <= 0.0 or keep_fraction > 1.0 or not np.any(flow_ok):
                continue
            cutoff = float(np.nanquantile(flow[flow_ok], keep_fraction))
            extra_mask = flow_ok & (flow <= cutoff)
            run_args = make_args(args, z_min=args.object_z_min)
            for z_col in z_cols:
                try:
                    row, _comp = fit_pair_for_subset(
                        cfg,
                        z_col,
                        run_args,
                        "flow_keep_low_abs_v",
                        f"keep_low_abs_v_{keep_fraction:.2f}",
                        {"z_min": args.object_z_min, "flow_keep_fraction": keep_fraction, "flow_abs_v_cut_kms": cutoff},
                        extra_mask=extra_mask,
                    )
                    fit_rows.append(row)
                    note_fit_progress(args, fit_rows, row)
                except Exception as exc:
                    errors.append(
                        {
                            "dataset": cfg.label,
                            "z_col": z_col,
                            "sensitivity": "flow_keep_low_abs_v",
                            "flow_keep_fraction": keep_fraction,
                            "status": "error",
                            "note": str(exc),
                        }
                    )

        for v_cut in parse_float_list(args.flow_v_cuts):
            if v_cut <= 0.0 or not np.any(flow_ok):
                continue
            extra_mask = flow_ok & (flow <= v_cut)
            run_args = make_args(args, z_min=args.object_z_min)
            for z_col in z_cols:
                try:
                    row, _comp = fit_pair_for_subset(
                        cfg,
                        z_col,
                        run_args,
                        "flow_abs_v_cut",
                        f"abs_v_le_{v_cut:g}_kms",
                        {"z_min": args.object_z_min, "flow_keep_fraction": np.nan, "flow_abs_v_cut_kms": v_cut},
                        extra_mask=extra_mask,
                    )
                    fit_rows.append(row)
                    note_fit_progress(args, fit_rows, row)
                except Exception as exc:
                    errors.append(
                        {
                            "dataset": cfg.label,
                            "z_col": z_col,
                            "sensitivity": "flow_abs_v_cut",
                            "flow_abs_v_cut_kms": v_cut,
                            "status": "error",
                            "note": str(exc),
                        }
                    )

        if "zHD" in raw.columns and "zHEL" in raw.columns:
            zhd = numeric(raw, "zHD")
            zhel = numeric(raw, "zHEL")
            zpair_ok = np.isfinite(zhd) & np.isfinite(zhel)
            run_args = make_args(args, z_min=args.object_z_min)
            for alpha in parse_float_list(args.blend_alpha_values):
                z_blend = redshift_blend(zhel, zhd, alpha, args.blend_mode)
                z_label = f"zHEL_to_zHD_{args.blend_mode}_alpha_{alpha:.3f}".replace(".", "p")
                try:
                    row, _comp = fit_pair_for_custom_z(
                        cfg,
                        z_label,
                        z_blend,
                        run_args,
                        "redshift_blend_hdhel",
                        f"alpha_{alpha:.3f}",
                        {
                            "z_min": args.object_z_min,
                            "flow_keep_fraction": np.nan,
                            "flow_abs_v_cut_kms": np.nan,
                            "blend_alpha_hdhel": alpha,
                            "blend_mode": args.blend_mode,
                        },
                        extra_mask=zpair_ok,
                    )
                    fit_rows.append(row)
                    note_fit_progress(args, fit_rows, row)
                except Exception as exc:
                    errors.append(
                        {
                            "dataset": cfg.label,
                            "z_col": z_label,
                            "sensitivity": "redshift_blend_hdhel",
                            "blend_alpha_hdhel": alpha,
                            "blend_mode": args.blend_mode,
                            "status": "error",
                            "note": str(exc),
                        }
                    )

        if args.leave_one_survey and "IDSURVEY" in correction.columns:
            run_args = make_args(args, z_min=args.object_z_min)
            surveys = correction.loc[correction["base_mask"], "IDSURVEY"].astype(str)
            for survey_id, count in surveys.value_counts().items():
                if int(count) < args.survey_drop_min_n:
                    continue
                extra_mask = correction["IDSURVEY"].astype(str).to_numpy() != str(survey_id)
                for z_col in z_cols:
                    try:
                        row, _comp = fit_pair_for_subset(
                            cfg,
                            z_col,
                            run_args,
                            "leave_one_survey",
                            f"drop_IDSURVEY_{survey_id}",
                            {
                                "z_min": args.object_z_min,
                                "flow_keep_fraction": np.nan,
                                "flow_abs_v_cut_kms": np.nan,
                                "dropped_IDSURVEY": survey_id,
                                "dropped_N": int(count),
                            },
                            extra_mask=extra_mask,
                        )
                        fit_rows.append(row)
                        note_fit_progress(args, fit_rows, row)
                    except Exception as exc:
                        errors.append(
                            {
                                "dataset": cfg.label,
                                "z_col": z_col,
                                "sensitivity": "leave_one_survey",
                                "dropped_IDSURVEY": survey_id,
                                "status": "error",
                                "note": str(exc),
                            }
                        )

    return pd.DataFrame(fit_rows), pd.concat(comp_frames, ignore_index=True) if comp_frames else pd.DataFrame(), pd.DataFrame(errors)


def build_correction_table(args: argparse.Namespace) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for cfg in dataset_configs(args):
        raw = read_table(cfg.data_path)
        frames.append(add_redshift_corrections(raw, cfg, args))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def merge_corrections(comp: pd.DataFrame, corrections: pd.DataFrame) -> pd.DataFrame:
    if comp.empty or corrections.empty:
        return comp
    keep = [c for c in corrections.columns if c not in {"CID", "IDSURVEY"}]
    merged = comp.merge(corrections[keep], on=["dataset", "source_index"], how="left", suffixes=("", "_corr"))
    return merged


def corr_pair(x: pd.Series, y: pd.Series) -> tuple[float, float, int]:
    xx = pd.to_numeric(x, errors="coerce")
    yy = pd.to_numeric(y, errors="coerce")
    mask = np.isfinite(xx) & np.isfinite(yy)
    n = int(mask.sum())
    if n < 8:
        return np.nan, np.nan, n
    if pd.Series(xx[mask]).nunique(dropna=True) < 2 or pd.Series(yy[mask]).nunique(dropna=True) < 2:
        return np.nan, np.nan, n
    pearson = float(np.corrcoef(xx[mask], yy[mask])[0, 1])
    if HAVE_SCIPY:
        spear = float(spearmanr(xx[mask], yy[mask]).statistic)
    else:
        spear = float(pd.Series(xx[mask]).rank().corr(pd.Series(yy[mask]).rank()))
    return pearson, spear, n


def build_correlation_audit(obj: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    targets = [
        "diag_chi2_delta_lcdm_minus_log",
        "delta_model_mu_log_minus_lcdm",
        "delta_residual_log_minus_lcdm",
        "abs_residual_improvement_log_minus_lcdm",
    ]
    for (dataset, z_col), part in obj.groupby(["dataset", "z_col"], dropna=False):
        for feature in FLOW_FEATURES:
            if feature not in part.columns:
                continue
            for target in targets:
                pearson, spearman, n = corr_pair(part[feature], part[target])
                if n:
                    rows.append(
                        {
                            "dataset": dataset,
                            "z_col": z_col,
                            "feature": feature,
                            "target": target,
                            "N": n,
                            "pearson": pearson,
                            "spearman": spearman,
                            "abs_spearman": abs(spearman) if np.isfinite(spearman) else np.nan,
                        }
                    )
    return pd.DataFrame(rows).sort_values("abs_spearman", ascending=False).reset_index(drop=True) if rows else pd.DataFrame()


def summarize_corrections(corrections: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    survey_rows: list[dict[str, object]] = []
    for dataset, part in corrections[corrections["base_mask"]].groupby("dataset", dropna=False):
        rows.append(correction_summary_row(part, {"dataset": dataset}))
        for survey_id, grp in part.groupby("IDSURVEY", dropna=False):
            row = correction_summary_row(grp, {"dataset": dataset, "IDSURVEY": survey_id})
            survey_rows.append(row)
    return pd.DataFrame(rows), pd.DataFrame(survey_rows)


def correction_summary_row(part: pd.DataFrame, keys: dict[str, object]) -> dict[str, object]:
    row = dict(keys)
    row["N"] = int(len(part))
    for col in ["zHD", "zCMB", "zHEL"]:
        if col in part.columns:
            vals = pd.to_numeric(part[col], errors="coerce").to_numpy(dtype=float)
            vals = vals[np.isfinite(vals)]
            row[f"{col}_median"] = float(np.median(vals)) if vals.size else np.nan
    for col in ["dz_hd_hel", "v_hd_hel_kms", "abs_v_hd_hel_kms", "dz_hd_cmb", "v_hd_cmb_kms", "abs_v_hd_cmb_kms", "dz_cmb_hel", "v_cmb_hel_kms", "abs_v_cmb_hel_kms", "MUERR_VPEC", "VPEC"]:
        if col in part.columns:
            vals = pd.to_numeric(part[col], errors="coerce").to_numpy(dtype=float)
            vals = vals[np.isfinite(vals)]
            if vals.size:
                row[f"{col}_mean"] = float(np.mean(vals))
                row[f"{col}_median"] = float(np.median(vals))
                row[f"{col}_p90"] = float(np.percentile(vals, 90))
                row[f"{col}_max"] = float(np.max(vals))
    return row


def build_flow_strata(obj: pd.DataFrame, bins: int) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if obj.empty or "abs_v_hd_hel_kms" not in obj.columns:
        return pd.DataFrame()
    for (dataset, z_col), part in obj.groupby(["dataset", "z_col"], dropna=False):
        valid = part[np.isfinite(part["abs_v_hd_hel_kms"])].copy()
        if valid.empty:
            continue
        try:
            valid["flow_bin"] = pd.qcut(valid["abs_v_hd_hel_kms"], q=min(bins, len(valid)), duplicates="drop").astype(str)
        except ValueError:
            valid["flow_bin"] = "all"
        for flow_bin, grp in valid.groupby("flow_bin", dropna=False):
            rows.append(
                {
                    "dataset": dataset,
                    "z_col": z_col,
                    "flow_bin": flow_bin,
                    "N": int(len(grp)),
                    "z_median": float(np.nanmedian(grp["z"])),
                    "abs_v_hd_hel_kms_median": float(np.nanmedian(grp["abs_v_hd_hel_kms"])),
                    "abs_v_hd_hel_kms_max": float(np.nanmax(grp["abs_v_hd_hel_kms"])),
                    "diag_chi2_delta_lcdm_minus_log_sum": float(grp["diag_chi2_delta_lcdm_minus_log"].sum(skipna=True)),
                    "diag_chi2_delta_lcdm_minus_log_mean": float(grp["diag_chi2_delta_lcdm_minus_log"].mean(skipna=True)),
                    "log_better_diag_fraction": float(grp["log_better_diag"].mean()),
                    "delta_model_mu_log_minus_lcdm_mean": float(grp["delta_model_mu_log_minus_lcdm"].mean(skipna=True)),
                }
            )
    return pd.DataFrame(rows)


def build_survey_flow_audit(obj: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if obj.empty:
        return pd.DataFrame()
    for (dataset, z_col, survey_id), grp in obj.groupby(["dataset", "z_col", "IDSURVEY"], dropna=False):
        rows.append(
            {
                "dataset": dataset,
                "z_col": z_col,
                "IDSURVEY": survey_id,
                "N": int(len(grp)),
                "z_median": float(np.nanmedian(grp["z"])),
                "abs_v_hd_hel_kms_median": float(np.nanmedian(grp["abs_v_hd_hel_kms"])) if "abs_v_hd_hel_kms" in grp.columns else np.nan,
                "abs_v_hd_hel_kms_p90": float(np.nanpercentile(grp["abs_v_hd_hel_kms"], 90)) if "abs_v_hd_hel_kms" in grp.columns else np.nan,
                "diag_chi2_delta_lcdm_minus_log_sum": float(grp["diag_chi2_delta_lcdm_minus_log"].sum(skipna=True)),
                "log_better_diag_fraction": float(grp["log_better_diag"].mean()),
            }
        )
    if not rows:
        return pd.DataFrame()
    out = pd.DataFrame(rows)
    out["abs_diag_delta"] = out["diag_chi2_delta_lcdm_minus_log_sum"].abs()
    return out.sort_values(["dataset", "z_col", "abs_diag_delta"], ascending=[True, True, False]).reset_index(drop=True)


def safe_name(*parts: object) -> str:
    raw = "_".join(str(p) for p in parts)
    return "".join(c if c.isalnum() or c in {"_", "-", "."} else "_" for c in raw)


def plot_zmin_sensitivity(fits: pd.DataFrame, outdir: Path) -> list[Path]:
    paths: list[Path] = []
    part = fits[fits["sensitivity"] == "z_min"].copy()
    if part.empty:
        return paths
    fig, ax = plt.subplots(figsize=(9, 5.2))
    for (dataset, z_col), grp in part.groupby(["dataset", "z_col"], dropna=False):
        grp = grp.sort_values("z_min")
        ax.plot(grp["z_min"], grp["delta_AIC_log_minus_lcdm"], marker="o", label=f"{dataset} {z_col}")
    ax.axhline(0.0, color="black", lw=0.8)
    ax.set_xlabel("minimum redshift cut")
    ax.set_ylabel("delta AIC: exact log - LCDM")
    ax.set_title("Local-flow sensitivity through low-z cuts")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    path = outdir / "local_flow_zmin_delta_aic.png"
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    paths.append(path)
    return paths


def plot_flow_keep_sensitivity(fits: pd.DataFrame, outdir: Path) -> list[Path]:
    paths: list[Path] = []
    part = fits[fits["sensitivity"] == "flow_keep_low_abs_v"].copy()
    if part.empty:
        return paths
    fig, ax = plt.subplots(figsize=(9, 5.2))
    for (dataset, z_col), grp in part.groupby(["dataset", "z_col"], dropna=False):
        grp = grp.sort_values("flow_keep_fraction")
        ax.plot(grp["flow_keep_fraction"], grp["delta_AIC_log_minus_lcdm"], marker="o", label=f"{dataset} {z_col}")
    ax.axhline(0.0, color="black", lw=0.8)
    ax.set_xlabel("fraction retained with smallest |zHD-zHEL| velocity correction")
    ax.set_ylabel("delta AIC: exact log - LCDM")
    ax.set_title("Sensitivity to high local-flow correction objects")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    path = outdir / "local_flow_keep_fraction_delta_aic.png"
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    paths.append(path)
    return paths


def plot_flow_vcut_sensitivity(fits: pd.DataFrame, outdir: Path) -> list[Path]:
    paths: list[Path] = []
    part = fits[fits["sensitivity"] == "flow_abs_v_cut"].copy()
    if part.empty:
        return paths
    fig, ax = plt.subplots(figsize=(9, 5.2))
    for (dataset, z_col), grp in part.groupby(["dataset", "z_col"], dropna=False):
        grp = grp.sort_values("flow_abs_v_cut_kms")
        ax.plot(grp["flow_abs_v_cut_kms"], grp["delta_AIC_log_minus_lcdm"], marker="o", label=f"{dataset} {z_col}")
    ax.axhline(0.0, color="black", lw=0.8)
    ax.set_xlabel("maximum |zHD-zHEL| velocity correction retained [km/s]")
    ax.set_ylabel("delta AIC: exact log - LCDM")
    ax.set_title("Sensitivity to absolute local-flow velocity cuts")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    path = outdir / "local_flow_vcut_delta_aic.png"
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    paths.append(path)
    return paths


def plot_blend_sensitivity(fits: pd.DataFrame, outdir: Path) -> list[Path]:
    paths: list[Path] = []
    part = fits[fits["sensitivity"] == "redshift_blend_hdhel"].copy()
    if part.empty:
        return paths
    fig, ax = plt.subplots(figsize=(9, 5.2))
    for (dataset, blend_mode), grp in part.groupby(["dataset", "blend_mode"], dropna=False):
        grp = grp.sort_values("blend_alpha_hdhel")
        ax.plot(grp["blend_alpha_hdhel"], grp["delta_AIC_log_minus_lcdm"], marker="o", label=f"{dataset} {blend_mode}")
    ax.axhline(0.0, color="black", lw=0.8)
    ax.set_xlabel("alpha: additive or multiplicative transition from zHEL to zHD")
    ax.set_ylabel("delta AIC: exact log - LCDM")
    ax.set_title("Preference transition from heliocentric to corrected redshift")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    path = outdir / "local_flow_redshift_blend_delta_aic.png"
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    paths.append(path)
    return paths


def plot_object_flow_scatter(obj: pd.DataFrame, outdir: Path) -> list[Path]:
    paths: list[Path] = []
    if obj.empty or "abs_v_hd_hel_kms" not in obj.columns:
        return paths
    for (dataset, z_col), part in obj.groupby(["dataset", "z_col"], dropna=False):
        valid = part[np.isfinite(part["abs_v_hd_hel_kms"])].copy()
        if valid.empty:
            continue
        fig, ax = plt.subplots(figsize=(8.5, 5))
        ax.scatter(valid["abs_v_hd_hel_kms"], valid["diag_chi2_delta_lcdm_minus_log"], s=8, alpha=0.35)
        ax.axhline(0.0, color="black", lw=0.8)
        ax.set_xlabel("approx |zHD-zHEL| velocity correction [km/s]")
        ax.set_ylabel("diag chi2 improvement: LCDM - log")
        ax.set_title(f"{dataset} {z_col}: row leverage vs local-flow correction")
        ax.grid(alpha=0.25)
        path = outdir / f"object_flow_scatter_{safe_name(dataset, z_col)}.png"
        fig.tight_layout()
        fig.savefig(path, dpi=160)
        plt.close(fig)
        paths.append(path)
    return paths


def plot_survey_flow(survey: pd.DataFrame, outdir: Path) -> list[Path]:
    paths: list[Path] = []
    if survey.empty:
        return paths
    for (dataset, z_col), part in survey.groupby(["dataset", "z_col"], dropna=False):
        show = part.sort_values("abs_diag_delta", ascending=False).head(16).sort_values("diag_chi2_delta_lcdm_minus_log_sum")
        fig, ax = plt.subplots(figsize=(9, max(4.8, 0.32 * len(show))))
        colors = np.where(show["diag_chi2_delta_lcdm_minus_log_sum"] >= 0, "#b44c43", "#3b6ea8")
        labels = show["IDSURVEY"].astype(str) + "  |v|~" + show["abs_v_hd_hel_kms_median"].round(0).astype(str)
        ax.barh(labels, show["diag_chi2_delta_lcdm_minus_log_sum"], color=colors)
        ax.axvline(0.0, color="black", lw=0.8)
        ax.set_xlabel("diag chi2 improvement: LCDM - log")
        ax.set_ylabel("IDSURVEY and median |vHD-HEL|")
        ax.set_title(f"{dataset} {z_col}: survey leverage with correction scale")
        ax.grid(axis="x", alpha=0.25)
        path = outdir / f"survey_flow_leverage_{safe_name(dataset, z_col)}.png"
        fig.tight_layout()
        fig.savefig(path, dpi=160)
        plt.close(fig)
        paths.append(path)
    return paths


def markdown_table(df: pd.DataFrame, cols: list[str], max_rows: int = 20, floatfmt: str = ".6g") -> str:
    if df.empty:
        return "_No rows._"
    show = df[[c for c in cols if c in df.columns]].head(max_rows).copy()
    return show.to_markdown(index=False, floatfmt=floatfmt)


def write_report(
    path: Path,
    fits: pd.DataFrame,
    correction_summary: pd.DataFrame,
    survey_correction: pd.DataFrame,
    object_diag: pd.DataFrame,
    corr: pd.DataFrame,
    strata: pd.DataFrame,
    survey_flow: pd.DataFrame,
    errors: pd.DataFrame,
    plots: list[Path],
    args: argparse.Namespace,
) -> None:
    lines: list[str] = []
    lines.append("# PLAMB Local-Flow / Redshift-Correction Audit")
    lines.append("")
    lines.append("This audit tests whether the exact logarithmic PLAMB kernel preference follows the redshift-correction machinery rather than a robust cosmological shape.")
    lines.append("")
    lines.append("## Configuration")
    lines.append("")
    lines.append(f"- datasets: `{args.dataset}`")
    lines.append(f"- z columns requested: `{args.z_col}`")
    lines.append(f"- H0 fixed in raw-MU kernel fits: `{args.H0_value}`")
    lines.append(f"- low-z cuts: `{args.z_min_values}`")
    lines.append(f"- flow keep fractions: `{args.flow_keep_fractions}` at z_min `{args.object_z_min}`")
    lines.append(f"- absolute flow velocity cuts: `{args.flow_v_cuts or 'not run'}` at z_min `{args.object_z_min}`")
    lines.append(f"- zHEL to zHD blend alphas: `{args.blend_alpha_values or 'not run'}` at z_min `{args.object_z_min}`")
    lines.append(f"- zHEL to zHD blend mode: `{args.blend_mode}`")
    lines.append(f"- leave-one-survey refits: `{args.leave_one_survey}`")
    lines.append("- velocity correction proxy: `c * (zHD - zHEL)/(1+zHEL)` where columns exist")
    lines.append("")
    lines.append("## Formal Low-z and Redshift-Basis Fits")
    lines.append("")
    lines.append("Negative delta AIC favors the exact log kernel; positive favors LCDM.")
    lines.append("")
    zfit = fits[fits["sensitivity"] == "z_min"].sort_values(["dataset", "z_col", "z_min"]) if not fits.empty else pd.DataFrame()
    lines.append(markdown_table(zfit, ["dataset", "z_col", "z_min", "N", "log_p", "lcdm_Om", "delta_AIC_log_minus_lcdm", "preferred_AIC"], max_rows=60))
    lines.append("")
    lines.append("## High-Correction Object Removal")
    lines.append("")
    lines.append("These rows keep only the fraction of objects with the smallest absolute `zHD-zHEL` velocity correction. If the preference changes sharply here, local-flow corrections are implicated.")
    lines.append("")
    flowfit = fits[fits["sensitivity"] == "flow_keep_low_abs_v"].sort_values(["dataset", "z_col", "flow_keep_fraction"]) if not fits.empty else pd.DataFrame()
    lines.append(markdown_table(flowfit, ["dataset", "z_col", "flow_keep_fraction", "flow_abs_v_cut_kms", "N", "delta_AIC_log_minus_lcdm", "preferred_AIC"], max_rows=60))
    lines.append("")
    lines.append("## Absolute Velocity-Cut Fits")
    lines.append("")
    vcut = fits[fits["sensitivity"] == "flow_abs_v_cut"].sort_values(["dataset", "z_col", "flow_abs_v_cut_kms"]) if not fits.empty else pd.DataFrame()
    lines.append(markdown_table(vcut, ["dataset", "z_col", "flow_abs_v_cut_kms", "N", "delta_AIC_log_minus_lcdm", "preferred_AIC"], max_rows=80))
    lines.append("")
    lines.append("## Redshift-Correction Blend Fits")
    lines.append("")
    lines.append("Here alpha=0 uses zHEL and alpha=1 uses zHD. Additive mode uses zHEL + alpha (zHD-zHEL); multiplicative mode uses (1+z) = (1+zHEL) [(1+zHD)/(1+zHEL)]^alpha. A smooth change in delta AIC across alpha means the apparent model preference is tied to the redshift-correction convention.")
    lines.append("")
    blend = fits[fits["sensitivity"] == "redshift_blend_hdhel"].sort_values(["dataset", "blend_alpha_hdhel"]) if not fits.empty else pd.DataFrame()
    lines.append(markdown_table(blend, ["dataset", "blend_mode", "blend_alpha_hdhel", "N", "log_p", "lcdm_Om", "delta_AIC_log_minus_lcdm", "preferred_AIC"], max_rows=80))
    lines.append("")
    lines.append("## Leave-One-Survey Formal Fits")
    lines.append("")
    loso = fits[fits["sensitivity"] == "leave_one_survey"].copy() if not fits.empty else pd.DataFrame()
    if not loso.empty:
        loso["abs_delta_AIC_log_minus_lcdm"] = pd.to_numeric(loso["delta_AIC_log_minus_lcdm"], errors="coerce").abs()
        loso = loso.sort_values(["dataset", "z_col", "abs_delta_AIC_log_minus_lcdm"], ascending=[True, True, False])
    lines.append(markdown_table(loso, ["dataset", "z_col", "dropped_IDSURVEY", "dropped_N", "N", "delta_AIC_log_minus_lcdm", "preferred_AIC"], max_rows=80))
    lines.append("")
    lines.append("## Redshift-Correction Scale")
    lines.append("")
    lines.append(markdown_table(correction_summary, ["dataset", "N", "zHD_median", "zCMB_median", "zHEL_median", "abs_v_hd_hel_kms_median", "abs_v_hd_hel_kms_p90", "abs_v_hd_hel_kms_max", "MUERR_VPEC_median", "VPEC_median"], max_rows=20))
    lines.append("")
    lines.append("## Dominant Survey Blocks")
    lines.append("")
    lines.append(markdown_table(survey_flow.sort_values("abs_diag_delta", ascending=False) if not survey_flow.empty else survey_flow, ["dataset", "z_col", "IDSURVEY", "N", "z_median", "abs_v_hd_hel_kms_median", "abs_v_hd_hel_kms_p90", "diag_chi2_delta_lcdm_minus_log_sum", "log_better_diag_fraction"], max_rows=30))
    lines.append("")
    lines.append("## Strongest Correction Correlations")
    lines.append("")
    lines.append(markdown_table(corr, ["dataset", "z_col", "feature", "target", "N", "pearson", "spearman", "abs_spearman"], max_rows=30))
    lines.append("")
    lines.append("## Row-Level Chi2 Correlations Only")
    lines.append("")
    lines.append("This table is the stricter diagnostic: it asks whether redshift-correction variables track which model actually improves the approximate per-object chi2.")
    lines.append("")
    diag_corr = (
        corr[corr["target"] == "diag_chi2_delta_lcdm_minus_log"].sort_values("abs_spearman", ascending=False)
        if not corr.empty and "target" in corr.columns
        else pd.DataFrame()
    )
    lines.append(markdown_table(diag_corr, ["dataset", "z_col", "feature", "N", "pearson", "spearman", "abs_spearman"], max_rows=30))
    lines.append("")
    lines.append("## Flow-Correction Strata")
    lines.append("")
    lines.append(markdown_table(strata, ["dataset", "z_col", "flow_bin", "N", "z_median", "abs_v_hd_hel_kms_median", "diag_chi2_delta_lcdm_minus_log_sum", "log_better_diag_fraction"], max_rows=60))
    lines.append("")
    lines.append("## Interpretation")
    lines.append("")
    lines.append("- A robust cosmological PLAMB preference should survive switching between `zHD`, `zCMB`, and `zHEL`, and should not disappear under conservative low-z cuts.")
    lines.append("- If the log-kernel preference is strongest in `zHEL` and weakens or flips under `z_min >= 0.03` or `0.05`, the result is local-flow/redshift-convention sensitive.")
    lines.append("- Survey blocks with high correction leverage should be treated as calibration diagnostics before they are treated as cosmology.")
    lines.append("- The diagonal row/block diagnostics identify where the signal lives; the full covariance rows above decide the formal ranking.")
    if not errors.empty:
        lines.append("")
        lines.append("## Errors")
        lines.append("")
        lines.append(markdown_table(errors, list(errors.columns), max_rows=50))
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
    parser.add_argument("--des-data", type=Path, default=DEFAULT_DES)
    parser.add_argument("--des-cov", type=Path, default=DEFAULT_DES_COV)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--dataset", default="all", help="all, PantheonPlusSH0ES, DES_Dovekie_raw, or comma-separated.")
    parser.add_argument("--z-col", default="all", help="z column, comma-separated list, or all for zHD,zCMB,zHEL where present.")
    parser.add_argument("--H0-value", type=float, default=67.5)
    parser.add_argument("--z-min", type=float, default=0.01)
    parser.add_argument("--z-min-values", default="0.01,0.02,0.03,0.05")
    parser.add_argument("--object-z-min", type=float, default=0.01)
    parser.add_argument("--flow-keep-fractions", default="1.0,0.9,0.8")
    parser.add_argument("--flow-v-cuts", default="", help="Comma-separated absolute |zHD-zHEL| velocity cuts in km/s.")
    parser.add_argument("--blend-alpha-values", default="", help="Comma-separated alpha values for z=(1-alpha)zHEL+alpha zHD.")
    parser.add_argument("--blend-mode", choices=["additive", "multiplicative"], default="additive")
    parser.add_argument("--leave-one-survey", action="store_true", help="Run formal leave-one-IDSURVEY refits at object-z-min.")
    parser.add_argument("--survey-drop-min-n", type=int, default=8)
    parser.add_argument("--max-z", type=float, default=None)
    parser.add_argument("--min-n", type=int, default=10)
    parser.add_argument("--keep-calibrators", action="store_true")
    parser.add_argument("--allow-precision-submatrix", action="store_true")
    parser.add_argument("--grid-steps", type=int, default=401)
    parser.add_argument("--maxiter", type=int, default=80)
    parser.add_argument("--flow-bins", type=int, default=8)
    parser.add_argument("--progress-every", type=int, default=25, help="Print one progress line after this many formal refits; 0 disables.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    plot_dir = args.outdir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    corrections = build_correction_table(args)
    fits, base_comp, errors = run_formal_refits(args)
    object_diag = merge_corrections(base_comp, corrections)
    corr = build_correlation_audit(object_diag)
    correction_summary, survey_correction = summarize_corrections(corrections)
    strata = build_flow_strata(object_diag, args.flow_bins)
    survey_flow = build_survey_flow_audit(object_diag)

    plots: list[Path] = []
    plots.extend(plot_zmin_sensitivity(fits, plot_dir))
    plots.extend(plot_flow_keep_sensitivity(fits, plot_dir))
    plots.extend(plot_flow_vcut_sensitivity(fits, plot_dir))
    plots.extend(plot_blend_sensitivity(fits, plot_dir))
    plots.extend(plot_object_flow_scatter(object_diag, plot_dir))
    plots.extend(plot_survey_flow(survey_flow, plot_dir))

    fit_path = args.outdir / "plamb_local_flow_fit_sensitivity.csv"
    corr_table_path = args.outdir / "plamb_local_flow_redshift_corrections.csv"
    correction_summary_path = args.outdir / "plamb_local_flow_correction_summary.csv"
    survey_correction_path = args.outdir / "plamb_local_flow_survey_correction_summary.csv"
    object_path = args.outdir / "plamb_local_flow_object_diagnostics.csv"
    correlation_path = args.outdir / "plamb_local_flow_correlation_audit.csv"
    strata_path = args.outdir / "plamb_local_flow_strata.csv"
    survey_flow_path = args.outdir / "plamb_local_flow_survey_leverage.csv"
    errors_path = args.outdir / "plamb_local_flow_errors.csv"
    report_path = args.outdir / "plamb_local_flow_redshift_audit.md"
    config_path = args.outdir / "plamb_local_flow_redshift_config.json"

    fits.to_csv(fit_path, index=False)
    corrections.to_csv(corr_table_path, index=False)
    correction_summary.to_csv(correction_summary_path, index=False)
    survey_correction.to_csv(survey_correction_path, index=False)
    object_diag.to_csv(object_path, index=False)
    corr.to_csv(correlation_path, index=False)
    strata.to_csv(strata_path, index=False)
    survey_flow.to_csv(survey_flow_path, index=False)
    errors.to_csv(errors_path, index=False)
    config_path.write_text(json.dumps(vars(args), indent=2, default=str), encoding="utf-8")
    write_report(report_path, fits, correction_summary, survey_correction, object_diag, corr, strata, survey_flow, errors, plots, args)

    print(f"Saved fit sensitivity: {fit_path}")
    print(f"Saved redshift corrections: {corr_table_path}")
    print(f"Saved correction summary: {correction_summary_path}")
    print(f"Saved survey correction summary: {survey_correction_path}")
    print(f"Saved object diagnostics: {object_path}")
    print(f"Saved correlation audit: {correlation_path}")
    print(f"Saved flow strata: {strata_path}")
    print(f"Saved survey leverage: {survey_flow_path}")
    print(f"Saved report: {report_path}")
    print(f"Saved plots: {plot_dir}")
    if not fits.empty:
        zfit = fits[fits["sensitivity"] == "z_min"].copy()
        for _, row in zfit.sort_values(["dataset", "z_col", "z_min"]).iterrows():
            print(
                f"{row['dataset']} {row['z_col']} z_min={float(row['z_min']):.3g}: "
                f"dAIC(log-LCDM)={float(row['delta_AIC_log_minus_lcdm']):.3f} "
                f"preferred={row['preferred_AIC']}"
            )


if __name__ == "__main__":
    main()
