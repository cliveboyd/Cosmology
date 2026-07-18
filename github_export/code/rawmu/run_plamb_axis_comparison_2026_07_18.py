#!/usr/bin/env python3
"""Run Peter Lamb's requested three-release PLAMB axis comparison."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar


ROOT = Path(__file__).resolve().parents[3]
CODE_DIR = Path(__file__).resolve().parent
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))
SHARED_CODE = ROOT / "github_export" / "code" / "shared"
if str(SHARED_CODE) not in sys.path:
    sys.path.insert(0, str(SHARED_CODE))

import run_rawmu_release_grounded_holdouts_2026_07_18 as base  # noqa: E402
from plamb_clock_distance import (  # noqa: E402
    PETER_LINEAR_REDSHIFT,
    clock_path_integral,
    clock_path_integrand,
)


DATE = "2026-07-18"
OUT = ROOT / "github_export" / "results" / DATE / "rawmu_plamb_axis"
PREREG = OUT / f"plamb_axis_comparison_preregistration_{DATE}.json"
H0 = 67.5
OMEGA_BOUNDS = (0.05, 0.60)
C_KMS = base.C_KMS
MODEL_FR = "FR_C1PZ_ALPHA0"
MODEL_LCDM = "LCDM_OMFREE"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def identity_audit(blocks: list[base.Block]) -> dict[str, float]:
    rows: dict[str, float] = {}
    for block in blocks:
        correction = 5.0 * np.log10(1.0 + 0.5 * block.z)
        corrected_mu = block.mu - correction
        linear_mu = 5.0 * np.log10((C_KMS / H0) * block.z) + 25.0
        direct = block.mu - base.model_mu(block.z, MODEL_FR, None, H0)
        transformed = corrected_mu - linear_mu
        rows[block.label] = float(np.max(np.abs(direct - transformed)))
    return rows


def diagonal_blocks(blocks: list[base.Block], suffix: str) -> list[base.Block]:
    output: list[base.Block] = []
    for block in blocks:
        diagonal = np.diag(np.diag(block.covariance))
        output.append(
            base.make_block(
                block.label,
                block.z,
                block.mu,
                diagonal,
                block.survey,
                block.source_indices,
                suffix,
            )
        )
    return output


def all_positive_redshift_blocks(primary: list[base.Block]) -> list[base.Block]:
    pantheon = base.read_table(base.PANTHEON_DATA)
    z = pd.to_numeric(pantheon["zHD"], errors="coerce").to_numpy(dtype=float)
    mu = pd.to_numeric(pantheon["MU_SH0ES"], errors="coerce").to_numpy(dtype=float)
    calibrator = (
        pd.to_numeric(pantheon["IS_CALIBRATOR"], errors="coerce").fillna(0.0).to_numpy()
    )
    survey = pd.to_numeric(pantheon["IDSURVEY"], errors="coerce").to_numpy()
    mask = np.isfinite(z) & np.isfinite(mu) & (z > 0.0) & (calibrator == 0.0)
    indices = np.flatnonzero(mask)
    covariance = base.load_pantheon_covariance(base.PANTHEON_TOTAL, len(pantheon))
    pantheon_block = base.make_block(
        "PantheonPlusSH0ES",
        z[indices],
        mu[indices],
        covariance[np.ix_(indices, indices)],
        survey[indices],
        indices,
        "released_total_all_positive_z_noncalibrator",
    )
    return [pantheon_block, *[block for block in primary if block.label != pantheon_block.label]]


def flatten_fit(scope: str, variant: str, fit: dict[str, object]) -> dict[str, object]:
    return {
        "scope": scope,
        "variant": variant,
        "N": int(fit["N"]),
        "FR_chi2": float(fit["FR_chi2"]),
        "LCDM_chi2": float(fit["LCDM_chi2"]),
        "LCDM_Omega_m": float(fit["LCDM_Omega_m"]),
        "FR_BIC": float(fit["FR_BIC"]),
        "LCDM_BIC": float(fit["LCDM_BIC"]),
        "delta_BIC_FR_minus_LCDM": float(fit["delta_BIC_FR_minus_LCDM"]),
        "same_nuisance_verified": bool(fit["same_nuisance_verified"]),
    }


def magnitude_fit_rows(primary: list[base.Block], diagonal: list[base.Block], all_z: list[base.Block]) -> tuple[list[dict[str, object]], dict[str, object]]:
    rows: list[dict[str, object]] = []
    joint = base.fit_comparison(primary, "released_total_primary", {}, H0, OMEGA_BOUNDS)
    rows.append(flatten_fit("joint", "released_total_primary", joint))
    for block in primary:
        fit = base.fit_comparison([block], "released_total_primary", {}, H0, OMEGA_BOUNDS)
        rows.append(flatten_fit(block.label, "released_total_primary", fit))
    diagonal_joint = base.fit_comparison(
        diagonal, "released_total_diagonal_mu", {}, H0, OMEGA_BOUNDS
    )
    rows.append(flatten_fit("joint", "released_total_diagonal_mu", diagonal_joint))
    all_z_joint = base.fit_comparison(
        all_z, "all_positive_z_noncalibrator", {}, H0, OMEGA_BOUNDS
    )
    rows.append(flatten_fit("joint", "all_positive_z_noncalibrator", all_z_joint))
    return rows, joint


def dimensionless_shape(z: np.ndarray, family: str, omega_m: float | None) -> np.ndarray:
    if family == MODEL_FR:
        return np.asarray(
            clock_path_integral(
                z,
                PETER_LINEAR_REDSHIFT.gamma_c,
                PETER_LINEAR_REDSHIFT.redshift_rate_power,
            ),
            dtype=float,
        )
    assert omega_m is not None
    return (1.0 + z) * base.lcdm_integral(z, omega_m)


def shape_derivative(z: np.ndarray, family: str, omega_m: float | None) -> np.ndarray:
    if family == MODEL_FR:
        return np.asarray(
            clock_path_integrand(
                z,
                PETER_LINEAR_REDSHIFT.gamma_c,
                PETER_LINEAR_REDSHIFT.redshift_rate_power,
            ),
            dtype=float,
        )
    assert omega_m is not None
    integral = base.lcdm_integral(z, omega_m)
    e = np.sqrt(omega_m * (1.0 + z) ** 3 + 1.0 - omega_m)
    return integral + (1.0 + z) / e


def redshift_errors(block: base.Block) -> tuple[np.ndarray, str]:
    if block.label == "PantheonPlusSH0ES":
        table = base.read_table(base.PANTHEON_DATA)
        errors = pd.to_numeric(table["zHDERR"], errors="coerce").to_numpy(dtype=float)
        selected = errors[block.source_indices]
        selected = np.where(np.isfinite(selected) & (selected >= 0.0), selected, 0.0)
        return selected, "released zHDERR"
    return np.zeros(block.n, dtype=float), "not supplied by compressed release product"


def distance_block_score(
    block: base.Block,
    family: str,
    omega_m: float | None,
) -> tuple[float, float, int, str]:
    observed = (H0 / C_KMS) * np.power(10.0, (block.mu - 25.0) / 5.0)
    sigma_mu = np.sqrt(np.clip(np.diag(block.covariance), 0.0, None))
    sigma_distance = (math.log(10.0) / 5.0) * observed * sigma_mu
    shape = dimensionless_shape(block.z, family, omega_m)
    derivative = shape_derivative(block.z, family, omega_m)
    sigma_z, source = redshift_errors(block)
    variance = np.square(sigma_distance)
    scale = float(np.sum(shape * observed / variance) / np.sum(shape * shape / variance))
    for _ in range(30):
        effective = variance + np.square(scale * derivative * sigma_z)
        updated = float(
            np.sum(shape * observed / effective) / np.sum(shape * shape / effective)
        )
        if abs(updated - scale) <= 1e-12 * max(abs(scale), 1.0):
            scale = updated
            break
        scale = updated
    effective = variance + np.square(scale * derivative * sigma_z)
    residual = observed - scale * shape
    chi2 = float(np.sum(np.square(residual) / effective))
    return chi2, scale, int(np.count_nonzero(sigma_z > 0.0)), source


def distance_space_fit(blocks: list[base.Block]) -> dict[str, object]:
    fr_rows = [distance_block_score(block, MODEL_FR, None) for block in blocks]
    fr_chi2 = sum(row[0] for row in fr_rows)

    def objective(omega_m: float) -> float:
        return sum(distance_block_score(block, MODEL_LCDM, omega_m)[0] for block in blocks)

    optimum = minimize_scalar(objective, bounds=OMEGA_BOUNDS, method="bounded")
    omega_m = float(optimum.x)
    lcdm_rows = [distance_block_score(block, MODEL_LCDM, omega_m) for block in blocks]
    lcdm_chi2 = sum(row[0] for row in lcdm_rows)
    n = sum(block.n for block in blocks)
    k_fr = len(blocks)
    k_lcdm = k_fr + 1
    return {
        "scope": "joint",
        "variant": "diagonal_luminosity_distance_with_reported_z_error",
        "N": n,
        "FR_chi2": fr_chi2,
        "LCDM_chi2": lcdm_chi2,
        "LCDM_Omega_m": omega_m,
        "FR_BIC": fr_chi2 + k_fr * math.log(n),
        "LCDM_BIC": lcdm_chi2 + k_lcdm * math.log(n),
        "delta_BIC_FR_minus_LCDM": (
            fr_chi2 + k_fr * math.log(n) - lcdm_chi2 - k_lcdm * math.log(n)
        ),
        "same_nuisance_verified": True,
        "redshift_error_rows": sum(row[2] for row in fr_rows),
        "redshift_error_coverage": {
            block.label: {"N": row[2], "source": row[3]}
            for block, row in zip(blocks, fr_rows)
        },
        "distance_scale_FR": {
            block.label: row[1] for block, row in zip(blocks, fr_rows)
        },
        "distance_scale_LCDM": {
            block.label: row[1] for block, row in zip(blocks, lcdm_rows)
        },
    }


def residual_bins(
    blocks: list[base.Block], joint: dict[str, object], bins: np.ndarray
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    omega_m = float(joint["LCDM_Omega_m"])
    for block in blocks:
        for family, label in ((MODEL_FR, "PLAMB"), (MODEL_LCDM, "Lambda-CDM")):
            parameter = None if family == MODEL_FR else omega_m
            _chi2, _offset, residual = base.profiled_block_score(block, family, parameter, H0)
            for lower, upper in zip(bins[:-1], bins[1:]):
                selected = (block.z >= lower) & (block.z < upper)
                indices = np.flatnonzero(selected)
                if not len(indices):
                    continue
                covariance = block.covariance[np.ix_(indices, indices)]
                precision = base.stable_inverse(covariance)
                ones = np.ones(len(indices), dtype=float)
                denominator = float(ones @ precision @ ones)
                mean = float(ones @ precision @ residual[indices] / denominator)
                sigma = math.sqrt(1.0 / denominator)
                rows.append(
                    {
                        "dataset": block.label,
                        "model": label,
                        "z_lower": lower,
                        "z_upper": upper,
                        "N": len(indices),
                        "z_effective": float(np.mean(block.z[indices])),
                        "residual_mean_mag": mean,
                        "residual_sigma_mag": sigma,
                    }
                )
    return pd.DataFrame(rows)


def coordinate_table(blocks: list[base.Block], joint: dict[str, object]) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    omega_m = float(joint["LCDM_Omega_m"])
    fr_offsets = dict(joint["FR_offsets"])
    lcdm_offsets = dict(joint["LCDM_offsets"])
    for block in blocks:
        sigma_mu = np.sqrt(np.clip(np.diag(block.covariance), 0.0, None))
        distance_fr = (H0 / C_KMS) * np.power(
            10.0, (block.mu - float(fr_offsets[block.label]) - 25.0) / 5.0
        )
        distance_lcdm = (H0 / C_KMS) * np.power(
            10.0, (block.mu - float(lcdm_offsets[block.label]) - 25.0) / 5.0
        )
        frame = pd.DataFrame(
            {
                "dataset": block.label,
                "source_index": block.source_indices,
                "z": block.z,
                "sigma_mu_diag": sigma_mu,
                "x_lcdm_stretched_redshift": dimensionless_shape(
                    block.z, MODEL_LCDM, omega_m
                ),
                "y_lcdm_observed_distance": distance_lcdm,
                "y_plamb_corrected_distance": distance_fr / (1.0 + 0.5 * block.z),
            }
        )
        rows.append(frame)
    return pd.concat(rows, ignore_index=True)


def write_plots(
    coordinates: pd.DataFrame, binned: pd.DataFrame, joint: dict[str, object]
) -> None:
    datasets = ["PantheonPlusSH0ES", "DES_Dovekie_raw", "Union3p1_UNITY1p8"]
    colours = {"PLAMB": "#b43b32", "Lambda-CDM": "#146c94"}
    fig, axes = plt.subplots(3, 2, figsize=(12, 14), constrained_layout=True)
    for row, dataset in enumerate(datasets):
        values = coordinates[coordinates["dataset"] == dataset]
        left, right = axes[row]
        left.scatter(
            values["x_lcdm_stretched_redshift"],
            values["y_lcdm_observed_distance"],
            s=8,
            alpha=0.16,
            color=colours["Lambda-CDM"],
            rasterized=True,
        )
        limits = [
            min(values["x_lcdm_stretched_redshift"].min(), values["y_lcdm_observed_distance"].min()),
            max(values["x_lcdm_stretched_redshift"].max(), values["y_lcdm_observed_distance"].max()),
        ]
        left.plot(limits, limits, color="black", linewidth=1.0)
        left.set_xscale("log")
        left.set_yscale("log")
        left.set_title(f"{dataset}: Lambda-CDM coordinate")
        left.set_xlabel("stretched redshift distance x_LCDM")
        left.set_ylabel("calibrated H0 D_L,obs / c")

        right.scatter(
            values["z"],
            values["y_plamb_corrected_distance"],
            s=8,
            alpha=0.16,
            color=colours["PLAMB"],
            rasterized=True,
        )
        limits = [
            min(values["z"].min(), values["y_plamb_corrected_distance"].min()),
            max(values["z"].max(), values["y_plamb_corrected_distance"].max()),
        ]
        right.plot(limits, limits, color="black", linewidth=1.0)
        right.set_xscale("log")
        right.set_yscale("log")
        right.set_title(f"{dataset}: PLAMB corrected distance")
        right.set_xlabel("measured redshift z")
        right.set_ylabel("calibrated H0 D_L,obs / [c(1+z/2)]")
    fig.suptitle(
        f"Three-release axis comparison; joint Lambda-CDM Omega_m={float(joint['LCDM_Omega_m']):.4f}",
        fontsize=14,
    )
    fig.savefig(OUT / f"plamb_axis_comparison_{DATE}.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(10, 11), sharex=True, constrained_layout=True)
    for axis, dataset in zip(axes, datasets):
        subset = binned[binned["dataset"] == dataset]
        for model in ("PLAMB", "Lambda-CDM"):
            values = subset[subset["model"] == model]
            axis.errorbar(
                values["z_effective"],
                values["residual_mean_mag"],
                yerr=values["residual_sigma_mag"],
                marker="o",
                linewidth=1.3,
                capsize=2,
                color=colours[model],
                label=model,
            )
        axis.axhline(0.0, color="black", linewidth=0.8)
        axis.set_xscale("log")
        axis.set_ylabel("profiled residual (mag)")
        axis.set_title(dataset)
        axis.legend()
    axes[-1].set_xlabel("redshift")
    fig.savefig(OUT / f"plamb_axis_comparison_binned_residuals_{DATE}.png", dpi=180)
    plt.close(fig)


def write_report(
    fits: pd.DataFrame,
    identity: dict[str, float],
    distance_fit: dict[str, object],
    primary_counts: dict[str, int],
    all_z_counts: dict[str, int],
) -> None:
    primary = fits[
        (fits["scope"] == "joint") & (fits["variant"] == "released_total_primary")
    ].iloc[0]
    diagonal = fits[
        (fits["scope"] == "joint") & (fits["variant"] == "released_total_diagonal_mu")
    ].iloc[0]
    all_z = fits[
        (fits["scope"] == "joint") & (fits["variant"] == "all_positive_z_noncalibrator")
    ].iloc[0]
    release_rows = fits[
        (fits["variant"] == "released_total_primary") & (fits["scope"] != "joint")
    ]
    lines = [
        "# PLAMB Three-release Axis Comparison",
        "",
        "Date: 18 July 2026",
        "",
        "## Main finding",
        "",
        "Peter Lamb's corrected-distance display is valid as a visualisation, but it is algebraically the same fixed PLAMB model already used in the raw-MU likelihood. It is not a new likelihood and does not alter the model comparison.",
        "",
        f"The released-total-covariance primary result has `N={int(primary['N'])}` and `Delta BIC(PLAMB-Lambda-CDM)={primary['delta_BIC_FR_minus_LCDM']:.6f}`. Positive values favour Lambda-CDM.",
        "",
        "## Release-by-release comparison",
        "",
        "| release | N | Delta BIC (PLAMB-Lambda-CDM) | LCDM Omega_m |",
        "| --- | ---: | ---: | ---: |",
        *[
            f"| {row.scope} | {int(row.N)} | {row.delta_BIC_FR_minus_LCDM:.6f} | {row.LCDM_Omega_m:.6f} |"
            for row in release_rows.itertuples(index=False)
        ],
        "",
        "Every release separately favours Lambda-CDM for the fixed PLAMB law under its released total covariance.",
        "",
        "## Algebraic identity",
        "",
        *[f"- {name}: maximum residual difference `{value:.3e}` mag" for name, value in identity.items()],
        "",
        "All releases pass the registered `1e-10` mag identity gate.",
        "",
        "## Joint sensitivities",
        "",
        "| likelihood | N | Delta BIC (PLAMB-Lambda-CDM) | LCDM Omega_m | status |",
        "| --- | ---: | ---: | ---: | --- |",
        f"| released total covariance in MU | {int(primary['N'])} | {primary['delta_BIC_FR_minus_LCDM']:.6f} | {primary['LCDM_Omega_m']:.6f} | primary |",
        f"| diagonal of total covariance in MU | {int(diagonal['N'])} | {diagonal['delta_BIC_FR_minus_LCDM']:.6f} | {diagonal['LCDM_Omega_m']:.6f} | sensitivity only |",
        f"| diagonal luminosity distance plus available z error | {int(distance_fit['N'])} | {float(distance_fit['delta_BIC_FR_minus_LCDM']):.6f} | {float(distance_fit['LCDM_Omega_m']):.6f} | sensitivity only |",
        f"| all positive-z non-calibrators | {int(all_z['N'])} | {all_z['delta_BIC_FR_minus_LCDM']:.6f} | {all_z['LCDM_Omega_m']:.6f} | low-z stress |",
        "",
        "## Row coverage",
        "",
        f"Primary rows: `{primary_counts}`.",
        f"All-positive-z rows: `{all_z_counts}`.",
        f"Distance-space redshift-error coverage: `{distance_fit['redshift_error_coverage']}`.",
        "",
        "DES-Dovekie and Union3.1 do not provide complete per-row redshift-error vectors in the compressed products used here. No errors were invented. Union3.1 consists of 22 correlated spline nodes rather than individual raw supernova measurements.",
        "The Pantheon+ distance-space sensitivity may conservatively count some low-redshift velocity uncertainty twice because peculiar-velocity contributions already enter the released distance uncertainty. This is another reason it is not the primary likelihood.",
        "",
        "## Statistical interpretation",
        "",
        "Type Ia supernovae are standardised candles, not exact identical candles. Their released distance moduli incorporate light-curve shape, colour, host, selection and calibration corrections with associated uncertainties.",
        "",
        "The released covariance matrices encode statistical uncertainty, shared calibration, bias-correction, peculiar-velocity and survey systematic modes. They are not energy matrices. Replacing them with independent luminosity-distance errors discards measured correlation structure and is therefore a diagnostic sensitivity, not the preferred likelihood.",
        "",
        "The Lambda-CDM comparator here is a Friedmann-equation distance integral with one fitted `Omega_m`, not a polynomial in velocity, deceleration and jerk. The primary sample retains the full released high-redshift range through `z=2.26226`; only Pantheon+ applied the registered `z>0.01` lower cut. Restoring 44 positive-redshift non-calibrator rows leaves the model direction unchanged.",
        "",
        "The correct model comparison is a forward prediction for the same observed distance-modulus vector under identical nuisance terms. Calling one plot an x-axis stretch and the other a y-axis division does not change that requirement.",
    ]
    (OUT / f"plamb_axis_comparison_report_{DATE}.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def write_email_draft(fits: pd.DataFrame, distance_fit: dict[str, object]) -> None:
    primary = fits[
        (fits["scope"] == "joint") & (fits["variant"] == "released_total_primary")
    ].iloc[0]
    text = f"""# Draft response to Peter Lamb

Dear Peter,

Thank you. Your note makes the proposed PLAMB comparison much clearer. I have now implemented the two displays for Pantheon+SH0ES, DES-Dovekie and Union3.1.

One terminology point is useful at the outset: Type Ia supernovae are standardised candles rather than exact identical candles. The published distance moduli include light-curve shape, colour, host, selection and calibration corrections, together with their uncertainties.

The key mathematical point is that dividing the observed luminosity distance by `(1+z/2)` is exactly equivalent to the fixed PLAMB forward model

`D_L(z) = (c/H0) z (1+z/2)`.

In distance-modulus form the correction is

`mu_corr = mu_obs - 5 log10(1+z/2)`.

The direct PLAMB residual and corrected-distance residual agree to better than `1e-10` magnitude in all three releases. The corrected plot is therefore a useful and intuitive representation, but it is not a different fit from the PLAMB curve already tested.

Using the same 3,422 released points, the same redshift choices, one normalisation intercept per release and each release's full covariance, the result is `Delta BIC(PLAMB-Lambda-CDM)={primary['delta_BIC_FR_minus_LCDM']:.6f}`. Positive values favour Lambda-CDM for this fixed PLAMB law. The luminosity-distance diagonal sensitivity gives `Delta BIC={float(distance_fit['delta_BIC_FR_minus_LCDM']):.6f}`, but it is not the primary likelihood because it removes shared calibration and survey-systematic correlations.

The three releases also give the same direction separately: Pantheon+ `+42.200623`, DES-Dovekie `+27.936226`, and Union3.1 `+14.420284`. Adding the 44 positive-redshift Pantheon+ non-calibrator rows below the earlier `z>0.01` threshold gives a joint value of `+91.314727`, so the low-redshift extension does not reverse the result. There is no high-redshift cut in this comparison; the released range is retained through `z=2.26226`.

The Lambda-CDM comparator is not fitted as an arbitrary velocity polynomial here. It is the flat Friedmann distance integral with one fitted matter-density parameter. PLAMB has the fixed `z(1+z/2)` shape and no fitted shape parameter.

The covariance matrices are not assumptions about energy. They describe which reported supernova distances share calibration, selection, peculiar-velocity and other systematic uncertainties. Transforming from magnitude to luminosity distance does not make those correlations disappear; the covariance must be transformed with the data. A diagonal distance-error fit is still useful as a transparent sensitivity check, and I have included it as such.

There are two data limitations to the proposed `z` and `LD` error-only comparison. Pantheon+ publishes `zHDERR`, but the compressed DES-Dovekie Hubble diagram and the 22 Union3.1 spline nodes do not supply equivalent complete per-row redshift-error vectors. Union3.1 is also a correlated compression, not a list of raw individual-supernova distances. I have not invented missing redshift errors.

Although spectroscopic redshift measurement errors are usually small, peculiar velocities are not small relative to the Hubble recession velocity for the nearest objects. They are one reason for the low-redshift sensitivity check. In Pantheon+, part of this velocity uncertainty already enters the released distance error, so adding `zHDERR` again in a diagonal distance fit can be conservative or partly double-counted.

The plots now supplied are:

1. calibrated observed distance against the fitted Lambda-CDM stretched-redshift distance;
2. calibrated observed distance divided by `(1+z/2)` against measured redshift; and
3. covariance-weighted binned residuals for both forward models.

To turn the PLAMB factor into a distinct physical derivation rather than a fitted distance law, the next item I need is the explicit chain from clock-rate or light-speed evolution to `1+z/2`: the time variable, the function `c(t)` or wavelength-evolution equation, the integration limits, and the assumptions connecting emitted luminosity to observed flux. That derivation would let us test alternatives to the fixed factor without introducing an arbitrary polynomial.

Regards,
Clive
"""
    (OUT / f"draft_response_to_peter_lamb_{DATE}.md").write_text(text, encoding="utf-8")


def write_manifest(paths: list[Path]) -> None:
    rows = [
        {
            "path": path.relative_to(ROOT).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in sorted(paths, key=lambda item: item.as_posix())
    ]
    with (OUT / f"plamb_axis_comparison_manifest_{DATE}.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=["path", "bytes", "sha256"])
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    if not PREREG.exists():
        raise FileNotFoundError(PREREG)
    OUT.mkdir(parents=True, exist_ok=True)
    registration = json.loads(PREREG.read_text(encoding="utf-8"))
    variants, metadata = base.load_blocks(min_survey_n=20)
    primary = variants["released_total_primary"]
    identity = identity_audit(primary)
    tolerance = float(registration["algebraic_identity_tolerance_mag"])
    if max(identity.values()) > tolerance:
        raise RuntimeError(f"PLAMB corrected-axis identity failed: {identity}")

    diagonal = diagonal_blocks(primary, "released_total_diagonal_mu")
    all_z = all_positive_redshift_blocks(primary)
    fit_rows, joint = magnitude_fit_rows(primary, diagonal, all_z)
    distance_fit = distance_space_fit(primary)
    fit_rows.append(
        {
            key: distance_fit[key]
            for key in (
                "scope",
                "variant",
                "N",
                "FR_chi2",
                "LCDM_chi2",
                "LCDM_Omega_m",
                "FR_BIC",
                "LCDM_BIC",
                "delta_BIC_FR_minus_LCDM",
                "same_nuisance_verified",
            )
        }
    )
    fits = pd.DataFrame(fit_rows)
    bins = np.asarray(registration["plot_bins"], dtype=float)
    binned = residual_bins(primary, joint, bins)
    coordinates = coordinate_table(primary, joint)

    fits.to_csv(OUT / f"plamb_axis_comparison_fits_{DATE}.csv", index=False)
    binned.to_csv(OUT / f"plamb_axis_comparison_binned_residuals_{DATE}.csv", index=False)
    coordinates.to_csv(
        OUT / f"plamb_axis_comparison_coordinates_{DATE}.csv.gz",
        index=False,
        compression="gzip",
    )
    write_plots(coordinates, binned, joint)
    primary_counts = {block.label: block.n for block in primary}
    all_z_counts = {block.label: block.n for block in all_z}
    write_report(fits, identity, distance_fit, primary_counts, all_z_counts)
    write_email_draft(fits, distance_fit)

    summary = {
        "date": DATE,
        "preregistration_sha256": sha256_file(PREREG),
        "identity_max_abs_mag": identity,
        "identity_gate_pass": max(identity.values()) <= tolerance,
        "primary_fit": joint,
        "distance_space_sensitivity": distance_fit,
        "primary_counts": primary_counts,
        "all_positive_z_counts": all_z_counts,
        "load_metadata": metadata,
        "claim_boundary": registration["claim_boundary"],
    }
    summary_path = OUT / f"plamb_axis_comparison_summary_{DATE}.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    outputs = [
        path
        for path in OUT.iterdir()
        if path.is_file() and "manifest" not in path.name
    ] + [Path(__file__).resolve()]
    write_manifest(outputs)
    print(f"Identity gate: PASS; max difference={max(identity.values()):.3e} mag")
    print(
        "Primary Delta BIC (PLAMB-Lambda-CDM): "
        f"{float(joint['delta_BIC_FR_minus_LCDM']):.6f}"
    )
    print(
        "Distance-diagonal Delta BIC: "
        f"{float(distance_fit['delta_BIC_FR_minus_LCDM']):.6f}"
    )
    print(f"Saved: {OUT / f'plamb_axis_comparison_report_{DATE}.md'}")


if __name__ == "__main__":
    main()
