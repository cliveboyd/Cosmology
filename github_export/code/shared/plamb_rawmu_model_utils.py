r"""Shared raw-MU distance-law helpers for PLAMB diagnostics."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

try:
    from scipy.optimize import minimize_scalar

    HAVE_SCIPY = True
except Exception:  # pragma: no cover
    minimize_scalar = None
    HAVE_SCIPY = False

from diagnose_pantheon_rawmu_fr import C_KMS
from diagnose_pantheon_survey_id_projection import solve_offsets


@dataclass(frozen=True)
class BaselineSpec:
    name: str
    family: str
    param_name: str | None
    fixed_value: float | None
    default_value: float
    bounds: tuple[float, float]
    k_shape: int


BASELINE_SPECS: list[BaselineSpec] = [
    BaselineSpec("PLAMB_BETA05_H0fixed675", "plamb", "beta", 0.5, 0.5, (-0.25, 1.25), 0),
    BaselineSpec("PLAMB_BETAfree_H0fixed675", "plamb", "beta", None, 0.5, (-0.25, 1.25), 1),
    BaselineSpec("LCDM_Om03_H0fixed675", "lcdm", "Om", 0.3, 0.3, (0.05, 0.6), 0),
    BaselineSpec("LCDM_Omfree_H0fixed675", "lcdm", "Om", None, 0.3, (0.05, 0.6), 1),
    BaselineSpec("COSMO_q0free_H0fixed675", "cosmographic", "q0", None, -0.55, (-2.0, 1.0), 1),
]


def specs_from_request(requested: str) -> list[BaselineSpec]:
    if requested.lower() == "all":
        return BASELINE_SPECS
    wanted = {item.strip() for item in requested.split(",") if item.strip()}
    specs = [spec for spec in BASELINE_SPECS if spec.name in wanted or spec.family in wanted]
    missing = wanted - {spec.name for spec in specs} - {spec.family for spec in specs}
    if missing:
        raise ValueError(f"Unknown baseline model(s): {', '.join(sorted(missing))}")
    return specs


def lcdm_integral(z: np.ndarray, om: float, n_grid: int = 4096) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    if om <= 0.0 or om >= 1.5 or not np.all(np.isfinite(z)):
        return np.full_like(z, np.inf, dtype=float)
    zmax = float(np.nanmax(z)) if z.size else 0.0
    if zmax <= 0.0:
        return np.zeros_like(z)
    grid = np.linspace(0.0, zmax, max(128, int(n_grid)))
    ez = np.sqrt(om * np.power(1.0 + grid, 3.0) + (1.0 - om))
    if np.any(ez <= 0.0) or not np.all(np.isfinite(ez)):
        return np.full_like(z, np.inf, dtype=float)
    inv_e = 1.0 / ez
    step = grid[1] - grid[0]
    cumulative = np.empty_like(grid)
    cumulative[0] = 0.0
    cumulative[1:] = np.cumsum(0.5 * (inv_e[1:] + inv_e[:-1]) * step)
    return np.interp(z, grid, cumulative)


def mu_for_baseline(z: np.ndarray, spec: BaselineSpec, h0: float, value: float) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    if spec.family == "plamb":
        d = (C_KMS / h0) * z * (1.0 + value * z)
    elif spec.family == "lcdm":
        d = (C_KMS / h0) * (1.0 + z) * lcdm_integral(z, value)
    elif spec.family == "cosmographic":
        # Second-order luminosity-distance expansion:
        # DL = c/H0 * z * [1 + (1-q0)z/2].
        d = (C_KMS / h0) * z * (1.0 + 0.5 * (1.0 - value) * z)
    else:
        raise ValueError(f"Unknown baseline family: {spec.family}")
    if np.any(d <= 0.0) or not np.all(np.isfinite(d)):
        return np.full_like(z, np.inf, dtype=float)
    return 5.0 * np.log10(d) + 25.0


def one_column(n: int) -> np.ndarray:
    return np.ones((n, 1), dtype=float)


def fit_global_offset_for_value(
    z: np.ndarray,
    mu: np.ndarray,
    precision: np.ndarray,
    spec: BaselineSpec,
    h0: float,
    value: float,
) -> tuple[float, float, np.ndarray]:
    model = mu_for_baseline(z, spec, h0, value)
    if not np.all(np.isfinite(model)):
        return float("nan"), float("inf"), np.full_like(mu, np.nan)
    offsets, profiled, chi2 = solve_offsets(mu - model, precision, one_column(len(z)))
    return float(offsets[0]), float(chi2), profiled


def score_with_offset_for_value(
    z: np.ndarray,
    mu: np.ndarray,
    precision: np.ndarray,
    spec: BaselineSpec,
    h0: float,
    value: float,
    offset: float,
) -> tuple[float, float, np.ndarray]:
    model = mu_for_baseline(z, spec, h0, value)
    if not np.all(np.isfinite(model)):
        return float("inf"), float("nan"), np.full_like(mu, np.nan)
    residual = mu - model - offset
    chi2 = float(residual @ precision @ residual)
    rms = float(np.sqrt(np.nanmean(np.square(residual)))) if residual.size else float("nan")
    return chi2, rms, residual


def fit_shape_with_global_offset(
    z: np.ndarray,
    mu: np.ndarray,
    precision: np.ndarray,
    spec: BaselineSpec,
    h0: float,
    grid_steps: int = 401,
) -> dict[str, object]:
    if spec.fixed_value is not None:
        value = spec.fixed_value
        offset, chi2, profiled = fit_global_offset_for_value(z, mu, precision, spec, h0, value)
        return {
            "model": spec.name,
            "family": spec.family,
            "param_name": spec.param_name,
            "param_value": value,
            "h0": h0,
            "global_offset_mag": offset,
            "chi2": chi2,
            "dof": max(len(z) - 1, 1),
            "profiled_rms_mag": float(np.sqrt(np.nanmean(np.square(profiled)))) if profiled.size else float("nan"),
            "method": "fixed_shape_global_offset",
            "success": np.isfinite(chi2),
        }

    def objective(value: float) -> float:
        _offset, chi2, _profiled = fit_global_offset_for_value(z, mu, precision, spec, h0, float(value))
        return chi2

    if HAVE_SCIPY:
        opt = minimize_scalar(objective, bounds=spec.bounds, method="bounded")
        value = float(opt.x)
        success = bool(opt.success)
        method = f"minimize_{spec.param_name}_global_offset"
    else:
        grid = np.linspace(spec.bounds[0], spec.bounds[1], grid_steps)
        vals = np.array([objective(v) for v in grid], dtype=float)
        value = float(grid[int(np.nanargmin(vals))])
        success = bool(np.isfinite(vals).any())
        method = f"grid_{spec.param_name}_global_offset"

    offset, chi2, profiled = fit_global_offset_for_value(z, mu, precision, spec, h0, value)
    return {
        "model": spec.name,
        "family": spec.family,
        "param_name": spec.param_name,
        "param_value": value,
        "h0": h0,
        "global_offset_mag": offset,
        "chi2": chi2,
        "dof": max(len(z) - 2, 1),
        "profiled_rms_mag": float(np.sqrt(np.nanmean(np.square(profiled)))) if profiled.size else float("nan"),
        "method": method,
        "success": bool(success and np.isfinite(chi2)),
    }
