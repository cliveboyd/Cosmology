r"""
Hilbert-Huang resonance diagnostic for the SU2/chirality investigation.

This standalone diagnostic is inspired by arXiv:2605.20333. It does not alter
the SN/BAO/Planck likelihood. Instead, it asks whether a nonlinear
time-frequency method can recover the resonance chain

    SU2 condensate Q(t) -> pseudoscalar chi_dot(t) -> U(1) helicity modes.

The default "synthetic" mode builds a controlled toy signal with the same
paper-level structure:

    * Q(t) is an oscillating isotropic SU2 condensate proxy.
    * Q^2 dQ/dt sources chi_dot(t), with first and third harmonics.
    * chi_dot(t) modulates U(1) helicity-mode frequencies.
    * The U(1) modes obey C_lambda'' + [k^2 - lambda k beta(t)] C_lambda = 0.

The script then applies a small built-in EMD implementation and Hilbert
transform to estimate instantaneous phase/frequency and coupling diagnostics:

    * Q -> chi phase locking.
    * chi -> U(1) 2:1 resonance phase locking.
    * U(1) energy growth and helicity asymmetry.
    * A compact resonance-coupling score.

An optional CSV mode can run the same HHT decomposition on an arbitrary series,
for example future model residuals ordered by log(1+z), lookback time, or
comoving distance.

Outputs:
    plamb_runs/diagnostics/hht_resonance/hht_resonance_synthetic_signals.csv
    plamb_runs/diagnostics/hht_resonance/hht_resonance_imf_summary.csv
    plamb_runs/diagnostics/hht_resonance/hht_resonance_coupling_summary.csv
    plamb_runs/diagnostics/hht_resonance/hht_resonance_summary.md
    plamb_runs/diagnostics/hht_resonance/hht_resonance_overview.png
    plamb_runs/diagnostics/hht_resonance/hht_resonance_hilbert_spectrum.png
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

try:
    from scipy.interpolate import PchipInterpolator
    from scipy.signal import hilbert as scipy_hilbert

    HAVE_SCIPY = True
except Exception:  # pragma: no cover - fallback for lean Python installs.
    PchipInterpolator = None
    scipy_hilbert = None
    HAVE_SCIPY = False


ROOT = Path(__file__).resolve().parent
OUTDIR = ROOT / "plamb_runs" / "diagnostics" / "hht_resonance"


@dataclass(frozen=True)
class SyntheticConfig:
    omega_q: float
    m_chi: float
    k_mode: float
    su2_chi_coupling: float
    u1_coupling: float
    bias: float
    damping: float
    duration_periods: float
    samples: int
    q0: float
    anharmonic: float
    noise: float
    seed: int
    max_imfs: int
    sift_max_iter: int
    sift_sd: float
    top_n: int
    no_plots: bool


@dataclass(frozen=True)
class CsvConfig:
    input_csv: Path
    x_col: str
    y_col: str
    ref_col: str | None
    label: str
    samples: int
    max_imfs: int
    sift_max_iter: int
    sift_sd: float
    no_plots: bool


def ensure_outdir() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)


def clean_signal(y: np.ndarray) -> np.ndarray:
    arr = np.asarray(y, dtype=float)
    if arr.ndim != 1:
        raise ValueError("expected one-dimensional signal")
    finite = np.isfinite(arr)
    if not finite.all():
        if finite.sum() < 2:
            raise ValueError("signal contains fewer than two finite values")
        x = np.arange(arr.size)
        arr = np.interp(x, x[finite], arr[finite])
    arr = arr - np.nanmean(arr)
    return arr


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def safe_slug(value: str) -> str:
    slug = "".join(ch if ch.isalnum() or ch in {"_", "-"} else "_" for ch in value.strip())
    slug = "_".join(part for part in slug.split("_") if part)
    return slug or "series"


def read_csv_columns(path: Path) -> dict[str, list[str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError(f"{path} has no header")
        cols = {name: [] for name in reader.fieldnames}
        for row in reader:
            for name in cols:
                cols[name].append(row.get(name, ""))
    return cols


def as_float_array(values: Iterable[object], label: str) -> np.ndarray:
    out: list[float] = []
    for item in values:
        try:
            out.append(float(item))
        except Exception as exc:
            raise ValueError(f"column {label!r} contains nonnumeric value {item!r}") from exc
    arr = np.asarray(out, dtype=float)
    if arr.size < 8:
        raise ValueError(f"column {label!r} needs at least 8 points for HHT")
    return arr


def extrema_indices(y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return local maxima and minima indices for a one-dimensional series."""
    dy = np.diff(y)
    if dy.size == 0:
        return np.array([], dtype=int), np.array([], dtype=int)

    # Replace zero slopes with nearest nonzero sign so plateaus do not erase extrema.
    signs = np.sign(dy)
    for i in range(1, signs.size):
        if signs[i] == 0:
            signs[i] = signs[i - 1]
    for i in range(signs.size - 2, -1, -1):
        if signs[i] == 0:
            signs[i] = signs[i + 1]

    change = np.diff(signs)
    maxima = np.where(change < 0)[0] + 1
    minima = np.where(change > 0)[0] + 1
    return maxima.astype(int), minima.astype(int)


def envelope_from_extrema(x: np.ndarray, y: np.ndarray, idx: np.ndarray) -> np.ndarray:
    """Interpolate an envelope through extrema, forcing endpoints for stability."""
    if idx.size == 0:
        return np.full_like(y, float(np.mean(y)))

    points = np.unique(np.concatenate(([0], idx, [len(y) - 1]))).astype(int)
    xp = x[points]
    yp = y[points]

    if points.size >= 4 and HAVE_SCIPY and PchipInterpolator is not None:
        return np.asarray(PchipInterpolator(xp, yp, extrapolate=True)(x), dtype=float)
    return np.interp(x, xp, yp)


def zero_crossing_count(y: np.ndarray) -> int:
    signs = np.signbit(y)
    return int(np.count_nonzero(signs[1:] != signs[:-1]))


def is_imf_like(y: np.ndarray) -> bool:
    maxima, minima = extrema_indices(y)
    extrema_count = maxima.size + minima.size
    zc = zero_crossing_count(y)
    return abs(extrema_count - zc) <= 1


def emd(
    y: np.ndarray,
    max_imfs: int = 6,
    sift_max_iter: int = 80,
    sift_sd: float = 0.08,
) -> tuple[list[np.ndarray], np.ndarray]:
    """Small empirical-mode decomposition implementation.

    It is intentionally conservative and dependency-light. For publication
    work, cross-check candidate features with EEMD/CEEMDAN from a specialist
    package, but this is sufficient for a first project diagnostic.
    """
    residue = clean_signal(y)
    x = np.arange(residue.size, dtype=float)
    imfs: list[np.ndarray] = []

    for _ in range(max_imfs):
        maxima, minima = extrema_indices(residue)
        if maxima.size + minima.size < 4:
            break

        h = residue.copy()
        for _sift in range(sift_max_iter):
            maxima, minima = extrema_indices(h)
            if maxima.size < 2 or minima.size < 2:
                break

            upper = envelope_from_extrema(x, h, maxima)
            lower = envelope_from_extrema(x, h, minima)
            mean_env = 0.5 * (upper + lower)
            h_next = h - mean_env

            denom = float(np.sum(h * h)) + 1e-30
            sd = float(np.sum((h - h_next) ** 2) / denom)
            h = h_next

            if sd < sift_sd and is_imf_like(h):
                break

        imfs.append(h)
        residue = residue - h

        maxima, minima = extrema_indices(residue)
        if maxima.size + minima.size < 4:
            break
        if np.std(residue) < 1e-8 * (np.std(y) + 1e-30):
            break

    return imfs, residue


def analytic_signal(y: np.ndarray) -> np.ndarray:
    if HAVE_SCIPY and scipy_hilbert is not None:
        return np.asarray(scipy_hilbert(y), dtype=complex)

    # FFT-based Hilbert transform fallback following the standard analytic
    # signal construction.
    n = y.size
    spectrum = np.fft.fft(y)
    h = np.zeros(n)
    if n % 2 == 0:
        h[0] = 1
        h[n // 2] = 1
        h[1 : n // 2] = 2
    else:
        h[0] = 1
        h[1 : (n + 1) // 2] = 2
    return np.fft.ifft(spectrum * h)


def hilbert_features(t: np.ndarray, y: np.ndarray) -> dict[str, np.ndarray]:
    analytic = analytic_signal(clean_signal(y))
    amp = np.abs(analytic)
    phase = np.unwrap(np.angle(analytic))
    dt = float(np.median(np.diff(t)))
    freq = np.gradient(phase, dt) / (2.0 * math.pi)
    return {
        "amplitude": amp,
        "phase": phase,
        "frequency": freq,
    }


def median_positive_frequency(freq: np.ndarray, amp: np.ndarray) -> float:
    mask = np.isfinite(freq) & np.isfinite(amp) & (freq > 0) & (amp > np.percentile(amp, 20))
    if mask.sum() < 3:
        return math.nan
    return float(np.median(freq[mask]))


def weighted_frequency(freq: np.ndarray, amp: np.ndarray) -> float:
    mask = np.isfinite(freq) & np.isfinite(amp) & (freq > 0)
    if mask.sum() < 3:
        return math.nan
    weights = amp[mask] ** 2
    total = float(np.sum(weights))
    if total <= 0:
        return math.nan
    return float(np.sum(freq[mask] * weights) / total)


def summarize_imfs(label: str, t: np.ndarray, signal: np.ndarray, imfs: list[np.ndarray], residue: np.ndarray) -> list[dict[str, object]]:
    total_energy = float(np.sum(clean_signal(signal) ** 2)) + 1e-30
    rows: list[dict[str, object]] = []
    for i, imf in enumerate(imfs, start=1):
        features = hilbert_features(t, imf)
        amp = features["amplitude"]
        freq = features["frequency"]
        rows.append(
            {
                "series": label,
                "component": f"IMF{i}",
                "component_index": i,
                "energy": float(np.sum(imf * imf)),
                "energy_fraction": float(np.sum(imf * imf) / total_energy),
                "median_amplitude": float(np.median(amp)),
                "max_amplitude": float(np.max(amp)),
                "median_positive_frequency": median_positive_frequency(freq, amp),
                "amplitude_weighted_frequency": weighted_frequency(freq, amp),
                "zero_crossings": zero_crossing_count(imf),
                "is_imf_like": is_imf_like(imf),
            }
        )
    rows.append(
        {
            "series": label,
            "component": "residue",
            "component_index": len(imfs) + 1,
            "energy": float(np.sum(residue * residue)),
            "energy_fraction": float(np.sum(residue * residue) / total_energy),
            "median_amplitude": math.nan,
            "max_amplitude": math.nan,
            "median_positive_frequency": math.nan,
            "amplitude_weighted_frequency": math.nan,
            "zero_crossings": zero_crossing_count(residue),
            "is_imf_like": False,
        }
    )
    return rows


def pick_imf_near(
    label: str,
    summaries: list[dict[str, object]],
    imfs_by_label: dict[str, list[np.ndarray]],
    target_freq: float,
) -> tuple[int, np.ndarray | None, dict[str, object] | None]:
    candidates = [
        row
        for row in summaries
        if row["series"] == label
        and str(row["component"]).startswith("IMF")
        and np.isfinite(float(row["amplitude_weighted_frequency"]))
    ]
    if not candidates:
        return -1, None, None

    def key(row: dict[str, object]) -> tuple[float, float]:
        freq = float(row["amplitude_weighted_frequency"])
        energy = float(row["energy_fraction"])
        if target_freq > 0 and freq > 0:
            offset = abs(math.log(freq / target_freq))
        else:
            offset = math.inf
        return offset, -energy

    best = min(candidates, key=key)
    idx = int(best["component_index"]) - 1
    imfs = imfs_by_label.get(label, [])
    if idx < 0 or idx >= len(imfs):
        return -1, None, best
    return idx + 1, imfs[idx], best


def phase_locking_value(phase_a: np.ndarray, phase_b: np.ndarray, m: int = 1, n: int = 1) -> float:
    mask = np.isfinite(phase_a) & np.isfinite(phase_b)
    if mask.sum() < 8:
        return math.nan
    z = np.exp(1j * (m * phase_a[mask] - n * phase_b[mask]))
    return float(abs(np.mean(z)))


def pearson_corr(a: np.ndarray, b: np.ndarray) -> float:
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 3:
        return math.nan
    aa = a[mask] - np.mean(a[mask])
    bb = b[mask] - np.mean(b[mask])
    denom = math.sqrt(float(np.sum(aa * aa) * np.sum(bb * bb))) + 1e-30
    return float(np.sum(aa * bb) / denom)


def growth_ratio(energy: np.ndarray) -> float:
    n = energy.size
    if n < 10:
        return math.nan
    early = float(np.median(energy[: max(5, n // 20)]))
    late = float(np.median(energy[-max(5, n // 20) :]))
    return late / max(early, 1e-300)


def finite_log_score(value: float, scale: float = 10.0) -> float:
    if not np.isfinite(value) or value <= 0:
        return 0.0
    return float(max(0.0, min(1.0, math.log10(value) / math.log10(scale))))


def resonance_score(
    q_chi_plv: float,
    chi_u1_plv: float,
    u1_growth: float,
    resonance_mismatch: float,
) -> float:
    plv1 = 0.0 if not np.isfinite(q_chi_plv) else max(0.0, min(1.0, q_chi_plv))
    plv2 = 0.0 if not np.isfinite(chi_u1_plv) else max(0.0, min(1.0, chi_u1_plv))
    growth = finite_log_score(u1_growth, scale=100.0)
    match = 0.0 if not np.isfinite(resonance_mismatch) else math.exp(-abs(resonance_mismatch))
    return float(0.30 * plv1 + 0.35 * plv2 + 0.25 * growth + 0.10 * match)


def simulate_u1_mode(t: np.ndarray, beta: np.ndarray, k_mode: float, helicity: int, damping: float) -> tuple[np.ndarray, np.ndarray]:
    """Integrate C'' + damping C' + [k^2 - helicity k beta(t)] C = 0."""
    if t.size != beta.size:
        raise ValueError("time and beta arrays must have same length")
    c = np.zeros_like(t)
    v = np.zeros_like(t)
    c[0] = 1.0e-4
    v[0] = 0.0

    def accel(i: int, ci: float, vi: float) -> float:
        omega2 = k_mode * k_mode - helicity * k_mode * beta[i]
        return -damping * vi - omega2 * ci

    for i in range(t.size - 1):
        dt = float(t[i + 1] - t[i])
        a1 = accel(i, c[i], v[i])
        k1c = v[i]
        k1v = a1

        c2 = c[i] + 0.5 * dt * k1c
        v2 = v[i] + 0.5 * dt * k1v
        a2 = -damping * v2 - (k_mode * k_mode - helicity * k_mode * 0.5 * (beta[i] + beta[i + 1])) * c2
        k2c = v2
        k2v = a2

        c3 = c[i] + 0.5 * dt * k2c
        v3 = v[i] + 0.5 * dt * k2v
        a3 = -damping * v3 - (k_mode * k_mode - helicity * k_mode * 0.5 * (beta[i] + beta[i + 1])) * c3
        k3c = v3
        k3v = a3

        c4 = c[i] + dt * k3c
        v4 = v[i] + dt * k3v
        a4 = -damping * v4 - (k_mode * k_mode - helicity * k_mode * beta[i + 1]) * c4
        k4c = v4
        k4v = a4

        c[i + 1] = c[i] + dt * (k1c + 2 * k2c + 2 * k3c + k4c) / 6.0
        v[i + 1] = v[i] + dt * (k1v + 2 * k2v + 2 * k3v + k4v) / 6.0

    return c, v


def build_synthetic_signals(cfg: SyntheticConfig) -> dict[str, np.ndarray]:
    rng = np.random.default_rng(cfg.seed)
    period = 2.0 * math.pi / cfg.omega_q
    t_max = cfg.duration_periods * period
    t = np.linspace(0.0, t_max, cfg.samples)
    dt = float(np.median(np.diff(t)))

    q = cfg.q0 * (
        np.cos(cfg.omega_q * t)
        + cfg.anharmonic * np.cos(3.0 * cfg.omega_q * t + 0.35)
    )
    q_dot = np.gradient(q, dt)
    source = q * q * q_dot

    # Approximate the forced scalar response from first and third harmonics.
    # The damping term keeps the response finite near m_chi ~= n omega_q.
    denom1 = math.sqrt((cfg.m_chi * cfg.m_chi - cfg.omega_q * cfg.omega_q) ** 2 + (0.08 * cfg.omega_q) ** 2)
    denom3 = math.sqrt((cfg.m_chi * cfg.m_chi - (3.0 * cfg.omega_q) ** 2) ** 2 + (0.08 * 3.0 * cfg.omega_q) ** 2)
    amp1 = cfg.su2_chi_coupling * cfg.omega_q * cfg.q0**3 / max(denom1, 1e-12)
    amp3 = 0.35 * cfg.su2_chi_coupling * 3.0 * cfg.omega_q * cfg.q0**3 / max(denom3, 1e-12)
    chi_dot = (
        amp1 * np.cos(cfg.omega_q * t + 0.25)
        + amp3 * np.cos(3.0 * cfg.omega_q * t - 0.6)
        + cfg.bias
    )
    beta = cfg.u1_coupling * chi_dot

    c_plus, v_plus = simulate_u1_mode(t, beta, cfg.k_mode, helicity=+1, damping=cfg.damping)
    c_minus, v_minus = simulate_u1_mode(t, beta, cfg.k_mode, helicity=-1, damping=cfg.damping)
    energy_plus = 0.5 * (v_plus * v_plus + cfg.k_mode * cfg.k_mode * c_plus * c_plus)
    energy_minus = 0.5 * (v_minus * v_minus + cfg.k_mode * cfg.k_mode * c_minus * c_minus)
    helicity_asymmetry = (energy_plus - energy_minus) / np.maximum(energy_plus + energy_minus, 1e-300)

    if cfg.noise > 0:
        q = q + cfg.noise * np.std(q) * rng.normal(size=q.size)
        chi_dot = chi_dot + cfg.noise * np.std(chi_dot) * rng.normal(size=chi_dot.size)
        c_plus = c_plus + cfg.noise * np.std(c_plus) * rng.normal(size=c_plus.size)
        c_minus = c_minus + cfg.noise * np.std(c_minus) * rng.normal(size=c_minus.size)

    return {
        "t": t,
        "Q": q,
        "Q_dot": q_dot,
        "Q2_Q_dot_source": source,
        "chi_dot": chi_dot,
        "beta": beta,
        "U1_plus": c_plus,
        "U1_minus": c_minus,
        "U1_plus_velocity": v_plus,
        "U1_minus_velocity": v_minus,
        "U1_plus_energy": energy_plus,
        "U1_minus_energy": energy_minus,
        "helicity_asymmetry": helicity_asymmetry,
    }


def signal_rows(signals: dict[str, np.ndarray]) -> list[dict[str, object]]:
    keys = [k for k in signals if k != "t"]
    rows: list[dict[str, object]] = []
    for i, t_value in enumerate(signals["t"]):
        row: dict[str, object] = {"t": float(t_value)}
        for key in keys:
            row[key] = float(signals[key][i])
        rows.append(row)
    return rows


def decompose_series(
    t: np.ndarray,
    series: dict[str, np.ndarray],
    max_imfs: int,
    sift_max_iter: int,
    sift_sd: float,
) -> tuple[dict[str, list[np.ndarray]], dict[str, np.ndarray], list[dict[str, object]]]:
    imfs_by_label: dict[str, list[np.ndarray]] = {}
    residue_by_label: dict[str, np.ndarray] = {}
    summary_rows: list[dict[str, object]] = []

    for label, values in series.items():
        imfs, residue = emd(values, max_imfs=max_imfs, sift_max_iter=sift_max_iter, sift_sd=sift_sd)
        imfs_by_label[label] = imfs
        residue_by_label[label] = residue
        summary_rows.extend(summarize_imfs(label, t, values, imfs, residue))

    return imfs_by_label, residue_by_label, summary_rows


def synthetic_coupling_rows(
    cfg: SyntheticConfig,
    signals: dict[str, np.ndarray],
    imfs_by_label: dict[str, list[np.ndarray]],
    summary_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    t = signals["t"]
    f_q = cfg.omega_q / (2.0 * math.pi)
    f_u = cfg.k_mode / (2.0 * math.pi)

    q_idx, q_imf, q_row = pick_imf_near("Q", summary_rows, imfs_by_label, f_q)
    chi_idx, chi_imf, chi_row = pick_imf_near("chi_dot", summary_rows, imfs_by_label, f_q)
    plus_idx, plus_imf, plus_row = pick_imf_near("U1_plus", summary_rows, imfs_by_label, f_u)
    minus_idx, minus_imf, minus_row = pick_imf_near("U1_minus", summary_rows, imfs_by_label, f_u)

    q_chi_plv = math.nan
    q_chi_env_corr = math.nan
    chi_plus_plv = math.nan
    chi_minus_plv = math.nan
    chi_plus_env_corr = math.nan
    chi_minus_env_corr = math.nan

    if q_imf is not None and chi_imf is not None:
        q_features = hilbert_features(t, q_imf)
        chi_features = hilbert_features(t, chi_imf)
        q_chi_plv = phase_locking_value(q_features["phase"], chi_features["phase"], 1, 1)
        q_chi_env_corr = pearson_corr(q_features["amplitude"], chi_features["amplitude"])

    if chi_imf is not None and plus_imf is not None:
        chi_features = hilbert_features(t, chi_imf)
        plus_features = hilbert_features(t, plus_imf)
        chi_plus_plv = phase_locking_value(plus_features["phase"], chi_features["phase"], 2, 1)
        chi_plus_env_corr = pearson_corr(chi_features["amplitude"], plus_features["amplitude"])

    if chi_imf is not None and minus_imf is not None:
        chi_features = hilbert_features(t, chi_imf)
        minus_features = hilbert_features(t, minus_imf)
        chi_minus_plv = phase_locking_value(minus_features["phase"], chi_features["phase"], 2, 1)
        chi_minus_env_corr = pearson_corr(chi_features["amplitude"], minus_features["amplitude"])

    plus_growth = growth_ratio(signals["U1_plus_energy"])
    minus_growth = growth_ratio(signals["U1_minus_energy"])
    max_growth = max(plus_growth, minus_growth)
    final_asymmetry = float(np.median(signals["helicity_asymmetry"][-max(8, signals["t"].size // 20) :]))
    resonance_mismatch = (2.0 * cfg.k_mode - cfg.omega_q) / cfg.omega_q
    raw_chirality_flag = abs(final_asymmetry) > 0.05
    intrinsic_chirality_flag = raw_chirality_flag and abs(cfg.bias) > 1e-12
    max_chi_u1_plv = max(chi_plus_plv, chi_minus_plv)
    score = resonance_score(q_chi_plv, max_chi_u1_plv, max_growth, resonance_mismatch)

    row = {
        "mode": "synthetic",
        "omega_Q": cfg.omega_q,
        "expected_f_Q_cycles_per_t": f_q,
        "m_chi": cfg.m_chi,
        "k_mode": cfg.k_mode,
        "expected_f_U1_cycles_per_t": f_u,
        "resonance_condition": "2 k ~= omega_Q for first Hill/Mathieu band",
        "resonance_mismatch_fraction": resonance_mismatch,
        "Q_selected_imf": q_idx,
        "chi_dot_selected_imf": chi_idx,
        "U1_plus_selected_imf": plus_idx,
        "U1_minus_selected_imf": minus_idx,
        "Q_selected_frequency": None if q_row is None else q_row["amplitude_weighted_frequency"],
        "chi_dot_selected_frequency": None if chi_row is None else chi_row["amplitude_weighted_frequency"],
        "U1_plus_selected_frequency": None if plus_row is None else plus_row["amplitude_weighted_frequency"],
        "U1_minus_selected_frequency": None if minus_row is None else minus_row["amplitude_weighted_frequency"],
        "Q_chi_1to1_PLV": q_chi_plv,
        "Q_chi_envelope_corr": q_chi_env_corr,
        "chi_U1_plus_1to2_PLV": chi_plus_plv,
        "chi_U1_minus_1to2_PLV": chi_minus_plv,
        "max_chi_U1_1to2_PLV": max_chi_u1_plv,
        "chi_U1_plus_envelope_corr": chi_plus_env_corr,
        "chi_U1_minus_envelope_corr": chi_minus_env_corr,
        "U1_plus_energy_growth_late_over_early": plus_growth,
        "U1_minus_energy_growth_late_over_early": minus_growth,
        "max_U1_energy_growth": max_growth,
        "median_final_helicity_asymmetry": final_asymmetry,
        "raw_final_helicity_asymmetry_abs_gt_0p05": raw_chirality_flag,
        "intrinsic_chirality_flag_requires_bias": intrinsic_chirality_flag,
        "bias": cfg.bias,
        "damping": cfg.damping,
        "resonance_coupling_score_0to1": score,
        "interpretation": interpret_score(score, intrinsic_chirality_flag, raw_chirality_flag, cfg.bias),
    }
    return [row]


def interpret_score(score: float, intrinsic_chirality_flag: bool, raw_chirality_flag: bool, bias: float) -> str:
    if score >= 0.70:
        strength = "strong HHT recovery of the synthetic resonance chain"
    elif score >= 0.45:
        strength = "moderate HHT recovery of the synthetic resonance chain"
    else:
        strength = "weak HHT recovery; inspect parameters or signal quality"

    if intrinsic_chirality_flag:
        chirality = "with bias-supported helicity imbalance"
    elif abs(bias) < 1e-12:
        if raw_chirality_flag:
            chirality = "with raw finite-phase asymmetry, but no intrinsic chirality in the zero-bias limit"
        else:
            chirality = "with no intrinsic chirality in the zero-bias limit"
    else:
        chirality = "without a resolved helicity imbalance at this setting"
    return f"{strength}, {chirality}"


def fmt_float(value: object, precision: int = 4) -> str:
    try:
        number = float(value)
    except Exception:
        return str(value)
    if not np.isfinite(number):
        return "nan"
    return f"{number:.{precision}g}"


def write_synthetic_summary(
    path: Path,
    cfg: SyntheticConfig,
    coupling_rows: list[dict[str, object]],
    imf_rows: list[dict[str, object]],
    plot_paths: list[Path],
) -> None:
    row = coupling_rows[0]
    lines = [
        "# HHT resonance diagnostic",
        "",
        "## Purpose",
        "",
        "This is a diagnostic side test for the SU2/chirality programme. It checks whether a Hilbert-Huang workflow can recover a controlled SU2 -> chi_dot -> U(1) resonance chain before applying the method to noisy cosmology residuals or parity observables.",
        "",
        "It is not a SN/BAO/Planck likelihood and it is not a full U(2) model.",
        "",
        "## Synthetic setup",
        "",
        f"- omega_Q: {fmt_float(cfg.omega_q, 6)}",
        f"- m_chi: {fmt_float(cfg.m_chi, 6)}",
        f"- k_mode: {fmt_float(cfg.k_mode, 6)}",
        f"- first-band mismatch (2k - omega_Q)/omega_Q: {fmt_float(row['resonance_mismatch_fraction'], 6)}",
        f"- bias: {fmt_float(cfg.bias, 6)}",
        f"- damping: {fmt_float(cfg.damping, 6)}",
        f"- samples: {cfg.samples}",
        f"- duration periods: {fmt_float(cfg.duration_periods, 6)}",
        "",
        "## Coupling result",
        "",
        "| metric | value |",
        "|---|---:|",
        f"| Q-chi 1:1 phase-locking value | {fmt_float(row['Q_chi_1to1_PLV'], 5)} |",
        f"| chi-U1 plus 1:2 phase-locking value | {fmt_float(row['chi_U1_plus_1to2_PLV'], 5)} |",
        f"| chi-U1 minus 1:2 phase-locking value | {fmt_float(row['chi_U1_minus_1to2_PLV'], 5)} |",
        f"| max chi-U1 1:2 phase-locking value | {fmt_float(row['max_chi_U1_1to2_PLV'], 5)} |",
        f"| U1 plus energy growth | {fmt_float(row['U1_plus_energy_growth_late_over_early'], 5)} |",
        f"| U1 minus energy growth | {fmt_float(row['U1_minus_energy_growth_late_over_early'], 5)} |",
        f"| raw median final helicity asymmetry | {fmt_float(row['median_final_helicity_asymmetry'], 5)} |",
        f"| intrinsic chirality flag | {row['intrinsic_chirality_flag_requires_bias']} |",
        f"| resonance score | {fmt_float(row['resonance_coupling_score_0to1'], 5)} |",
        "",
        f"Interpretation: {row['interpretation']}.",
        "",
        "## Selected IMF components",
        "",
        "| series | selected IMF | selected frequency | expected frequency |",
        "|---|---:|---:|---:|",
        f"| Q | {row['Q_selected_imf']} | {fmt_float(row['Q_selected_frequency'], 5)} | {fmt_float(row['expected_f_Q_cycles_per_t'], 5)} |",
        f"| chi_dot | {row['chi_dot_selected_imf']} | {fmt_float(row['chi_dot_selected_frequency'], 5)} | {fmt_float(row['expected_f_Q_cycles_per_t'], 5)} |",
        f"| U1_plus | {row['U1_plus_selected_imf']} | {fmt_float(row['U1_plus_selected_frequency'], 5)} | {fmt_float(row['expected_f_U1_cycles_per_t'], 5)} |",
        f"| U1_minus | {row['U1_minus_selected_imf']} | {fmt_float(row['U1_minus_selected_frequency'], 5)} | {fmt_float(row['expected_f_U1_cycles_per_t'], 5)} |",
        "",
        "## Notes for project use",
        "",
        "- A high score here only proves that HHT can detect the synthetic resonance pattern.",
        "- Applying this to cosmology residuals needs surrogate tests: shuffled redshift/order controls, LCDM-vs-SU2 comparisons, and repeated noise realizations.",
        "- In the exact zero-bias periodic case, the SU2 -> chi -> U(1) mechanism is not intrinsically chiral; a raw finite-time helicity asymmetry should be treated as a phase-window effect unless a bias or other symmetry-breaking ingredient is present.",
        "- The natural next data-facing targets are residual series ordered by log(1+z), lookback time, chiral magnetic-field diagnostics, or parity-sensitive 4PCF/GW observables.",
        "",
        "## Output files",
        "",
        "- hht_resonance_synthetic_signals.csv",
        "- hht_resonance_imf_summary.csv",
        "- hht_resonance_coupling_summary.csv",
    ]
    for plot_path in plot_paths:
        lines.append(f"- {plot_path.name}")
    lines.append("")

    path.write_text("\n".join(lines), encoding="utf-8")


def write_csv_mode_summary(
    path: Path,
    cfg: CsvConfig,
    imf_rows: list[dict[str, object]],
    plot_paths: list[Path],
    tag: str,
) -> None:
    dominant = sorted(
        [r for r in imf_rows if str(r["component"]).startswith("IMF")],
        key=lambda r: float(r["energy_fraction"]),
        reverse=True,
    )[:5]

    lines = [
        "# HHT CSV diagnostic",
        "",
        f"Input: `{cfg.input_csv}`",
        f"x column: `{cfg.x_col}`",
        f"y column: `{cfg.y_col}`",
        f"reference column: `{cfg.ref_col or ''}`",
        "",
        "## Dominant IMF components",
        "",
        "| rank | series | component | energy fraction | weighted frequency | median positive frequency |",
        "|---:|---|---|---:|---:|---:|",
    ]
    if dominant:
        for i, row in enumerate(dominant, start=1):
            lines.append(
                f"| {i} | {row['series']} | {row['component']} | {fmt_float(row['energy_fraction'], 5)} | {fmt_float(row['amplitude_weighted_frequency'], 5)} | {fmt_float(row['median_positive_frequency'], 5)} |"
            )
    else:
        lines.append("| - | no IMF extracted | residue only | 1 | nan | nan |")
    lines.extend(
        [
            "",
            "## Caution",
            "",
            "HHT on cosmology residuals is exploratory. Treat any oscillatory component as a candidate feature only after comparing against LCDM/SU2 residuals, redshift shuffles, bootstrap resamples, and survey-window/systematic controls.",
            "",
            "## Output files",
            "",
            f"- hht_resonance_csv_{tag}_imf_summary.csv",
        ]
    )
    for plot_path in plot_paths:
        lines.append(f"- {plot_path.name}")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def maybe_plot_synthetic(
    signals: dict[str, np.ndarray],
    imfs_by_label: dict[str, list[np.ndarray]],
    cfg: SyntheticConfig,
) -> list[Path]:
    if cfg.no_plots:
        return []
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return []

    t = signals["t"]
    plot_paths: list[Path] = []

    overview = OUTDIR / "hht_resonance_overview.png"
    fig, axes = plt.subplots(4, 1, figsize=(11, 9), sharex=True)
    axes[0].plot(t, signals["Q"], lw=1.0)
    axes[0].set_ylabel("Q")
    axes[1].plot(t, signals["chi_dot"], lw=1.0, color="tab:green")
    axes[1].set_ylabel("chi_dot")
    axes[2].plot(t, signals["U1_plus_energy"], lw=1.0, label="U1 +")
    axes[2].plot(t, signals["U1_minus_energy"], lw=1.0, label="U1 -")
    axes[2].set_yscale("log")
    axes[2].set_ylabel("U1 energy")
    axes[2].legend(loc="upper left", fontsize=8)
    axes[3].plot(t, signals["helicity_asymmetry"], lw=1.0, color="tab:red")
    axes[3].axhline(0.0, color="black", lw=0.8)
    axes[3].set_ylabel("helicity asym")
    axes[3].set_xlabel("dimensionless time")
    fig.suptitle("Synthetic SU2 -> chi -> U1 resonance diagnostic")
    fig.tight_layout()
    fig.savefig(overview, dpi=160)
    plt.close(fig)
    plot_paths.append(overview)

    spectrum = OUTDIR / "hht_resonance_hilbert_spectrum.png"
    fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
    for ax, label in zip(axes, ["Q", "chi_dot", "U1_plus"]):
        imfs = imfs_by_label.get(label, [])
        for imf in imfs[: min(4, len(imfs))]:
            features = hilbert_features(t, imf)
            amp = features["amplitude"]
            freq = features["frequency"]
            mask = np.isfinite(freq) & (freq > 0)
            if mask.sum() == 0:
                continue
            ax.scatter(t[mask], freq[mask], c=amp[mask], s=3, cmap="viridis", alpha=0.55)
        ax.set_ylabel(f"{label}\nfreq")
        ax.set_ylim(bottom=0)
    axes[-1].set_xlabel("dimensionless time")
    fig.suptitle("HHT instantaneous-frequency traces")
    fig.tight_layout()
    fig.savefig(spectrum, dpi=160)
    plt.close(fig)
    plot_paths.append(spectrum)

    return plot_paths


def maybe_plot_csv(t: np.ndarray, series: dict[str, np.ndarray], imfs_by_label: dict[str, list[np.ndarray]], cfg: CsvConfig, tag: str) -> list[Path]:
    if cfg.no_plots:
        return []
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return []

    paths: list[Path] = []
    overview = OUTDIR / f"hht_resonance_csv_{tag}_overview.png"
    labels = list(series)
    fig, axes = plt.subplots(len(labels), 1, figsize=(11, max(4, 3 * len(labels))), sharex=True)
    if len(labels) == 1:
        axes = [axes]
    for ax, label in zip(axes, labels):
        ax.plot(t, series[label], lw=1.0, label=label)
        for imf in imfs_by_label.get(label, [])[:3]:
            ax.plot(t, imf, lw=0.8, alpha=0.65)
        ax.set_ylabel(label)
    axes[-1].set_xlabel(cfg.x_col)
    fig.suptitle("HHT CSV decomposition")
    fig.tight_layout()
    fig.savefig(overview, dpi=160)
    plt.close(fig)
    paths.append(overview)
    return paths


def run_synthetic(cfg: SyntheticConfig) -> None:
    ensure_outdir()
    signals = build_synthetic_signals(cfg)
    write_csv(OUTDIR / "hht_resonance_synthetic_signals.csv", signal_rows(signals))

    t = signals["t"]
    series = {
        "Q": signals["Q"],
        "chi_dot": signals["chi_dot"],
        "U1_plus": signals["U1_plus"],
        "U1_minus": signals["U1_minus"],
        "helicity_asymmetry": signals["helicity_asymmetry"],
    }
    imfs_by_label, _residue_by_label, imf_rows = decompose_series(
        t,
        series,
        max_imfs=cfg.max_imfs,
        sift_max_iter=cfg.sift_max_iter,
        sift_sd=cfg.sift_sd,
    )
    coupling_rows = synthetic_coupling_rows(cfg, signals, imfs_by_label, imf_rows)
    write_csv(OUTDIR / "hht_resonance_imf_summary.csv", imf_rows)
    write_csv(OUTDIR / "hht_resonance_coupling_summary.csv", coupling_rows)
    plot_paths = maybe_plot_synthetic(signals, imfs_by_label, cfg)
    write_synthetic_summary(OUTDIR / "hht_resonance_summary.md", cfg, coupling_rows, imf_rows, plot_paths)

    row = coupling_rows[0]
    print(f"Saved signals: {OUTDIR / 'hht_resonance_synthetic_signals.csv'}")
    print(f"Saved IMF summary: {OUTDIR / 'hht_resonance_imf_summary.csv'}")
    print(f"Saved coupling summary: {OUTDIR / 'hht_resonance_coupling_summary.csv'}")
    print(f"Saved summary: {OUTDIR / 'hht_resonance_summary.md'}")
    for path in plot_paths:
        print(f"Saved plot: {path}")
    print(
        "Resonance score="
        f"{float(row['resonance_coupling_score_0to1']):.3f}, "
        f"Q-chi PLV={float(row['Q_chi_1to1_PLV']):.3f}, "
        f"chi-U1+ PLV={float(row['chi_U1_plus_1to2_PLV']):.3f}, "
        f"max chi-U1 PLV={float(row['max_chi_U1_1to2_PLV']):.3f}, "
        f"growth+={float(row['U1_plus_energy_growth_late_over_early']):.3g}, "
        f"raw final asym={float(row['median_final_helicity_asymmetry']):.3g}, "
        f"intrinsic chirality={row['intrinsic_chirality_flag_requires_bias']}"
    )


def run_csv_mode(cfg: CsvConfig) -> None:
    ensure_outdir()
    cols = read_csv_columns(cfg.input_csv)
    if cfg.x_col not in cols:
        raise ValueError(f"x column {cfg.x_col!r} not found in {cfg.input_csv}")
    if cfg.y_col not in cols:
        raise ValueError(f"y column {cfg.y_col!r} not found in {cfg.input_csv}")
    if cfg.ref_col and cfg.ref_col not in cols:
        raise ValueError(f"reference column {cfg.ref_col!r} not found in {cfg.input_csv}")

    x_raw = as_float_array(cols[cfg.x_col], cfg.x_col)
    y_raw = as_float_array(cols[cfg.y_col], cfg.y_col)
    order = np.argsort(x_raw)
    x_raw = x_raw[order]
    y_raw = y_raw[order]
    unique = np.concatenate(([True], np.diff(x_raw) > 0))
    x_raw = x_raw[unique]
    y_raw = y_raw[unique]

    samples = min(cfg.samples, x_raw.size) if cfg.samples > 0 else x_raw.size
    if samples < 8:
        samples = x_raw.size
    t = np.linspace(float(x_raw[0]), float(x_raw[-1]), samples)
    y = np.interp(t, x_raw, y_raw)
    series = {cfg.label: y}

    if cfg.ref_col:
        ref_raw = as_float_array(cols[cfg.ref_col], cfg.ref_col)[order][unique]
        series[f"{cfg.ref_col}_reference"] = np.interp(t, x_raw, ref_raw)

    imfs_by_label, _residue_by_label, imf_rows = decompose_series(
        t,
        series,
        max_imfs=cfg.max_imfs,
        sift_max_iter=cfg.sift_max_iter,
        sift_sd=cfg.sift_sd,
    )
    tag = safe_slug(cfg.label or cfg.y_col)
    imf_path = OUTDIR / f"hht_resonance_csv_{tag}_imf_summary.csv"
    summary_path = OUTDIR / f"hht_resonance_csv_{tag}_summary.md"
    write_csv(imf_path, imf_rows)
    plot_paths = maybe_plot_csv(t, series, imfs_by_label, cfg, tag)
    write_csv_mode_summary(summary_path, cfg, imf_rows, plot_paths, tag)

    # Keep a latest-run copy for convenience while preserving the label-specific files.
    write_csv(OUTDIR / "hht_resonance_csv_imf_summary.csv", imf_rows)
    write_csv_mode_summary(OUTDIR / "hht_resonance_csv_summary.md", cfg, imf_rows, plot_paths, tag)

    print(f"Saved CSV IMF summary: {imf_path}")
    print(f"Saved CSV summary: {summary_path}")
    print(f"Updated latest-run copies: {OUTDIR / 'hht_resonance_csv_imf_summary.csv'}, {OUTDIR / 'hht_resonance_csv_summary.md'}")
    for path in plot_paths:
        print(f"Saved plot: {path}")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Diagnose SU2 -> chi -> U(1) resonance coupling with HHT/EMD.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--mode", choices=["synthetic", "csv"], default="synthetic")

    parser.add_argument("--omega-q", type=float, default=1.0, help="Dimensionless SU2 condensate angular frequency.")
    parser.add_argument("--m-chi", type=float, default=1.05, help="Dimensionless pseudoscalar mass proxy.")
    parser.add_argument("--k-mode", type=float, default=0.5, help="Dimensionless U(1) mode wavenumber.")
    parser.add_argument("--su2-chi-coupling", type=float, default=0.12, help="Strength of Q^2 Q_dot sourcing chi.")
    parser.add_argument("--u1-coupling", type=float, default=0.08, help="Strength of chi_dot modulation of U(1) modes.")
    parser.add_argument("--bias", type=float, default=0.0, help="Constant chi_dot bias; nonzero values can break helicity symmetry.")
    parser.add_argument("--damping", type=float, default=0.0, help="Uniform U(1) damping term in the toy mode equation.")
    parser.add_argument("--duration-periods", type=float, default=80.0)
    parser.add_argument("--samples", type=int, default=8192)
    parser.add_argument("--q0", type=float, default=1.0)
    parser.add_argument("--anharmonic", type=float, default=0.05)
    parser.add_argument("--noise", type=float, default=0.0)
    parser.add_argument("--seed", type=int, default=260520333)

    parser.add_argument("--input-csv", type=Path, default=None, help="CSV file for mode=csv.")
    parser.add_argument("--x-col", default="x", help="Independent variable column for mode=csv.")
    parser.add_argument("--y-col", default="y", help="Signal/residual column for mode=csv.")
    parser.add_argument("--ref-col", default=None, help="Optional reference signal column for mode=csv.")
    parser.add_argument("--label", default="series", help="Series label for mode=csv.")

    parser.add_argument("--max-imfs", type=int, default=6)
    parser.add_argument("--sift-max-iter", type=int, default=80)
    parser.add_argument("--sift-sd", type=float, default=0.08)
    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument("--no-plots", action="store_true")
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()

    if args.mode == "synthetic":
        if args.samples < 128:
            raise ValueError("--samples must be at least 128 in synthetic mode")
        cfg = SyntheticConfig(
            omega_q=args.omega_q,
            m_chi=args.m_chi,
            k_mode=args.k_mode,
            su2_chi_coupling=args.su2_chi_coupling,
            u1_coupling=args.u1_coupling,
            bias=args.bias,
            damping=args.damping,
            duration_periods=args.duration_periods,
            samples=args.samples,
            q0=args.q0,
            anharmonic=args.anharmonic,
            noise=args.noise,
            seed=args.seed,
            max_imfs=args.max_imfs,
            sift_max_iter=args.sift_max_iter,
            sift_sd=args.sift_sd,
            top_n=args.top_n,
            no_plots=args.no_plots,
        )
        run_synthetic(cfg)
        return

    if args.input_csv is None:
        raise ValueError("--input-csv is required in mode=csv")
    cfg = CsvConfig(
        input_csv=args.input_csv,
        x_col=args.x_col,
        y_col=args.y_col,
        ref_col=args.ref_col,
        label=args.label,
        samples=args.samples,
        max_imfs=args.max_imfs,
        sift_max_iter=args.sift_max_iter,
        sift_sd=args.sift_sd,
        no_plots=args.no_plots,
    )
    run_csv_mode(cfg)


if __name__ == "__main__":
    main()
