#!/usr/bin/env python3
r"""Run a literature-defined high-temperature SU(2)-CMB benchmark in LINX.

This is deliberately separate from the repository's late-time SU2/SU2R
effective-fluid branch.  It implements the asymptotic deconfining SU(2)-CMB
thermodynamic limit described by Hahn, Hofmann, and Kramer:

    g_YM / g_gamma = 8 / 2 = 4
    nu_CMB          = (1 / 4)^(1 / 3)
    T_gamma(a)      = nu_CMB T0 / a        at T >> T0

At BBN temperatures the six additional gauge polarisations are treated as
ideal, ultra-relativistic bosonic modes sharing the electromagnetic-sector
temperature.  They affect H(T), electromagnetic heat capacity, neutrino
decoupling, and hence the standard LINX weak rates.  Nuclear masses, binding
energies, Q-values, and PRIMAT reaction rates remain standard with q_i = 0,
because the cited SU(2)-CMB thermodynamic model does not predict their shifts.

This is a physical thermodynamic benchmark, not a complete gauge-matter BBN
calculation.  In particular, SU(2)-specific radiative corrections to charged
leptons and nuclei are not available in the cited model or LINX interface.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCRIPT_PATH = Path(__file__).resolve()
ROOT = SCRIPT_PATH.parents[3]
TASK_ROOT = Path(r"C:\Users\clive\Documents\Codex\2026-07-14\ok-continu")
DEFAULT_LINX_PATH = ROOT / "plamb_runs" / "deps" / "LINX"
DEFAULT_OUTDIR = ROOT / "plamb_runs" / "diagnostics" / "su2cmb_linx_thermal_bridge_20260717"
DEFAULT_OUTPUTS = TASK_ROOT / "outputs"
RUN_DATE = "2026-07-17"

G_GAMMA = 2.0
G_YM = 8.0
RHO_YM_OVER_RHO_GAMMA_AT_FIXED_T = G_YM / G_GAMMA
EXTRA_RHO_OVER_RHO_GAMMA = RHO_YM_OVER_RHO_GAMMA_AT_FIXED_T - 1.0
NU_CMB = (G_GAMMA / G_YM) ** (1.0 / 3.0)

OBSERVED = {
    "D_H_1e5": (2.508, 0.029),
    "Yp_mass": (0.245, 0.003),
    "Li7_H_1e10": (1.45, 0.25),
}

# The SU(2)-CMB angular-power-spectrum paper used this prior range for omega_b.
SU2CMB_OMEGAB_H2_RANGE = (0.014, 0.027)
LINX_OMEGAB_H2_REFERENCE = 0.02242

SOURCES = {
    "su2cmb_temperature_relation": "https://arxiv.org/abs/1712.08561",
    "su2cmb_cosmological_model": "https://arxiv.org/abs/1810.01253",
    "linx": "https://arxiv.org/abs/2408.14538",
    "linx_repository": "https://github.com/cgiovanetti/LINX",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--linx-path", type=Path, default=DEFAULT_LINX_PATH)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--outputs", type=Path, default=DEFAULT_OUTPUTS)
    parser.add_argument("--network", choices=["key", "small"], default="key")
    parser.add_argument("--eta-min", type=float, default=0.08)
    parser.add_argument("--eta-max", type=float, default=1.22)
    parser.add_argument("--eta-n", type=int, default=58)
    parser.add_argument("--sampling-ntop", type=int, default=150)
    parser.add_argument("--rtol", type=float, default=1e-5)
    parser.add_argument("--atol", type=float, default=1e-9)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_value(repo: Path, *args: str) -> str:
    try:
        result = subprocess.run(
            ["git", "-C", str(repo), *args],
            check=True,
            capture_output=True,
            text=True,
        )
        return result.stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def import_linx(linx_path: Path) -> SimpleNamespace:
    if not linx_path.exists():
        raise FileNotFoundError(f"LINX source path not found: {linx_path}")
    sys.path.insert(0, str(linx_path))

    import equinox as eqx  # type: ignore
    import jax  # type: ignore
    import jax.numpy as jnp  # type: ignore
    from diffrax import Event, ODETerm, PIDController, SaveAt, Tsit5, diffeqsolve  # type: ignore
    from linx import const, thermo  # type: ignore
    from linx.abundances import AbundanceModel  # type: ignore
    from linx.background import BackgroundModel  # type: ignore
    from linx.nuclear import NuclearRates  # type: ignore

    jax.config.update("jax_enable_x64", True)
    return SimpleNamespace(
        eqx=eqx,
        jax=jax,
        jnp=jnp,
        Event=Event,
        ODETerm=ODETerm,
        PIDController=PIDController,
        SaveAt=SaveAt,
        Tsit5=Tsit5,
        diffeqsolve=diffeqsolve,
        const=const,
        thermo=thermo,
        AbundanceModel=AbundanceModel,
        BackgroundModel=BackgroundModel,
        NuclearRates=NuclearRates,
    )


def make_su2cmb_background(api: SimpleNamespace):
    eqx = api.eqx
    jax = api.jax
    jnp = api.jnp
    thermo = api.thermo
    const = api.const
    rho_massless_be_v = jax.vmap(thermo.rho_massless_BE, in_axes=(0, None, None))
    rho_massless_fd_v = jax.vmap(thermo.rho_massless_FD, in_axes=(0, None, None))

    class SU2CMBAsymptoticBackground(eqx.Module):
        """LINX-compatible asymptotic deconfining SU(2)-CMB background."""

        decoupled: bool
        use_FD: bool
        collision_me: bool
        LO: bool
        NLO: bool
        max_steps: int
        throw: bool

        def __init__(self, max_steps: int = 512, throw: bool = True) -> None:
            self.decoupled = False
            self.use_FD = True
            self.collision_me = True
            self.LO = True
            self.NLO = True
            self.max_steps = max_steps
            self.throw = throw

        @eqx.filter_jit
        def __call__(
            self,
            T_start=const.T_start,
            T_end=const.T_end,
            me=const.me,
            rtol=1e-8,
            atol=1e-10,
            solver=api.Tsit5(),
        ):
            lna_init = 0.0
            T_g_init = T_start
            T_nu_init = T_start
            initial_g_star = 10.75 + (G_YM - G_GAMMA)
            t0 = (1.5 / T_start * initial_g_star ** (-0.25)) ** 2
            y0 = (lna_init, T_g_init, T_nu_init)

            def temperature_event(t, y, args, **kwargs):
                del t, args, kwargs
                return y[1] < T_end

            solution = api.diffeqsolve(
                api.ODETerm(self.dY),
                solver,
                args=(lna_init, me),
                t0=t0,
                t1=jnp.inf,
                dt0=None,
                y0=y0,
                saveat=api.SaveAt(steps=True),
                event=api.Event(temperature_event),
                stepsize_controller=api.PIDController(rtol=rtol, atol=atol),
                max_steps=self.max_steps,
                throw=self.throw,
            )

            a_vec = jnp.exp(solution.ys[0])
            rho_g_vec = rho_massless_be_v(solution.ys[1], 0.0, G_GAMMA)
            rho_nu_vec = rho_massless_fd_v(solution.ys[2], 0.0, 2.0)

            last_step = jnp.max(
                jnp.argwhere(solution.ys[1] < T_end, size=self.max_steps)[:, 0]
            )
            t_vec = jnp.where(solution.ts == jnp.inf, solution.ts[last_step], solution.ts)
            a_vec = jnp.where(a_vec == jnp.inf, a_vec[last_step], a_vec)
            rho_g_vec = jnp.where(rho_g_vec == jnp.inf, rho_g_vec[last_step], rho_g_vec)
            rho_nu_vec = jnp.where(rho_nu_vec == jnp.inf, rho_nu_vec[last_step], rho_nu_vec)

            # At BBN temperatures the published high-T entropy relation applies.
            final_a = NU_CMB * const.T0CMB / solution.ys[1][last_step]
            a_vec = a_vec * final_a / a_vec[-1]

            rho_extra_vec = EXTRA_RHO_OVER_RHO_GAMMA * rho_g_vec
            p_extra_vec = rho_extra_vec / 3.0
            T_g_vec = thermo.T_g(rho_g_vec)
            rho_total_vec = (
                thermo.rho_EM_std_v(T_g_vec) + 3.0 * rho_nu_vec + rho_extra_vec
            )
            neff_vec = thermo.N_eff(rho_total_vec, rho_g_vec)
            return (
                t_vec,
                a_vec,
                rho_g_vec,
                rho_nu_vec,
                rho_extra_vec,
                p_extra_vec,
                neff_vec,
            )

        @eqx.filter_jit
        def dY(self, t, y, args):
            del t
            lna, T_g, T_nu = y
            lna_init, me = args
            del lna_init

            rho_em = thermo.rho_EM_std(T_g, me=me, LO=self.LO, NLO=self.NLO)
            rho_plus_p_em = thermo.rho_plus_p_EM_std(
                T_g, me=me, LO=self.LO, NLO=self.NLO
            )
            drho_em_dT = thermo.drho_EM_dT_g_std(
                T_g, me=me, LO=self.LO, NLO=self.NLO
            )

            rho_g = thermo.rho_massless_BE(T_g, 0.0, G_GAMMA)
            rho_su2_extra = EXTRA_RHO_OVER_RHO_GAMMA * rho_g
            rho_plus_p_su2_extra = (4.0 / 3.0) * rho_su2_extra
            drho_su2_extra_dT = 4.0 * rho_su2_extra / T_g

            rho_nu = 3.0 * thermo.rho_nue_std(T_nu)
            rho_plus_p_nu = (4.0 / 3.0) * rho_nu
            drho_nu_dT = 3.0 * thermo.drho_nue_dT_nue_std(T_nu)

            H = thermo.Hubble(rho_em + rho_su2_extra + rho_nu)
            C_nue, C_numu, _, _ = thermo.collision_terms_std(
                T_g,
                T_nu,
                T_nu,
                me=me,
                decoupled=self.decoupled,
                use_FD=self.use_FD,
                collision_me=self.collision_me,
            )

            drho_coupled_dt = (
                -3.0 * H * (rho_plus_p_em + rho_plus_p_su2_extra)
                - C_nue
                - 2.0 * C_numu
            )
            drho_nu_dt = -3.0 * H * rho_plus_p_nu + C_nue + 2.0 * C_numu
            dT_g_dt = drho_coupled_dt / (drho_em_dT + drho_su2_extra_dT)
            dT_nu_dt = drho_nu_dt / drho_nu_dT
            return H, dT_g_dt, dT_nu_dt

    return SU2CMBAsymptoticBackground()


def unique_trajectory(api: SimpleNamespace, model: str, background: tuple[Any, ...]) -> pd.DataFrame:
    t_vec, a_vec, rho_g, rho_nu, rho_np, p_np, neff_vec = background
    del p_np
    t = np.asarray(t_vec, dtype=float)
    a = np.asarray(a_vec, dtype=float)
    rho_g_np = np.asarray(rho_g, dtype=float)
    rho_nu_np = np.asarray(rho_nu, dtype=float)
    rho_extra_np = np.asarray(rho_np, dtype=float)
    neff = np.asarray(neff_vec, dtype=float)
    T_g = np.asarray(api.thermo.T_g(api.jnp.asarray(rho_g_np)), dtype=float)
    T_nu = np.asarray(api.thermo.T_nu(api.jnp.asarray(rho_nu_np)), dtype=float)
    rho_total = np.asarray(
        api.thermo.rho_EM_std_v(api.jnp.asarray(T_g))
        + 3.0 * api.jnp.asarray(rho_nu_np)
        + api.jnp.asarray(rho_extra_np),
        dtype=float,
    )
    H = np.asarray(api.thermo.Hubble(api.jnp.asarray(rho_total)), dtype=float)
    frame = pd.DataFrame(
        {
            "model": model,
            "time_s": t,
            "scale_factor": a,
            "T_gamma_MeV": T_g,
            "T_nu_MeV": T_nu,
            "T_nu_over_T_gamma": T_nu / T_g,
            "H_s_inv": H,
            "rho_gamma_MeV4": rho_g_np,
            "rho_extra_MeV4": rho_extra_np,
            "neff_linx_diagnostic": neff,
            "a_Tgamma_over_T0": a * T_g / float(api.const.T0CMB),
        }
    )
    frame = frame.replace([np.inf, -np.inf], np.nan).dropna()
    frame = frame.drop_duplicates(subset=["time_s", "T_gamma_MeV"], keep="first")
    return frame.sort_values("T_gamma_MeV", ascending=False).reset_index(drop=True)


def abundance_values(y: Any) -> dict[str, float]:
    values = np.asarray([float(value) for value in tuple(y)], dtype=float)
    proton = values[1]
    deuterium = values[2]
    tritium = values[3] if len(values) > 3 else 0.0
    helium3 = values[4] if len(values) > 4 else 0.0
    helium4 = values[5] if len(values) > 5 else 0.0
    lithium7 = values[6] if len(values) > 6 else 0.0
    beryllium7 = values[7] if len(values) > 7 else 0.0
    lithium6 = values[8] if len(values) > 8 else 0.0
    return {
        "Yp_mass": 4.0 * helium4,
        "D_H_1e5": deuterium / proton * 1e5,
        "He3_H_1e5": (tritium + helium3) / proton * 1e5,
        "Li7_H_1e10": (lithium7 + beryllium7) / proton * 1e10,
        "Li6_H_1e14": lithium6 / proton * 1e14,
    }


def score_abundances(row: dict[str, Any]) -> None:
    chi2 = 0.0
    for key, (mean, sigma) in OBSERVED.items():
        z = (float(row[key]) - mean) / sigma
        row[f"z_{key}"] = z
        row[f"chi2_{key}"] = z * z
        chi2 += z * z
    row["chi2_abundances_measurement"] = chi2
    row["chi2_D_He"] = row["chi2_D_H_1e5"] + row["chi2_Yp_mass"]
    row["D_He_each_2sigma"] = (
        abs(row["z_D_H_1e5"]) <= 2.0 and abs(row["z_Yp_mass"]) <= 2.0
    )
    row["D_He_Li_each_2sigma"] = (
        row["D_He_each_2sigma"] and abs(row["z_Li7_H_1e10"]) <= 2.0
    )
    row["stellar_depletion_to_plateau"] = row["Li7_H_1e10"] / OBSERVED["Li7_H_1e10"][0]


def run_abundance(
    api: SimpleNamespace,
    model: Any,
    reaction_names: list[str],
    background: tuple[Any, ...],
    eta_fac: float,
    sampling_ntop: int,
    rtol: float,
    atol: float,
) -> dict[str, float]:
    t_vec, a_vec, rho_g, rho_nu, rho_np, p_np, _ = background
    q = api.jnp.zeros(len(reaction_names))
    y = model(
        api.jnp.asarray(rho_g),
        api.jnp.asarray(rho_nu),
        api.jnp.asarray(rho_np),
        api.jnp.asarray(p_np),
        t_vec=api.jnp.asarray(t_vec),
        a_vec=api.jnp.asarray(a_vec),
        eta_fac=api.jnp.asarray(eta_fac),
        tau_n_fac=api.jnp.asarray(1.0),
        nuclear_rates_q=q,
        sampling_nTOp=sampling_ntop,
        rtol=rtol,
        atol=atol,
    )
    return abundance_values(y)


def selected_summary(scan: pd.DataFrame) -> pd.DataFrame:
    su2 = scan.loc[scan["model"].eq("SU2CMB_asymptotic")].copy()
    prior_low = SU2CMB_OMEGAB_H2_RANGE[0] / LINX_OMEGAB_H2_REFERENCE
    prior_high = SU2CMB_OMEGAB_H2_RANGE[1] / LINX_OMEGAB_H2_REFERENCE
    selectors = {
        "SBBN_reference": scan["model"].eq("SBBN"),
        "SU2CMB_fixed_reference_eta": su2["eta_fac"].sub(1.0).abs().eq(
            su2["eta_fac"].sub(1.0).abs().min()
        ),
        "SU2CMB_entropy_compensated_eta": su2["eta_fac"].sub(0.25).abs().eq(
            su2["eta_fac"].sub(0.25).abs().min()
        ),
        "SU2CMB_best_D_He_any_eta": su2["chi2_D_He"].eq(su2["chi2_D_He"].min()),
        "SU2CMB_best_all_abundances_any_eta": su2["chi2_abundances_measurement"].eq(
            su2["chi2_abundances_measurement"].min()
        ),
        "SU2CMB_best_in_published_omega_b_prior": (
            su2["eta_fac"].between(prior_low, prior_high)
            & su2["chi2_abundances_measurement"].eq(
                su2.loc[su2["eta_fac"].between(prior_low, prior_high), "chi2_abundances_measurement"].min()
            )
        ),
    }
    rows: list[pd.Series] = []
    for label, selector in selectors.items():
        source = scan if label == "SBBN_reference" else su2
        matches = source.loc[selector]
        if matches.empty:
            continue
        row = matches.iloc[0].copy()
        row["summary_case"] = label
        rows.append(row)
    summary = pd.DataFrame(rows)
    columns = ["summary_case"] + [column for column in scan.columns if column != "summary_case"]
    return summary[columns]


def interpolate_log(x: np.ndarray, xp: np.ndarray, fp: np.ndarray) -> np.ndarray:
    order = np.argsort(xp)
    return np.exp(np.interp(np.log(x), np.log(xp[order]), np.log(fp[order])))


def write_plot(trajectory: pd.DataFrame, scan: pd.DataFrame, path: Path) -> None:
    standard = trajectory.loc[trajectory["model"].eq("SBBN")]
    su2 = trajectory.loc[trajectory["model"].eq("SU2CMB_asymptotic")]
    H_standard_at_su2 = interpolate_log(
        su2["T_gamma_MeV"].to_numpy(),
        standard["T_gamma_MeV"].to_numpy(),
        standard["H_s_inv"].to_numpy(),
    )
    su2_scan = scan.loc[scan["model"].eq("SU2CMB_asymptotic")].sort_values("eta_fac")

    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.2), constrained_layout=True)
    axes[0, 0].semilogx(
        su2["T_gamma_MeV"], su2["H_s_inv"].to_numpy() / H_standard_at_su2, color="#b23a48"
    )
    axes[0, 0].axhline(1.0, color="black", linewidth=0.8)
    axes[0, 0].set_xlabel(r"$T_\gamma$ [MeV]")
    axes[0, 0].set_ylabel(r"$H_{\rm SU2CMB}/H_{\rm SBBN}$")
    axes[0, 0].set_title("Expansion history")

    for model, colour in (("SBBN", "#2878b5"), ("SU2CMB_asymptotic", "#b23a48")):
        subset = trajectory.loc[trajectory["model"].eq(model)]
        axes[0, 1].semilogx(
            subset["T_gamma_MeV"], subset["a_Tgamma_over_T0"], label=model, color=colour
        )
    axes[0, 1].axhline(NU_CMB, color="#b23a48", linewidth=0.8, linestyle="--")
    axes[0, 1].set_xlabel(r"$T_\gamma$ [MeV]")
    axes[0, 1].set_ylabel(r"$aT_\gamma/T_0$")
    axes[0, 1].set_title("Temperature-scale-factor relation")
    axes[0, 1].legend(frameon=False)

    for model, colour in (("SBBN", "#2878b5"), ("SU2CMB_asymptotic", "#b23a48")):
        subset = trajectory.loc[trajectory["model"].eq(model)]
        axes[1, 0].semilogx(
            subset["T_gamma_MeV"], subset["T_nu_over_T_gamma"], label=model, color=colour
        )
    axes[1, 0].set_xlabel(r"$T_\gamma$ [MeV]")
    axes[1, 0].set_ylabel(r"$T_\nu/T_\gamma$")
    axes[1, 0].set_title("Neutrino-temperature response")

    abundance_styles = (
        ("D_H_1e5", "D/H x 1e5", "#2878b5"),
        ("Li7_H_1e10", "Li-7/H x 1e10", "#b23a48"),
    )
    for key, label, colour in abundance_styles:
        mean, sigma = OBSERVED[key]
        axes[1, 1].plot(su2_scan["eta_fac"], su2_scan[key], label=label, color=colour)
        axes[1, 1].axhspan(mean - 2 * sigma, mean + 2 * sigma, color=colour, alpha=0.10)
        axes[1, 1].axhline(mean, color=colour, linewidth=0.8, linestyle="--")
    axes[1, 1].axvspan(
        SU2CMB_OMEGAB_H2_RANGE[0] / LINX_OMEGAB_H2_REFERENCE,
        SU2CMB_OMEGAB_H2_RANGE[1] / LINX_OMEGAB_H2_REFERENCE,
        color="#666666",
        alpha=0.09,
        label="Published omega_b prior",
    )
    axes[1, 1].set_xlabel(r"Present baryon normalisation $\eta_{\rm fac}$")
    axes[1, 1].set_ylabel("Abundance in displayed units")
    axes[1, 1].set_title("Baryon-normalisation scan")
    axes[1, 1].legend(frameon=False, fontsize=8)

    fig.suptitle("Literature-defined SU(2)-CMB to LINX thermodynamic bridge", fontsize=13)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def fmt(value: Any, digits: int = 5) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "NA"
    return f"{number:.{digits}g}" if math.isfinite(number) else "NA"


def write_report(path: Path, trajectory: pd.DataFrame, summary: pd.DataFrame, scan: pd.DataFrame) -> None:
    by_case = summary.set_index("summary_case")
    standard = by_case.loc["SBBN_reference"]
    fixed = by_case.loc["SU2CMB_fixed_reference_eta"]
    compensated = by_case.loc["SU2CMB_entropy_compensated_eta"]
    best_dhe = by_case.loc["SU2CMB_best_D_He_any_eta"]
    best_all = by_case.loc["SU2CMB_best_all_abundances_any_eta"]
    best_prior = by_case.loc["SU2CMB_best_in_published_omega_b_prior"]
    standard_traj = trajectory.loc[trajectory["model"].eq("SBBN")]
    su2_traj = trajectory.loc[trajectory["model"].eq("SU2CMB_asymptotic")]
    target_T = 1.0
    H_std = interpolate_log(
        np.asarray([target_T]),
        standard_traj["T_gamma_MeV"].to_numpy(),
        standard_traj["H_s_inv"].to_numpy(),
    )[0]
    H_su2 = interpolate_log(
        np.asarray([target_T]),
        su2_traj["T_gamma_MeV"].to_numpy(),
        su2_traj["H_s_inv"].to_numpy(),
    )[0]
    joint_passes = int(scan.loc[scan["model"].eq("SU2CMB_asymptotic"), "D_He_Li_each_2sigma"].sum())

    lines = [
        "# SU(2)-CMB to LINX Thermodynamic Bridge",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "## Identity and Scope",
        "",
        "This calculation is a literature-defined high-temperature SU(2)-CMB benchmark. It is not the repository's late-time SU2/SU2R effective dark-sector fluid, and no equivalence between those models is assumed.",
        "",
        "The implemented asymptotic relations are",
        "",
        r"\[",
        r"\begin{aligned}",
        r"g_\gamma                         &= 2, \\",
        r"g_{\rm YM}                       &= 8, \\",
        r"\xi                              &= g_{\rm YM}/g_\gamma = 4, \\",
        r"\nu_{\rm CMB}                    &= \xi^{-1/3}=(1/4)^{1/3}, \\",
        r"T_\gamma(a)                      &= \nu_{\rm CMB}T_0/a, \\",
        r"\rho_{\rm YM}(T)                 &= \xi\rho_\gamma(T), \\",
        r"H^2(T)                           &= \frac{8\pi G}{3}\left(\rho_{\rm EM}+\rho_\nu+\rho_{\rm YM,extra}\right).",
        r"\end{aligned}",
        r"\]",
        "",
        "The six additional gauge polarisations are ideal ultra-relativistic modes at BBN temperature. Their energy density and heat capacity enter the LINX background integration. Standard neutrino-electron collision terms then predict a changed T_nu/T_gamma history, which is passed to the standard LINX neutron-proton weak rates.",
        "",
        "All PRIMAT nuclear-rate pulls are fixed to zero. Nuclear masses, binding energies and Q-values are standard because the cited SU(2)-CMB thermodynamic model supplies no predicted modifications to them.",
        "",
        "## Background Readout",
        "",
        f"- High-temperature coefficient: `nu_CMB = {NU_CMB:.9f}`.",
        f"- Expansion ratio at T_gamma = 1 MeV: `H_SU2CMB/H_SBBN = {H_su2 / H_std:.6g}`.",
        f"- Final LINX diagnostic N_eff: SBBN `{fmt(standard_traj.iloc[-1]['neff_linx_diagnostic'])}`; SU2-CMB `{fmt(su2_traj.iloc[-1]['neff_linx_diagnostic'])}`.",
        "",
        "## Abundance Readout",
        "",
        "| Case | eta_fac | post-e-annihilation baryon factor | D/H x 1e5 | Yp | Li-7/H x 1e10 | chi2 D+He | chi2 all | depletion |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    labels = {
        "SBBN_reference": "SBBN reference",
        "SU2CMB_fixed_reference_eta": "SU2-CMB, eta_fac=1",
        "SU2CMB_entropy_compensated_eta": "SU2-CMB, eta_fac near 1/4",
        "SU2CMB_best_D_He_any_eta": "SU2-CMB, best D+He",
        "SU2CMB_best_all_abundances_any_eta": "SU2-CMB, best D+He+Li",
        "SU2CMB_best_in_published_omega_b_prior": "SU2-CMB, best in published omega_b prior",
    }
    for case, row in by_case.iterrows():
        lines.append(
            f"| {labels[case]} | {fmt(row['eta_fac'])} | {fmt(row['post_eann_baryon_density_factor'])} | "
            f"{fmt(row['D_H_1e5'])} | {fmt(row['Yp_mass'])} | {fmt(row['Li7_H_1e10'])} | "
            f"{fmt(row['chi2_D_He'])} | {fmt(row['chi2_abundances_measurement'])} | "
            f"{fmt(row['stellar_depletion_to_plateau'])} |"
        )
    lines.extend(
        [
            "",
            "## Gate Interpretation",
            "",
            f"- Simultaneous D, He and Li two-sigma passes in the SU2-CMB eta scan: `{joint_passes}`.",
            f"- At fixed present reference eta, SU2-CMB predicts D/H x 1e5 `{fmt(fixed['D_H_1e5'])}`, Yp `{fmt(fixed['Yp_mass'])}`, and Li/H x 1e10 `{fmt(fixed['Li7_H_1e10'])}`.",
            f"- Cancelling the fourfold entropy-induced baryon-density increase requires eta_fac near one quarter; that control gives D/H x 1e5 `{fmt(compensated['D_H_1e5'])}` and Li/H x 1e10 `{fmt(compensated['Li7_H_1e10'])}`.",
            f"- The unrestricted best D+He row occurs at eta_fac `{fmt(best_dhe['eta_fac'])}`; the unrestricted direct-lithium best row occurs at eta_fac `{fmt(best_all['eta_fac'])}`.",
            f"- Inside the published SU2-CMB omega_b prior, the best abundance chi-squared is `{fmt(best_prior['chi2_abundances_measurement'])}` at eta_fac `{fmt(best_prior['eta_fac'])}`.",
            "",
            "## Claim Boundary",
            "",
            "This bridge implements the published asymptotic SU(2)-CMB equation-of-state limit and propagates it through LINX background and standard weak-rate machinery. It is not a complete first-principles SU(2) gauge-matter BBN calculation: the cited model does not provide SU(2)-specific electron radiative corrections, nuclear binding shifts, Q-value shifts, or reaction-rate changes. Those inputs remain standard and are not fitted.",
            "",
            "The repository's SU2/SU2R late-time fluid must remain scientifically separate unless an action-level derivation explicitly connects it to SU(2)-CMB thermodynamics.",
            "",
            "## Sources",
            "",
            f"- High-temperature SU(2)-CMB relation: {SOURCES['su2cmb_temperature_relation']}",
            f"- Eight-mode cosmological implementation and omega_b prior: {SOURCES['su2cmb_cosmological_model']}",
            f"- LINX BBN framework: {SOURCES['linx']}",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    cli = parse_args()
    cli.outdir.mkdir(parents=True, exist_ok=True)
    cli.outputs.mkdir(parents=True, exist_ok=True)
    completion_marker = cli.outdir / "su2cmb_linx_thermal_bridge_complete.json"
    if completion_marker.exists() and not cli.force:
        print(f"Complete result already exists: {completion_marker}")
        return 0

    print("Loading LINX", flush=True)
    api = import_linx(cli.linx_path)
    print(f"JAX devices: {api.jax.devices()}", flush=True)
    network_name = "key_PRIMAT_2023" if cli.network == "key" else "small_PRIMAT_2023"
    rates = api.NuclearRates(nuclear_net=network_name)
    abundance_model = api.AbundanceModel(rates, throw=False)
    reaction_names = list(rates.reactions_names)

    print("Computing SBBN background", flush=True)
    standard_background = api.BackgroundModel()(api.jnp.asarray(0.0))
    print("Computing asymptotic SU2-CMB background", flush=True)
    su2_background = make_su2cmb_background(api)()

    trajectory = pd.concat(
        [
            unique_trajectory(api, "SBBN", standard_background),
            unique_trajectory(api, "SU2CMB_asymptotic", su2_background),
        ],
        ignore_index=True,
    )

    rows: list[dict[str, Any]] = []
    print("Computing SBBN abundance reference", flush=True)
    standard_row: dict[str, Any] = {
        "model": "SBBN",
        "eta_fac": 1.0,
        "omega_b_h2": LINX_OMEGAB_H2_REFERENCE,
        "post_eann_baryon_density_factor": 1.0,
    }
    standard_row.update(
        run_abundance(
            api,
            abundance_model,
            reaction_names,
            standard_background,
            1.0,
            cli.sampling_ntop,
            cli.rtol,
            cli.atol,
        )
    )
    score_abundances(standard_row)
    rows.append(standard_row)

    eta_values = np.unique(
        np.concatenate(
            [
                np.linspace(cli.eta_min, cli.eta_max, cli.eta_n),
                np.asarray([0.25, 1.0]),
                np.asarray(SU2CMB_OMEGAB_H2_RANGE) / LINX_OMEGAB_H2_REFERENCE,
            ]
        )
    )
    for index, eta_fac in enumerate(eta_values, start=1):
        print(f"SU2-CMB abundance {index}/{len(eta_values)} eta_fac={eta_fac:.8g}", flush=True)
        row: dict[str, Any] = {
            "model": "SU2CMB_asymptotic",
            "eta_fac": float(eta_fac),
            "omega_b_h2": float(eta_fac * LINX_OMEGAB_H2_REFERENCE),
            "post_eann_baryon_density_factor": float(eta_fac / NU_CMB**3),
        }
        row.update(
            run_abundance(
                api,
                abundance_model,
                reaction_names,
                su2_background,
                float(eta_fac),
                cli.sampling_ntop,
                cli.rtol,
                cli.atol,
            )
        )
        score_abundances(row)
        rows.append(row)

    scan = pd.DataFrame(rows)
    summary = selected_summary(scan)
    trajectory_path = cli.outdir / "su2cmb_linx_background_trajectory.csv"
    scan_path = cli.outdir / "su2cmb_linx_eta_scan.csv"
    summary_path = cli.outdir / "su2cmb_linx_summary.csv"
    plot_path = cli.outdir / "su2cmb_linx_thermal_bridge.png"
    report_path = cli.outdir / "su2cmb_linx_thermal_bridge_report.md"
    config_path = cli.outdir / "su2cmb_linx_thermal_bridge_config.json"
    audit_path = cli.outdir / "su2cmb_linx_thermal_bridge_audit.json"

    trajectory.to_csv(trajectory_path, index=False)
    scan.to_csv(scan_path, index=False)
    summary.to_csv(summary_path, index=False)
    write_plot(trajectory, scan, plot_path)
    write_report(report_path, trajectory, summary, scan)

    config = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "python": sys.version,
        "platform": platform.platform(),
        "argv": sys.argv,
        "network": network_name,
        "sampling_nTOp": cli.sampling_ntop,
        "rtol": cli.rtol,
        "atol": cli.atol,
        "g_gamma": G_GAMMA,
        "g_YM": G_YM,
        "nu_CMB": NU_CMB,
        "extra_rho_over_rho_gamma_at_fixed_T": EXTRA_RHO_OVER_RHO_GAMMA,
        "eta_scan": eta_values.tolist(),
        "published_SU2CMB_omega_b_h2_prior": list(SU2CMB_OMEGAB_H2_RANGE),
        "nuclear_rate_pulls": "all fixed to zero",
        "nuclear_masses_binding_energies_Q_values": "standard LINX/PRIMAT inputs",
        "weak_rates": "standard LINX rates evaluated on the SU2-CMB Tnu/Tgamma trajectory",
        "sources": SOURCES,
        "identity_boundary": "SU2-CMB thermodynamic benchmark; not late-time SU2/SU2R fluid",
    }
    config_path.write_text(json.dumps(config, indent=2), encoding="utf-8")

    audit = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "programme": str(SCRIPT_PATH),
        "programme_sha256": sha256(SCRIPT_PATH),
        "cosmology_repo_commit": git_value(ROOT, "rev-parse", "HEAD"),
        "linx_repo_commit": git_value(cli.linx_path, "rev-parse", "HEAD"),
        "files": {},
    }
    for path in (trajectory_path, scan_path, summary_path, plot_path, report_path, config_path):
        audit["files"][path.name] = {"bytes": path.stat().st_size, "sha256": sha256(path)}
    audit_path.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    completion_marker.write_text(
        json.dumps(
            {
                "completed": datetime.now().isoformat(timespec="seconds"),
                "programme_sha256": sha256(SCRIPT_PATH),
                "report": str(report_path),
                "report_sha256": sha256(report_path),
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    for source in (
        trajectory_path,
        scan_path,
        summary_path,
        plot_path,
        report_path,
        config_path,
        audit_path,
        completion_marker,
    ):
        destination = cli.outputs / f"{source.stem}_{RUN_DATE}{source.suffix}"
        shutil.copy2(source, destination)
        print(f"Saved: {destination}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
