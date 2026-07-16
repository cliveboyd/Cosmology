r"""Full-network LINX BBN scan for FR/no-expansion Lithium gates.

This script wraps the upstream LINX BBN code cloned under
``plamb_runs/deps/LINX``.  It is intended as the full-network successor to the
lightweight Lithium scoping emulator:

* eta_BBN / eta_CMB is mapped to LINX ``eta_fac``;
* an early-clock/Hubble-rate proxy is mapped to initial Delta N_eff through
  S = H_BBN/H_SBBN and Delta N_eff = (43/7)(S^2 - 1);
* neutron/weak-sector sensitivity is mapped to ``tau_n_fac``;
* Li/Be nuclear-rate controls are mapped to selected LINX ``nuclear_rates_q``
  entries, with Gaussian pull penalties.

Important limitation:
LINX exposes rate and weak-sector nuisance controls, but it does not expose a
fully self-consistent public API for varying the neutron-proton mass difference
or deuteron binding energy through every weak rate, NSE abundance, Q-value, and
tabulated nuclear rate.  This script therefore records tau/rate proxy gates
rather than pretending to have performed a true fundamental-constant scan.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import platform
import sys
import time
from dataclasses import dataclass
from datetime import UTC, datetime
from itertools import product
from pathlib import Path

import numpy as np


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[3]
DEFAULT_LINX_PATH = REPO_ROOT / "plamb_runs" / "deps" / "LINX"
DEFAULT_OUTDIR = REPO_ROOT / "plamb_runs" / "diagnostics" / "bbn_lithium_linx_fr_network_20260716"


@dataclass(frozen=True)
class ObservedAbundances:
    """Observed anchors in the same units used by the output table."""

    yD: float = 2.508  # 1e5 D/H; PDG D/H x 1e6 = 25.08
    yD_sigma: float = 0.029
    Yp: float = 0.245
    Yp_sigma: float = 0.003
    yLi: float = 1.45  # 1e10 Li/H
    yLi_sigma: float = 0.25
    eta_fac_cmb: float = 1.0
    eta_fac_sigma: float = 0.04 / 6.12


@dataclass(frozen=True)
class ScanPoint:
    point_id: str
    family: str
    eta_fac: float
    clock_factor: float
    tau_n_fac: float
    q_values: dict[str, float]


SOURCE_NOTES = {
    "linx_repository": "https://github.com/cgiovanetti/LINX",
    "linx_docs": "https://linx.readthedocs.io/",
    "pdg_2025_bbn": "https://ccwww.kek.jp/pdg/2025/reviews/rpp2025-rev-bbang-nucleosynthesis.pdf",
    "fields_lithium_review": "https://ned.ipac.caltech.edu/level5/Sept15/Fields/Fields2.html",
}

TARGET_REACTIONS = (
    "He3aBe7g",  # 3He(alpha,gamma)7Be, dominant high-eta production path
    "Be7nLi7p",  # 7Be(n,p)7Li conversion
    "Li7paa",  # 7Li(p,alpha)alpha destruction
    "Li7paag",
    "Be7naa",
    "Be7daap",
    "Be7pB8g",
    "Li7daan",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--linx-path", type=Path, default=DEFAULT_LINX_PATH)
    parser.add_argument("--network", choices=["key", "small", "full"], default="full")
    parser.add_argument("--preset", choices=["smoke", "overnight"], default="overnight")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--max-points", type=int, default=None)
    parser.add_argument("--n-random", type=int, default=900)
    parser.add_argument("--seed", type=int, default=20260716)
    parser.add_argument("--rtol", type=float, default=1e-5)
    parser.add_argument("--atol", type=float, default=1e-9)
    parser.add_argument("--sampling-key", type=int, default=120)
    parser.add_argument("--sampling-small", type=int, default=90)
    parser.add_argument("--sampling-full", type=int, default=45)
    parser.add_argument("--summary-every", type=int, default=25)
    return parser.parse_args()


def import_linx(linx_path: Path):
    if not linx_path.exists():
        raise FileNotFoundError(f"LINX source path not found: {linx_path}")
    sys.path.insert(0, str(linx_path))
    import jax  # type: ignore
    import jax.numpy as jnp  # type: ignore

    jax.config.update("jax_enable_x64", True)

    from linx.background import BackgroundModel  # type: ignore
    from linx.abundances import AbundanceModel  # type: ignore
    from linx.nuclear import NuclearRates  # type: ignore
    import linx.const as const  # type: ignore

    return jax, jnp, BackgroundModel, AbundanceModel, NuclearRates, const


def clock_to_delta_neff(clock_factor: float) -> float:
    return (43.0 / 7.0) * (clock_factor * clock_factor - 1.0)


def q_label(q_values: dict[str, float]) -> str:
    if not q_values:
        return "q0"
    parts = [f"{name}{value:+.2f}" for name, value in sorted(q_values.items()) if abs(value) > 1e-12]
    return "_".join(parts) if parts else "q0"


def make_point_id(family: str, eta_fac: float, clock_factor: float, tau_n_fac: float, q_values: dict[str, float]) -> str:
    return f"{family}|eta{eta_fac:.6f}|S{clock_factor:.6f}|tau{tau_n_fac:.7f}|{q_label(q_values)}"


def add_point(points: list[ScanPoint], family: str, eta_fac: float, clock_factor: float, tau_n_fac: float, q_values: dict[str, float]) -> None:
    point_id = make_point_id(family, eta_fac, clock_factor, tau_n_fac, q_values)
    points.append(
        ScanPoint(
            point_id=point_id,
            family=family,
            eta_fac=float(eta_fac),
            clock_factor=float(clock_factor),
            tau_n_fac=float(tau_n_fac),
            q_values={k: float(v) for k, v in q_values.items() if abs(float(v)) > 1e-12},
        )
    )


def build_points(preset: str, n_random: int, seed: int) -> list[ScanPoint]:
    tau_sigma = 0.6 / 879.4
    points: list[ScanPoint] = []

    if preset == "smoke":
        add_point(points, "baseline", 1.0, 1.0, 1.0, {})
        add_point(points, "eta_low", 0.96, 1.0, 1.0, {})
        add_point(points, "clock_high", 1.0, 1.03, 1.0, {})
        add_point(points, "tau_high", 1.0, 1.0, 1.0 + 2 * tau_sigma, {})
        add_point(points, "li_prod_low", 1.0, 1.0, 1.0, {"He3aBe7g": -2.0})
        add_point(points, "li_dest_high", 1.0, 1.0, 1.0, {"Li7paa": 2.0})
        return points

    add_point(points, "baseline", 1.0, 1.0, 1.0, {})

    eta_grid = np.linspace(0.88, 1.08, 11)
    clock_grid = np.linspace(0.92, 1.08, 9)
    tau_grid = np.array([1.0 - 2.0 * tau_sigma, 1.0, 1.0 + 2.0 * tau_sigma])
    for eta_fac, clock_factor, tau_n_fac in product(eta_grid, clock_grid, tau_grid):
        add_point(points, "eta_clock_tau_grid", eta_fac, clock_factor, tau_n_fac, {})

    for reaction in TARGET_REACTIONS:
        for q in [-4.0, -3.0, -2.0, -1.0, 1.0, 2.0, 3.0, 4.0]:
            add_point(points, "individual_li_be_rate_pull", 1.0, 1.0, 1.0, {reaction: q})

    for q_prod, q_conv, q_dest in product([-5.0, -4.0, -3.0, -2.0, 0.0], [0.0, 1.5, 3.0], [0.0, 2.0, 4.0]):
        add_point(
            points,
            "li_be_combo_ladder",
            1.0,
            1.0,
            1.0,
            {"He3aBe7g": q_prod, "Be7nLi7p": q_conv, "Li7paa": q_dest},
        )

    for eta_fac, clock_factor, q_prod, q_dest in product(
        [0.94, 0.98, 1.0, 1.02, 1.06],
        [0.96, 1.0, 1.04],
        [-4.0, -2.0, 0.0],
        [0.0, 2.0, 4.0],
    ):
        add_point(
            points,
            "fr_li_joint_ladder",
            eta_fac,
            clock_factor,
            1.0,
            {"He3aBe7g": q_prod, "Li7paa": q_dest},
        )

    rng = np.random.default_rng(seed)
    for i in range(n_random):
        q_values = {}
        for reaction in TARGET_REACTIONS:
            value = float(np.clip(rng.normal(0.0, 1.35), -4.0, 4.0))
            if abs(value) > 0.25:
                q_values[reaction] = value
        add_point(
            points,
            "random_targeted_li_be",
            float(rng.uniform(0.90, 1.08)),
            float(rng.uniform(0.93, 1.07)),
            float(1.0 + np.clip(rng.normal(0.0, tau_sigma), -3.0 * tau_sigma, 3.0 * tau_sigma)),
            q_values,
        )

    dedup: dict[str, ScanPoint] = {}
    for point in points:
        dedup[point.point_id] = point
    return list(dedup.values())


def normalise_points_for_reactions(points: list[ScanPoint], reaction_names: list[str]) -> list[ScanPoint]:
    """Drop q pulls unavailable in the selected network and deduplicate rows."""

    available = set(reaction_names)
    dedup: dict[str, ScanPoint] = {}
    for point in points:
        q_values = {name: value for name, value in point.q_values.items() if name in available}
        point_id = make_point_id(point.family, point.eta_fac, point.clock_factor, point.tau_n_fac, q_values)
        dedup[point_id] = ScanPoint(
            point_id=point_id,
            family=point.family,
            eta_fac=point.eta_fac,
            clock_factor=point.clock_factor,
            tau_n_fac=point.tau_n_fac,
            q_values=q_values,
        )
    return list(dedup.values())


def output_fieldnames(reaction_names: list[str]) -> list[str]:
    return [
        "point_id",
        "family",
        "status",
        "message",
        "elapsed_s",
        "eta_fac",
        "omega_b_h2",
        "eta10_approx",
        "clock_factor",
        "delta_neff_initial",
        "neff_final",
        "tau_n_fac",
        "q_pull_chi2",
        "n_nonzero_q",
        "Yp_mass",
        "D_H_1e5",
        "He3_H_1e5",
        "Li7_H_1e10",
        "Li6_H_1e14",
        "stellar_depletion_needed",
        "chi2_D",
        "chi2_He4",
        "chi2_Li_measurement",
        "chi2_Li_lower_bound",
        "chi2_eta_cmb",
        "chi2_tau",
        "chi2_total_li_measurement",
        "chi2_total_li_lower_bound",
        "run_utc",
    ] + [f"q_{name}" for name in reaction_names]


def existing_done_ids(path: Path) -> set[str]:
    if not path.exists():
        return set()
    done: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("status") == "ok":
                done.add(row["point_id"])
    return done


def append_row(path: Path, fieldnames: list[str], row: dict[str, object]) -> None:
    exists = path.exists()
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def as_float(row: dict[str, str], key: str, default: float = math.nan) -> float:
    try:
        return float(row.get(key, default))
    except (TypeError, ValueError):
        return default


def write_interim_summary(outdir: Path, rows: list[dict[str, str]], observed: ObservedAbundances, completed: int, total: int) -> None:
    ok = [row for row in rows if row.get("status") == "ok"]
    summary_path = outdir / "bbn_lithium_linx_fr_network_interim_summary.csv"
    report_path = outdir / "bbn_lithium_linx_fr_network_interim_report.md"
    if not ok:
        return

    ranking_measurement = sorted(ok, key=lambda row: as_float(row, "chi2_total_li_measurement"))
    ranking_lower = sorted(ok, key=lambda row: as_float(row, "chi2_total_li_lower_bound"))
    families = {}
    for row in ok:
        fam = row["family"]
        families.setdefault(fam, {"n": 0, "best_measurement": math.inf, "best_lower": math.inf})
        families[fam]["n"] += 1
        families[fam]["best_measurement"] = min(families[fam]["best_measurement"], as_float(row, "chi2_total_li_measurement"))
        families[fam]["best_lower"] = min(families[fam]["best_lower"], as_float(row, "chi2_total_li_lower_bound"))

    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = ["family", "n_ok", "best_chi2_li_measurement", "best_chi2_li_lower_bound"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for family, stats in sorted(families.items()):
            writer.writerow(
                {
                    "family": family,
                    "n_ok": stats["n"],
                    "best_chi2_li_measurement": stats["best_measurement"],
                    "best_chi2_li_lower_bound": stats["best_lower"],
                }
            )

    def row_line(row: dict[str, str]) -> str:
        return (
            f"| {row['family']} | {as_float(row, 'chi2_total_li_measurement'):.4g} | "
            f"{as_float(row, 'chi2_total_li_lower_bound'):.4g} | {as_float(row, 'eta_fac'):.5g} | "
            f"{as_float(row, 'clock_factor'):.5g} | {as_float(row, 'tau_n_fac'):.7g} | "
            f"{as_float(row, 'D_H_1e5'):.5g} | {as_float(row, 'Yp_mass'):.5g} | "
            f"{as_float(row, 'Li7_H_1e10'):.5g} | {as_float(row, 'stellar_depletion_needed'):.5g} |"
        )

    lines = [
        "# LINX BBN Lithium FR Scan",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        f"Progress: {completed} / {total} requested points completed.",
        "",
        "## Current Best Rows, Li As Measurement",
        "",
        "| family | chi2_meas | chi2_lower | eta_fac | S | tau_n_fac | D/H x1e5 | Yp | Li7/H x1e10 | Li depletion |",
        "|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|",
    ]
    lines.extend(row_line(row) for row in ranking_measurement[:12])
    lines.extend(
        [
            "",
            "## Current Best Rows, Li As Lower Bound",
            "",
            "| family | chi2_meas | chi2_lower | eta_fac | S | tau_n_fac | D/H x1e5 | Yp | Li7/H x1e10 | Li depletion |",
            "|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|",
        ]
    )
    lines.extend(row_line(row) for row in ranking_lower[:12])
    lines.extend(
        [
            "",
            "## Claim Boundary",
            "",
            "This is a selected-network LINX scan of eta, early-clock/Delta-Neff, neutron-lifetime, and selected Li/Be rate-pull proxies. It is not a self-consistent fundamental-constant variation scan.",
            "",
            f"Observed anchors: D/H x1e5={observed.yD}, Yp={observed.Yp}, Li/H x1e10={observed.yLi}.",
        ]
    )
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


class LinxRunner:
    def __init__(
        self,
        linx_path: Path,
        network: str,
        rtol: float,
        atol: float,
        sampling_key: int,
        sampling_small: int,
        sampling_full: int,
    ) -> None:
        (
            self.jax,
            self.jnp,
            self.BackgroundModel,
            self.AbundanceModel,
            self.NuclearRates,
            self.const,
        ) = import_linx(linx_path)
        self.network = network
        self.rtol = rtol
        self.atol = atol
        self.sampling_key = sampling_key
        self.sampling_small = sampling_small
        self.sampling_full = sampling_full
        self.background_model = self.BackgroundModel()
        if network == "key":
            self.key_rates = self.NuclearRates(nuclear_net="key_PRIMAT_2023")
            self.key_model = self.AbundanceModel(self.key_rates, throw=False)
            self.reaction_names = list(self.key_rates.reactions_names)
        elif network == "small":
            self.small_rates = self.NuclearRates(nuclear_net="small_PRIMAT_2023")
            self.small_model = self.AbundanceModel(self.small_rates, throw=False)
            self.reaction_names = list(self.small_rates.reactions_names)
        else:
            self.small_rates = self.NuclearRates(nuclear_net="small_PRIMAT_2023")
            self.full_rates = self.NuclearRates(nuclear_net="full_PRIMAT_2023")
            self.small_model = self.AbundanceModel(self.small_rates, throw=False)
            self.full_model = self.AbundanceModel(self.full_rates, throw=False)
            self.reaction_names = list(self.full_rates.reactions_names)
        self.background_cache: dict[float, tuple[object, ...]] = {}

    def background(self, delta_neff: float) -> tuple[object, ...]:
        key = round(float(delta_neff), 8)
        if key not in self.background_cache:
            self.background_cache[key] = self.background_model(self.jnp.asarray(delta_neff))
        return self.background_cache[key]

    def q_vector(self, point: ScanPoint, names: list[str]) -> object:
        arr = np.zeros(len(names), dtype=float)
        for name, value in point.q_values.items():
            if name in names:
                arr[names.index(name)] = value
        return self.jnp.asarray(arr)

    def run_point(self, point: ScanPoint) -> dict[str, float]:
        jnp = self.jnp
        delta_neff = clock_to_delta_neff(point.clock_factor)
        t_vec, a_vec, rho_g, rho_nu, rho_np, p_np, neff_vec = self.background(delta_neff)
        if self.network in {"key", "small"}:
            q = self.q_vector(point, self.reaction_names)
            model = self.key_model if self.network == "key" else self.small_model
            sampling = self.sampling_key if self.network == "key" else self.sampling_small
            y = model(
                jnp.array(rho_g),
                jnp.array(rho_nu),
                jnp.array(rho_np),
                jnp.array(p_np),
                t_vec=jnp.array(t_vec),
                a_vec=jnp.array(a_vec),
                eta_fac=jnp.asarray(point.eta_fac),
                tau_n_fac=jnp.asarray(point.tau_n_fac),
                nuclear_rates_q=q,
                sampling_nTOp=sampling,
                rtol=self.rtol,
                atol=self.atol,
            )
        else:
            q_small = self.q_vector(point, list(self.small_rates.reactions_names))
            q_full = self.q_vector(point, self.reaction_names)
            y_switch = self.small_model(
                jnp.array(rho_g),
                jnp.array(rho_nu),
                jnp.array(rho_np),
                jnp.array(p_np),
                t_vec=jnp.array(t_vec),
                a_vec=jnp.array(a_vec),
                eta_fac=jnp.asarray(point.eta_fac),
                tau_n_fac=jnp.asarray(point.tau_n_fac),
                nuclear_rates_q=q_small,
                T_end=self.const.T_switch,
                sampling_nTOp=self.sampling_small,
                rtol=self.rtol,
                atol=self.atol,
            )
            pad = tuple(0.0 for _ in range(self.full_rates.max_i_species - len(tuple(y_switch))))
            y = self.full_model(
                jnp.array(rho_g),
                jnp.array(rho_nu),
                jnp.array(rho_np),
                jnp.array(p_np),
                t_vec=jnp.array(t_vec),
                a_vec=jnp.array(a_vec),
                eta_fac=jnp.asarray(point.eta_fac),
                tau_n_fac=jnp.asarray(point.tau_n_fac),
                nuclear_rates_q=q_full,
                Y_i=tuple(y_switch) + pad,
                T_start=self.const.T_switch,
                sampling_nTOp=self.sampling_full,
                rtol=self.rtol,
                atol=self.atol,
            )

        y_np = np.array([float(v) for v in tuple(y)], dtype=float)
        y_p = y_np[1]
        y_d = y_np[2]
        y_t = y_np[3] if len(y_np) > 3 else 0.0
        y_he3 = y_np[4] if len(y_np) > 4 else 0.0
        y_he4 = y_np[5] if len(y_np) > 5 else 0.0
        y_li7 = y_np[6] if len(y_np) > 6 else 0.0
        y_be7 = y_np[7] if len(y_np) > 7 else 0.0
        y_li6 = y_np[8] if len(y_np) > 8 else 0.0

        return {
            "neff_final": float(neff_vec[-1]),
            "Yp_mass": 4.0 * y_he4,
            "D_H_1e5": y_d / y_p * 1e5,
            "He3_H_1e5": (y_t + y_he3) / y_p * 1e5,
            "Li7_H_1e10": (y_li7 + y_be7) / y_p * 1e10,
            "Li6_H_1e14": y_li6 / y_p * 1e14,
        }


def score_row(row: dict[str, object], point: ScanPoint, observed: ObservedAbundances, tau_sigma: float) -> None:
    pred_yD = float(row["D_H_1e5"])
    pred_Yp = float(row["Yp_mass"])
    pred_yLi = float(row["Li7_H_1e10"])
    q_pull_chi2 = sum(value * value for value in point.q_values.values())
    chi2_D = ((pred_yD - observed.yD) / observed.yD_sigma) ** 2
    chi2_He4 = ((pred_Yp - observed.Yp) / observed.Yp_sigma) ** 2
    chi2_Li_measurement = ((pred_yLi - observed.yLi) / observed.yLi_sigma) ** 2
    chi2_Li_lower_bound = (max(0.0, observed.yLi - pred_yLi) / observed.yLi_sigma) ** 2
    chi2_eta = ((point.eta_fac - observed.eta_fac_cmb) / observed.eta_fac_sigma) ** 2
    chi2_tau = ((point.tau_n_fac - 1.0) / tau_sigma) ** 2
    prior = q_pull_chi2 + chi2_eta + chi2_tau
    row.update(
        {
            "q_pull_chi2": q_pull_chi2,
            "n_nonzero_q": len(point.q_values),
            "stellar_depletion_needed": pred_yLi / observed.yLi if observed.yLi > 0 else math.nan,
            "chi2_D": chi2_D,
            "chi2_He4": chi2_He4,
            "chi2_Li_measurement": chi2_Li_measurement,
            "chi2_Li_lower_bound": chi2_Li_lower_bound,
            "chi2_eta_cmb": chi2_eta,
            "chi2_tau": chi2_tau,
            "chi2_total_li_measurement": chi2_D + chi2_He4 + chi2_Li_measurement + prior,
            "chi2_total_li_lower_bound": chi2_D + chi2_He4 + chi2_Li_lower_bound + prior,
        }
    )


def write_config(outdir: Path, cli: argparse.Namespace, points: list[ScanPoint], runner: LinxRunner) -> None:
    config = {
        "script": str(SCRIPT_PATH),
        "repo_root": str(REPO_ROOT),
        "linx_path": str(cli.linx_path),
        "generated": datetime.now().isoformat(timespec="seconds"),
        "python": sys.version,
        "platform": platform.platform(),
        "argv": sys.argv,
        "network": cli.network,
        "preset": cli.preset,
        "n_points_requested": len(points),
        "target_reactions_requested": TARGET_REACTIONS,
        "reaction_names_available": runner.reaction_names,
        "sources": SOURCE_NOTES,
        "claim_boundary": (
            "Full-network LINX rate/clock/tau proxy scan; not a complete "
            "fundamental-constant variation scan for np mass difference or deuteron binding energy."
        ),
    }
    (outdir / "bbn_lithium_linx_fr_network_config.json").write_text(
        json.dumps(config, indent=2, sort_keys=True),
        encoding="utf-8",
    )


def write_manifest(outdir: Path) -> None:
    rows = []
    for path in sorted(outdir.glob("*")):
        if path.is_file():
            rows.append(
                {
                    "file": path.name,
                    "bytes": path.stat().st_size,
                    "last_modified": datetime.fromtimestamp(path.stat().st_mtime).isoformat(timespec="seconds"),
                }
            )
    with (outdir / "manifest.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["file", "bytes", "last_modified"])
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    cli = parse_args()
    cli.outdir.mkdir(parents=True, exist_ok=True)
    rows_path = cli.outdir / "bbn_lithium_linx_fr_network_points.csv"
    if rows_path.exists() and cli.force and not cli.resume:
        rows_path.unlink()
    if rows_path.exists() and not cli.resume and not cli.force:
        raise FileExistsError(f"{rows_path} exists; pass --resume or --force.")

    observed = ObservedAbundances()
    runner = LinxRunner(
        cli.linx_path,
        cli.network,
        cli.rtol,
        cli.atol,
        cli.sampling_key,
        cli.sampling_small,
        cli.sampling_full,
    )
    points = normalise_points_for_reactions(build_points(cli.preset, cli.n_random, cli.seed), runner.reaction_names)
    if cli.max_points is not None:
        points = points[: cli.max_points]
    write_config(cli.outdir, cli, points, runner)
    fieldnames = output_fieldnames(runner.reaction_names)
    done = existing_done_ids(rows_path) if cli.resume else set()
    tau_sigma = 0.6 / 879.4

    print(
        f"[{datetime.now().isoformat(timespec='seconds')}] Starting LINX BBN scan: "
        f"network={cli.network}, preset={cli.preset}, points={len(points)}, already_done={len(done)}",
        flush=True,
    )
    start_all = time.time()
    completed = len(done)

    for index, point in enumerate(points, start=1):
        if point.point_id in done:
            continue
        print(
            f"[{datetime.now().isoformat(timespec='seconds')}] Point {index}/{len(points)} "
            f"{point.family} eta={point.eta_fac:.5f} S={point.clock_factor:.5f} "
            f"tau={point.tau_n_fac:.7f} q={q_label(point.q_values)}",
            flush=True,
        )
        start = time.time()
        base_row: dict[str, object] = {
            "point_id": point.point_id,
            "family": point.family,
            "eta_fac": point.eta_fac,
            "omega_b_h2": point.eta_fac * 0.02242,
            "eta10_approx": point.eta_fac * 6.12,
            "clock_factor": point.clock_factor,
            "delta_neff_initial": clock_to_delta_neff(point.clock_factor),
            "tau_n_fac": point.tau_n_fac,
            "run_utc": datetime.now(UTC).isoformat(timespec="seconds").replace("+00:00", "Z"),
        }
        for name in runner.reaction_names:
            base_row[f"q_{name}"] = point.q_values.get(name, 0.0)
        try:
            pred = runner.run_point(point)
            row = {**base_row, **pred, "status": "ok", "message": ""}
            score_row(row, point, observed, tau_sigma)
        except Exception as exc:  # pragma: no cover - long-scan resilience
            row = {
                **base_row,
                "status": "failed",
                "message": repr(exc),
                "q_pull_chi2": sum(value * value for value in point.q_values.values()),
                "n_nonzero_q": len(point.q_values),
            }
        row["elapsed_s"] = time.time() - start
        append_row(rows_path, fieldnames, row)
        completed += 1
        if completed % cli.summary_every == 0 or completed == len(points):
            rows = read_rows(rows_path)
            write_interim_summary(cli.outdir, rows, observed, completed, len(points))
            write_manifest(cli.outdir)
            print(
                f"[{datetime.now().isoformat(timespec='seconds')}] Interim summary written "
                f"after {completed}/{len(points)} points; elapsed={(time.time() - start_all) / 3600:.3f} h",
                flush=True,
            )

    rows = read_rows(rows_path)
    write_interim_summary(cli.outdir, rows, observed, completed, len(points))
    write_manifest(cli.outdir)
    print(
        f"[{datetime.now().isoformat(timespec='seconds')}] Complete: {completed}/{len(points)} "
        f"points in {(time.time() - start_all) / 3600:.3f} h. Output: {cli.outdir}",
        flush=True,
    )


if __name__ == "__main__":
    main()
