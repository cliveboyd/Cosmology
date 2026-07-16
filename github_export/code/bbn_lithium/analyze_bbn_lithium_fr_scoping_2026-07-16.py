r"""Lithium-problem scoping gates for FR/no-expansion ideas.

This is a deliberately lightweight BBN emulator, not a replacement for a full
nuclear-network code such as PRIMAT, AlterBBN, PArthENoPE, or LINX.  Its role is
to make the first falsifiable question explicit:

    Can a baryon-density shift, an early clock/expansion-rate shift, or an
    Li-only nuisance move Li-7 without breaking D/H, He-4, and CMB eta?

Reference values are pinned to the 2025 PDG Big Bang Nucleosynthesis review:

* D/H x 1e6 = 25.08 +/- 0.29
* Yp(He-4) = 0.245 +/- 0.003
* Li/H = (1.45 +/- 0.25) x 1e-10, treated both as a measurement and as a
  lower-bound stress case because stellar depletion may have occurred
* eta10_CMB = 6.12 +/- 0.04

Approximate emulator choices:

* yD = 1e5 D/H is normalised to the PDG D/H concordance point and scales as
  eta_D^-1.6 over eta10 roughly 4..8, matching standard BBN fitting behaviour.
* Yp uses a local linear eta dependence plus Delta Yp ~= 0.013 Delta N_eff.
* yLi = 1e10 Li/H is normalised to the Planck/CMB standard-BBN prediction
  4.72 and scales as eta_Li^2 in the high-eta branch.
* The early clock factor S = H_BBN/H_SBBN maps to
  Delta N_eff = (43/7)(S^2 - 1).  The eta_D and eta_Li shifts follow the
  usual Steigman/Kneller-style scoping fits.

The key claim boundary is printed into the reports: this analysis can reject
simple one-knob explanations and prioritise follow-up, but any positive region
must be rerun with a full BBN network before being treated as physical.
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import sys
from dataclasses import dataclass
from datetime import datetime
from itertools import product
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[3]
DEFAULT_OUTDIR = (
    REPO_ROOT
    / "plamb_runs"
    / "diagnostics"
    / "bbn_lithium_fr_scoping_20260716"
)


@dataclass(frozen=True)
class ObservedAbundances:
    """Observed light-element and CMB-anchor values in analysis units."""

    yD: float = 2.508  # 1e5 D/H; PDG D/H x 1e6 = 25.08
    yD_sigma_obs: float = 0.029
    Yp: float = 0.245
    Yp_sigma_obs: float = 0.003
    yLi: float = 1.45  # 1e10 Li/H
    yLi_sigma_obs: float = 0.25
    eta10_cmb: float = 6.12
    eta10_cmb_sigma: float = 0.04
    omega_b_h2_cmb: float = 0.02237
    omega_b_h2_cmb_sigma: float = 0.00015


@dataclass(frozen=True)
class TheoryFloors:
    """Approximate theory/systematic floors for the emulator likelihood."""

    yD_frac: float = 0.03
    Yp_abs: float = 0.0005
    yLi_frac: float = 0.15


@dataclass(frozen=True)
class ModelSpec:
    name: str
    description: str
    vary_eta: bool
    vary_clock: bool
    vary_li_suppression: bool
    vary_stellar_depletion: bool
    eta_cmb_prior: bool
    li_lower_bound_mode: bool

    @property
    def k(self) -> int:
        return int(self.vary_eta) + int(self.vary_clock) + int(self.vary_li_suppression) + int(
            self.vary_stellar_depletion
        )


MODELS: tuple[ModelSpec, ...] = (
    ModelSpec(
        "sbbn_cmb",
        "Standard BBN emulator at CMB eta; no Li adjustment.",
        False,
        False,
        False,
        False,
        True,
        False,
    ),
    ModelSpec(
        "eta_only_free",
        "Fit eta10 only, with no CMB prior; exposes D/Li eta tension.",
        True,
        False,
        False,
        False,
        False,
        False,
    ),
    ModelSpec(
        "eta_only_with_cmb_prior",
        "Fit eta10 with the CMB eta prior retained.",
        True,
        False,
        False,
        False,
        True,
        False,
    ),
    ModelSpec(
        "clock_only_cmb_eta",
        "Fit early clock factor S at fixed CMB eta.",
        False,
        True,
        False,
        False,
        True,
        False,
    ),
    ModelSpec(
        "eta_clock_with_cmb_prior",
        "Fit eta10 and early clock factor S with CMB eta prior.",
        True,
        True,
        False,
        False,
        True,
        False,
    ),
    ModelSpec(
        "eta_clock_free",
        "Fit eta10 and early clock factor S without CMB prior.",
        True,
        True,
        False,
        False,
        False,
        False,
    ),
    ModelSpec(
        "li_production_suppression_cmb",
        "Fit a Li-only production suppression factor at fixed CMB eta.",
        False,
        False,
        True,
        False,
        True,
        False,
    ),
    ModelSpec(
        "stellar_depletion_cmb",
        "Fit post-BBN stellar depletion of surface Li at fixed CMB eta.",
        False,
        False,
        False,
        True,
        True,
        False,
    ),
    ModelSpec(
        "li_as_lower_bound_cmb",
        "Stress case: treat stellar Li as a lower bound, not a measurement.",
        False,
        False,
        False,
        False,
        True,
        True,
    ),
)


SOURCE_NOTES = {
    "pdg_2025_bbn": "https://ccwww.kek.jp/pdg/2025/reviews/rpp2025-rev-bbang-nucleosynthesis.pdf",
    "fields_lithium_review": "https://ned.ipac.caltech.edu/level5/Sept15/Fields/Fields2.html",
    "steigman_fit_context": "https://ned.ipac.caltech.edu/level5/Sept07/Steigman/Steigman2.html",
    "arxiv_recent_lithium_review": "https://arxiv.org/abs/2304.08032",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--force", action="store_true", help="Overwrite existing output files.")
    parser.add_argument("--smoke", action="store_true", help="Run a small grid for fast validation.")
    parser.add_argument("--eta-min", type=float, default=4.0)
    parser.add_argument("--eta-max", type=float, default=8.0)
    parser.add_argument("--eta-n", type=int, default=321)
    parser.add_argument("--clock-min", type=float, default=0.85)
    parser.add_argument("--clock-max", type=float, default=1.15)
    parser.add_argument("--clock-n", type=int, default=241)
    parser.add_argument("--li-supp-min", type=float, default=0.15)
    parser.add_argument("--li-supp-max", type=float, default=1.20)
    parser.add_argument("--li-supp-n", type=int, default=211)
    parser.add_argument("--depletion-min", type=float, default=1.0)
    parser.add_argument("--depletion-max", type=float, default=5.0)
    parser.add_argument("--depletion-n", type=int, default=241)
    return parser.parse_args()


def ensure_outdir(path: Path, force: bool) -> None:
    path.mkdir(parents=True, exist_ok=True)
    sentinel = path / "bbn_lithium_fr_scoping_summary.csv"
    if sentinel.exists() and not force:
        raise FileExistsError(f"{sentinel} exists; rerun with --force to overwrite.")


def grid_values(cli: argparse.Namespace) -> dict[str, np.ndarray]:
    if cli.smoke:
        return {
            "eta10": np.linspace(5.4, 6.6, 13),
            "clock": np.linspace(0.95, 1.05, 11),
            "li_suppression": np.linspace(0.25, 1.0, 16),
            "stellar_depletion": np.linspace(1.0, 4.0, 16),
        }
    return {
        "eta10": np.linspace(cli.eta_min, cli.eta_max, cli.eta_n),
        "clock": np.linspace(cli.clock_min, cli.clock_max, cli.clock_n),
        "li_suppression": np.linspace(cli.li_supp_min, cli.li_supp_max, cli.li_supp_n),
        "stellar_depletion": np.linspace(cli.depletion_min, cli.depletion_max, cli.depletion_n),
    }


def delta_neff_from_clock(clock_factor: float) -> float:
    return (43.0 / 7.0) * (clock_factor * clock_factor - 1.0)


def emulate_abundances(
    eta10: float,
    clock_factor: float,
    li_suppression: float,
    stellar_depletion: float,
) -> dict[str, float]:
    """Return BBN-emulator predictions in yD=1e5 D/H and yLi=1e10 Li/H units."""

    eta_d = eta10 - 6.0 * (clock_factor - 1.0)
    eta_li = eta10 - 3.0 * (clock_factor - 1.0)
    if eta_d <= 0 or eta_li <= 0 or li_suppression <= 0 or stellar_depletion <= 0:
        return {
            "yD": math.nan,
            "Yp": math.nan,
            "yLi_bbn": math.nan,
            "yLi_surface": math.nan,
            "delta_neff": math.nan,
        }

    delta_neff = delta_neff_from_clock(clock_factor)
    yD = 2.508 * (6.040 / eta_d) ** 1.6
    Yp = 0.2471 + 0.0016 * (eta10 - 6.040) + 0.013 * delta_neff
    yLi_bbn = 4.72 * (eta_li / 6.120) ** 2.0 * li_suppression
    yLi_surface = yLi_bbn / stellar_depletion
    return {
        "yD": float(yD),
        "Yp": float(Yp),
        "yLi_bbn": float(yLi_bbn),
        "yLi_surface": float(yLi_surface),
        "delta_neff": float(delta_neff),
    }


def chi_square(
    obs: ObservedAbundances,
    floors: TheoryFloors,
    model: ModelSpec,
    eta10: float,
    clock_factor: float,
    li_suppression: float,
    stellar_depletion: float,
) -> dict[str, float]:
    pred = emulate_abundances(eta10, clock_factor, li_suppression, stellar_depletion)
    if not all(np.isfinite(v) for v in pred.values()):
        return {"chi2_total": math.inf, **pred}

    sigma_yD = math.hypot(obs.yD_sigma_obs, floors.yD_frac * obs.yD)
    sigma_Yp = math.hypot(obs.Yp_sigma_obs, floors.Yp_abs)
    sigma_yLi = math.hypot(obs.yLi_sigma_obs, floors.yLi_frac * pred["yLi_surface"])

    chi2_yD = ((pred["yD"] - obs.yD) / sigma_yD) ** 2
    chi2_Yp = ((pred["Yp"] - obs.Yp) / sigma_Yp) ** 2
    if model.li_lower_bound_mode:
        # Treat the plateau value as a lower bound on primordial Li.  A model
        # above it is not penalised; underproduction is still disfavoured.
        residual = max(0.0, obs.yLi - pred["yLi_surface"])
        chi2_yLi = (residual / sigma_yLi) ** 2
    else:
        chi2_yLi = ((pred["yLi_surface"] - obs.yLi) / sigma_yLi) ** 2

    chi2_eta_cmb = 0.0
    if model.eta_cmb_prior:
        chi2_eta_cmb = ((eta10 - obs.eta10_cmb) / obs.eta10_cmb_sigma) ** 2

    chi2_total = chi2_yD + chi2_Yp + chi2_yLi + chi2_eta_cmb
    return {
        "chi2_total": float(chi2_total),
        "chi2_yD": float(chi2_yD),
        "chi2_Yp": float(chi2_Yp),
        "chi2_yLi": float(chi2_yLi),
        "chi2_eta_cmb": float(chi2_eta_cmb),
        "sigma_yD": float(sigma_yD),
        "sigma_Yp": float(sigma_Yp),
        "sigma_yLi": float(sigma_yLi),
        **pred,
    }


def values_for_model(model: ModelSpec, values: dict[str, np.ndarray], obs: ObservedAbundances) -> Iterable[tuple[float, float, float, float]]:
    eta_values = values["eta10"] if model.vary_eta else np.array([obs.eta10_cmb])
    clock_values = values["clock"] if model.vary_clock else np.array([1.0])
    li_values = values["li_suppression"] if model.vary_li_suppression else np.array([1.0])
    depletion_values = values["stellar_depletion"] if model.vary_stellar_depletion else np.array([1.0])
    return product(eta_values, clock_values, li_values, depletion_values)


def run_grid(
    obs: ObservedAbundances,
    floors: TheoryFloors,
    values: dict[str, np.ndarray],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary_rows: list[dict[str, object]] = []
    near_rows: list[dict[str, object]] = []

    for model in MODELS:
        best: dict[str, object] | None = None
        threshold_rows: list[dict[str, object]] = []
        n_eval = 0
        for eta10, clock_factor, li_suppression, stellar_depletion in values_for_model(model, values, obs):
            n_eval += 1
            score = chi_square(
                obs,
                floors,
                model,
                float(eta10),
                float(clock_factor),
                float(li_suppression),
                float(stellar_depletion),
            )
            row = {
                "model": model.name,
                "eta10": float(eta10),
                "clock_factor": float(clock_factor),
                "li_suppression": float(li_suppression),
                "stellar_depletion": float(stellar_depletion),
                **score,
            }
            if best is None or row["chi2_total"] < best["chi2_total"]:
                best = row
            if row["chi2_total"] < 12.0:
                threshold_rows.append(row)

        if best is None:
            raise RuntimeError(f"No grid rows evaluated for {model.name}")

        n_data = 4 if model.eta_cmb_prior else 3
        dof = max(n_data - model.k, 1)
        bic = float(best["chi2_total"] + model.k * math.log(n_data))
        summary_rows.append(
            {
                "model": model.name,
                "description": model.description,
                "k_parameters": model.k,
                "n_data_terms": n_data,
                "dof_nominal": dof,
                "n_grid_evaluated": n_eval,
                "best_chi2": float(best["chi2_total"]),
                "best_chi2_per_dof": float(best["chi2_total"]) / dof,
                "bic_like": bic,
                "best_eta10": best["eta10"],
                "best_clock_factor": best["clock_factor"],
                "best_delta_neff": best["delta_neff"],
                "best_li_suppression": best["li_suppression"],
                "best_stellar_depletion": best["stellar_depletion"],
                "pred_yD": best["yD"],
                "pred_Yp": best["Yp"],
                "pred_yLi_bbn": best["yLi_bbn"],
                "pred_yLi_surface": best["yLi_surface"],
                "chi2_yD": best["chi2_yD"],
                "chi2_Yp": best["chi2_Yp"],
                "chi2_yLi": best["chi2_yLi"],
                "chi2_eta_cmb": best["chi2_eta_cmb"],
                "n_rows_chi2_lt_12": len(threshold_rows),
            }
        )

        threshold_rows = sorted(threshold_rows, key=lambda row: row["chi2_total"])[:2000]
        near_rows.extend(threshold_rows)

    summary = pd.DataFrame(summary_rows).sort_values("best_chi2")
    near = pd.DataFrame(near_rows).sort_values(["model", "chi2_total"])
    return summary, near


def eta_inference_table(obs: ObservedAbundances) -> pd.DataFrame:
    """Back out eta proxies for each abundance under the S=1 emulator."""

    eta_from_d = 6.040 * (obs.yD / 2.508) ** (-1.0 / 1.6)
    eta_from_li = 6.120 * (obs.yLi / 4.72) ** 0.5
    eta_li_sigma = eta_from_li * 0.5 * obs.yLi_sigma_obs / obs.yLi
    eta_d_sigma = eta_from_d * (1.0 / 1.6) * obs.yD_sigma_obs / obs.yD
    rows = [
        {
            "probe": "D/H",
            "eta10_inferred": eta_from_d,
            "eta10_sigma_obs_only": eta_d_sigma,
            "note": "D/H baryometer; matches PDG BBN concordance eta by construction.",
        },
        {
            "probe": "Li/H",
            "eta10_inferred": eta_from_li,
            "eta10_sigma_obs_only": eta_li_sigma,
            "note": "Lithium as direct measurement; lower-bound interpretation weakens this tension.",
        },
        {
            "probe": "CMB",
            "eta10_inferred": obs.eta10_cmb,
            "eta10_sigma_obs_only": obs.eta10_cmb_sigma,
            "note": "Planck TT,TE,EE+lowE+lensing PDG value.",
        },
    ]
    return pd.DataFrame(rows)


def markdown_table(df: pd.DataFrame, cols: list[str], max_rows: int | None = None) -> str:
    if df.empty:
        return "_No rows._"
    out = df[cols].copy()
    if max_rows is not None:
        out = out.head(max_rows)
    return out.to_markdown(index=False, floatfmt=".6g")


def write_report(
    path: Path,
    summary: pd.DataFrame,
    eta_table: pd.DataFrame,
    obs: ObservedAbundances,
    floors: TheoryFloors,
    values: dict[str, np.ndarray],
    smoke: bool,
) -> None:
    best_sbbn = summary[summary["model"] == "sbbn_cmb"].iloc[0]
    best_clock = summary[summary["model"] == "clock_only_cmb_eta"].iloc[0]
    best_eta_clock = summary[summary["model"] == "eta_clock_with_cmb_prior"].iloc[0]
    best_deplete = summary[summary["model"] == "stellar_depletion_cmb"].iloc[0]
    best_li_supp = summary[summary["model"] == "li_production_suppression_cmb"].iloc[0]
    best_lower = summary[summary["model"] == "li_as_lower_bound_cmb"].iloc[0]

    lines: list[str] = []
    lines.append("# BBN Lithium FR Scoping Gate")
    lines.append("")
    lines.append(f"Generated: {datetime.now().isoformat(timespec='seconds')}")
    lines.append("")
    lines.append("## Claim Boundary")
    lines.append("")
    lines.append(
        "This is a first-pass emulator gate. It is suitable for ranking simple "
        "FR/no-expansion hypotheses and identifying vetoes, but any apparently "
        "viable region must be rerun with a full BBN nuclear network."
    )
    lines.append("")
    lines.append("## Input Anchors")
    lines.append("")
    lines.append(f"- D/H x 1e6 = {obs.yD * 10:.2f} +/- {obs.yD_sigma_obs * 10:.2f}")
    lines.append(f"- Yp(He-4) = {obs.Yp:.3f} +/- {obs.Yp_sigma_obs:.3f}")
    lines.append(f"- Li/H x 1e10 = {obs.yLi:.2f} +/- {obs.yLi_sigma_obs:.2f}")
    lines.append(f"- eta10_CMB = {obs.eta10_cmb:.2f} +/- {obs.eta10_cmb_sigma:.2f}")
    lines.append(
        f"- Theory floors: D {floors.yD_frac:.0%}, He {floors.Yp_abs:.4f}, "
        f"Li {floors.yLi_frac:.0%}"
    )
    lines.append("")
    lines.append("## Executive Readout")
    lines.append("")
    li_ratio = best_sbbn["pred_yLi_surface"] / obs.yLi
    lines.append(
        f"- Standard BBN at CMB eta predicts Li/H x 1e10 = "
        f"{best_sbbn['pred_yLi_surface']:.3g}, which is {li_ratio:.2f} times "
        "the plateau value used here."
    )
    lines.append(
        f"- The best clock-only FR-style shift at fixed CMB eta has "
        f"S = {best_clock['best_clock_factor']:.4g}, Delta N_eff = "
        f"{best_clock['best_delta_neff']:.4g}, and chi2 = {best_clock['best_chi2']:.3g}."
    )
    lines.append(
        f"- Allowing eta plus clock with the CMB prior gives eta10 = "
        f"{best_eta_clock['best_eta10']:.4g}, S = "
        f"{best_eta_clock['best_clock_factor']:.4g}, and chi2 = "
        f"{best_eta_clock['best_chi2']:.3g}."
    )
    lines.append(
        f"- A Li-only production suppression would need factor "
        f"{best_li_supp['best_li_suppression']:.3g}; a stellar-depletion "
        f"interpretation needs surface depletion factor "
        f"{best_deplete['best_stellar_depletion']:.3g}."
    )
    lines.append(
        f"- If Li is treated only as a lower bound, SBBN+CMB is no longer "
        f"penalised by Li; best chi2 = {best_lower['best_chi2']:.3g}. "
        "That is an astrophysical/systematic resolution, not an FR signal."
    )
    lines.append("")
    lines.append("## Model Ranking")
    lines.append("")
    lines.append(
        markdown_table(
            summary,
            [
                "model",
                "best_chi2",
                "best_eta10",
                "best_clock_factor",
                "best_delta_neff",
                "best_li_suppression",
                "best_stellar_depletion",
                "pred_yD",
                "pred_Yp",
                "pred_yLi_surface",
            ],
        )
    )
    lines.append("")
    lines.append("## Eta Inference Cross-Check")
    lines.append("")
    lines.append(
        markdown_table(
            eta_table,
            ["probe", "eta10_inferred", "eta10_sigma_obs_only", "note"],
        )
    )
    lines.append("")
    lines.append("## Interpretation")
    lines.append("")
    lines.append(
        "The scoping gate favours the conventional reading: deuterium and CMB eta "
        "pin the baryon density near eta10 ~ 6, while lithium alone prefers a much "
        "lower eta if treated as an undepleted primordial measurement. Simple "
        "FR-style clock changes do not move Li independently; they also move D/H "
        "and He-4, so deuterium is the hard veto."
    )
    lines.append("")
    lines.append(
        "The useful next experiment is therefore a full-network BBN run with "
        "explicit FR parameters mapped to dimensionless nuclear quantities: "
        "eta_BBN/eta_CMB, H_BBN/H_SBBN or equivalent clock factor, neutron-proton "
        "mass difference, deuteron binding energy, and the 3He(alpha,gamma)7Be / "
        "7Be(n,p)7Li / 7Li(p,alpha)4He reaction controls."
    )
    lines.append("")
    lines.append("## Grid")
    lines.append("")
    lines.append(f"- Smoke mode: {smoke}")
    for key, arr in values.items():
        lines.append(f"- {key}: {len(arr)} values from {arr.min():.6g} to {arr.max():.6g}")
    lines.append("")
    lines.append("## Sources")
    lines.append("")
    for key, url in SOURCE_NOTES.items():
        lines.append(f"- {key}: {url}")
    lines.append("")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_manifest(outdir: Path) -> pd.DataFrame:
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
    manifest = pd.DataFrame(rows)
    manifest.to_csv(outdir / "manifest.csv", index=False)
    return manifest


def write_config(path: Path, cli: argparse.Namespace, obs: ObservedAbundances, floors: TheoryFloors) -> None:
    config = {
        "script": str(SCRIPT_PATH),
        "repo_root": str(REPO_ROOT),
        "generated": datetime.now().isoformat(timespec="seconds"),
        "python": sys.version,
        "platform": platform.platform(),
        "argv": sys.argv,
        "smoke": bool(cli.smoke),
        "observed_abundances": obs.__dict__,
        "theory_floors": floors.__dict__,
        "sources": SOURCE_NOTES,
    }
    path.write_text(json.dumps(config, indent=2, sort_keys=True), encoding="utf-8")


def main() -> None:
    cli = parse_args()
    ensure_outdir(cli.outdir, cli.force)
    obs = ObservedAbundances()
    floors = TheoryFloors()
    values = grid_values(cli)

    summary, near = run_grid(obs, floors, values)
    eta_table = eta_inference_table(obs)

    summary.to_csv(cli.outdir / "bbn_lithium_fr_scoping_summary.csv", index=False)
    near.to_csv(cli.outdir / "bbn_lithium_fr_scoping_near_grid.csv", index=False)
    eta_table.to_csv(cli.outdir / "bbn_lithium_eta_inference.csv", index=False)
    write_report(
        cli.outdir / "bbn_lithium_fr_scoping_report.md",
        summary,
        eta_table,
        obs,
        floors,
        values,
        cli.smoke,
    )
    write_config(cli.outdir / "bbn_lithium_fr_scoping_config.json", cli, obs, floors)
    manifest = write_manifest(cli.outdir)

    print(f"Wrote {cli.outdir}")
    print(summary[["model", "best_chi2", "best_eta10", "best_clock_factor", "pred_yD", "pred_Yp", "pred_yLi_surface"]].to_string(index=False))
    print(f"Manifest rows: {len(manifest)}")


if __name__ == "__main__":
    main()
