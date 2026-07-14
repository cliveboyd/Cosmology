#!/usr/bin/env python3
r"""
Build a compact SU2/electroweak field registry for the cosmology programme.

This captures the 3d matter-field transformation table supplied in the SU2
discussion:

    U(1)_Y x SU(2)_L x SU(3)_c
    Phi       1/2   2   1
    Y_tau       0   1   1
    W_tau^a     0   3   1
    G_tau^A     0   1   8

Interpretation:
    * Phi is the Higgs doublet.
    * Y_tau, W_tau^a, and G_tau^A are temporal gauge components that behave as
      scalar fields in the three-dimensional finite-temperature/effective
      theory. They are not late-time 3D matter-density fields by themselves.
    * This is a theory/metadata layer for future electroweak-plasma,
      baryogenesis, magnetic-helicity, and parity diagnostics. It does not
      alter the current SN/BAO/Planck SU2/SU2R likelihoods.

Outputs:
    plamb_runs/diagnostics/su2_field_registry/su2_3d_field_registry.csv
    plamb_runs/diagnostics/su2_field_registry/su2_3d_field_registry.json
    plamb_runs/diagnostics/su2_field_registry/su2_3d_field_registry.md
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUTDIR = ROOT / "plamb_runs" / "diagnostics" / "su2_field_registry"


@dataclass(frozen=True)
class GaugeField:
    field_id: str
    display_name: str
    u1_y: str
    su2_l: str
    su3_c: str
    representation_type: str
    is_temporal_component: bool
    three_d_role: str
    cosmology_relevance: str
    current_program_status: str
    future_hook: str
    notes: str


def build_registry() -> list[GaugeField]:
    return [
        GaugeField(
            field_id="Phi",
            display_name="Phi",
            u1_y="1/2",
            su2_l="2",
            su3_c="1",
            representation_type="Higgs fundamental doublet; color singlet",
            is_temporal_component=False,
            three_d_role="3d scalar Higgs field in the electroweak effective theory",
            cosmology_relevance="Electroweak symmetry breaking, sphalerons, finite-temperature phase structure, and possible bridge to baryogenesis diagnostics.",
            current_program_status="documented_theory_metadata",
            future_hook="diagnose_electroweak_3d_effective_fields.py",
            notes="This is the only listed field with non-zero hypercharge and SU(2)_L doublet structure.",
        ),
        GaugeField(
            field_id="Y_tau",
            display_name="Y_tau",
            u1_y="0",
            su2_l="1",
            su3_c="1",
            representation_type="hypercharge temporal gauge component; gauge singlet under the listed 3d representations",
            is_temporal_component=True,
            three_d_role="temporal U(1)_Y gauge component represented as a scalar mode after dimensional reduction",
            cosmology_relevance="Can enter finite-temperature plasma screening and chemical-potential bookkeeping, but is not a direct late-time dark-sector density parameter.",
            current_program_status="documented_theory_metadata",
            future_hook="diagnose_electroweak_temporal_components.py",
            notes="Neutral under U(1)_Y in this representation table because it is a gauge-field component rather than a charged matter multiplet.",
        ),
        GaugeField(
            field_id="W_tau_a",
            display_name="W_tau^a",
            u1_y="0",
            su2_l="3",
            su3_c="1",
            representation_type="SU(2)_L adjoint temporal gauge component",
            is_temporal_component=True,
            three_d_role="temporal weak gauge component acting as an adjoint scalar in the 3d effective theory",
            cosmology_relevance="Most directly relevant to the SU2/chiral programme; it is the electroweak temporal component that could connect finite-temperature SU(2)_L physics, sphalerons, and parity/chirality diagnostics.",
            current_program_status="documented_theory_metadata",
            future_hook="diagnose_su2_temporal_adjoint.py",
            notes="The SU(2)_L entry 3 means adjoint representation, not three ordinary spatial matter axes.",
        ),
        GaugeField(
            field_id="G_tau_A",
            display_name="G_tau^A",
            u1_y="0",
            su2_l="1",
            su3_c="8",
            representation_type="SU(3)_c adjoint temporal gauge component",
            is_temporal_component=True,
            three_d_role="temporal gluon gauge component acting as an adjoint scalar in the 3d effective theory",
            cosmology_relevance="Relevant to QCD-sector screening and early-universe plasma modelling; currently a background/theory constraint rather than an active SU2 distance-fit variable.",
            current_program_status="documented_theory_metadata",
            future_hook="diagnose_qcd_temporal_adjoint.py",
            notes="The SU(3)_c entry 8 is the color-adjoint representation.",
        ),
    ]


def write_csv(path: Path, rows: list[GaugeField]) -> None:
    dicts = [asdict(row) for row in rows]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dicts[0].keys()))
        writer.writeheader()
        writer.writerows(dicts)


def write_json(path: Path, rows: list[GaugeField]) -> None:
    path.write_text(json.dumps([asdict(row) for row in rows], indent=2), encoding="utf-8")


def write_markdown(path: Path, rows: list[GaugeField]) -> None:
    lines: list[str] = []
    lines.append("# SU2 3d field registry")
    lines.append("")
    lines.append("This registry captures the transformation table for the 3d fields under `U(1)_Y x SU(2)_L x SU(3)_c`.")
    lines.append("")
    lines.append("| field | U(1)_Y | SU(2)_L | SU(3)_c | temporal? | role |")
    lines.append("|---|---:|---:|---:|---|---|")
    for row in rows:
        temporal = "yes" if row.is_temporal_component else "no"
        lines.append(
            f"| `{row.display_name}` | {row.u1_y} | {row.su2_l} | {row.su3_c} | "
            f"{temporal} | {row.three_d_role} |"
        )
    lines.append("")
    lines.append("## Interpretation for the cosmology programme")
    lines.append("")
    lines.append("- This table belongs in the theory layer of the SU2 programme, not directly inside the late-time SN/BAO/Planck likelihood.")
    lines.append("- `Y_tau`, `W_tau^a`, and `G_tau^A` are temporal gauge components that become scalar degrees of freedom in a 3d finite-temperature/effective description.")
    lines.append("- `W_tau^a` is the most relevant entry for the SU2/chiral line because it is the SU(2)_L adjoint temporal component.")
    lines.append("- These fields can inform future modules for electroweak plasma, sphalerons, magnetic helicity, baryogenesis, and parity-sensitive observables.")
    lines.append("- They should not be read as ordinary present-day 3D matter fields unless a separate mapping from early-universe field content to late-time density/structure observables is defined.")
    lines.append("")
    lines.append("## Future hooks")
    lines.append("")
    lines.append("| field | future hook | status |")
    lines.append("|---|---|---|")
    for row in rows:
        lines.append(f"| `{row.display_name}` | `{row.future_hook}` | {row.current_program_status} |")
    lines.append("")
    lines.append("## Programme boundary")
    lines.append("")
    lines.append("Current SU2/SU2R runs remain phenomenological distance fits with parameters such as `Omega_chi0`, `su2_fraction`, `w0_chi`, and `wa_chi`. This registry documents the particle/gauge-theory context that may motivate later physical priors or diagnostics.")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the SU2/electroweak 3d field registry.")
    parser.add_argument("--outdir", type=Path, default=OUTDIR)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rows = build_registry()
    csv_path = args.outdir / "su2_3d_field_registry.csv"
    json_path = args.outdir / "su2_3d_field_registry.json"
    md_path = args.outdir / "su2_3d_field_registry.md"
    write_csv(csv_path, rows)
    write_json(json_path, rows)
    write_markdown(md_path, rows)
    print(f"Saved CSV: {csv_path}")
    print(f"Saved JSON: {json_path}")
    print(f"Saved notes: {md_path}")
    print("Most relevant temporal component for SU2 follow-up: W_tau^a")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
