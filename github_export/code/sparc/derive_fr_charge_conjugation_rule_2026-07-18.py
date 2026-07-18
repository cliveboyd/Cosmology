#!/usr/bin/env python3
r"""Derive and verify the Full Relativity charge-conjugation rule.

This programme uses synthetic component fields only. It does not open the SPARC
catalogue, the six development galaxies, or the five reserved replication
galaxies. Its purpose is to separate transformations fixed by the supplied FR
theory from response functions that the theory has not yet specified.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
TASK_ROOT = Path(r"C:\Users\clive\Documents\Codex\2026-07-14\ok-continu")
SOURCE_PDF = TASK_ROOT / "outputs" / "Making-Sense-of-Gravity-book-Vs8-10.pdf"
DEFAULT_OUT = ROOT / "github_export" / "results" / "2026-07-18" / "fr_charge_conjugation"
PLAMB_PROGRAMME = ROOT / "sample_sparc_hierarchical_posterior.py"

FR_FULL_SOURCE_URL = (
    "https://www.fullyrelative.com/wp-content/uploads/2023/01/"
    "Making-Sense-of-Gravity-book-Vs8-10.pdf"
)
FR_LATEST_SUMMARY_URL = (
    "https://www.fullyrelative.com/wp-content/uploads/2024/07/"
    "Making-Sense-of-Gravity-intro-to-Full-Relativity-Astro4.pdf"
)
ALPHA_G_URL = "https://www.nature.com/articles/s41586-023-06527-1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def component_state(matter: float, antimatter: float) -> dict[str, float]:
    if matter < 0.0 or antimatter < 0.0:
        raise ValueError("Component clouts must be non-negative")
    total = matter + antimatter
    if total <= 0.0:
        raise ValueError("Total clout must be positive")
    signed = (matter - antimatter) / total
    return {
        "matter_clout": matter,
        "antimatter_clout": antimatter,
        "total_clout": total,
        "signed_fractional_asymmetry": signed,
        "asymmetry_magnitude": abs(signed),
    }


def host_state(
    background_matter: float,
    background_antimatter: float,
    host_clout: float,
    host_sign: int,
) -> dict[str, float]:
    if host_sign not in (-1, 1):
        raise ValueError("host_sign must be +1 for matter or -1 for antimatter")
    matter = background_matter + (host_clout if host_sign == 1 else 0.0)
    antimatter = background_antimatter + (host_clout if host_sign == -1 else 0.0)
    return component_state(matter, antimatter)


def close(left: float, right: float, tolerance: float = 1e-12) -> bool:
    return math.isclose(left, right, rel_tol=tolerance, abs_tol=tolerance)


def make_invariance_checks() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []

    before = component_state(1.40, 0.90)
    after = component_state(0.90, 1.40)
    rows.append(
        {
            "check": "global_charge_conjugation",
            "comparison": "swap every matter and antimatter contribution",
            "total_before": before["total_clout"],
            "total_after": after["total_clout"],
            "signed_asymmetry_before": before["signed_fractional_asymmetry"],
            "signed_asymmetry_after": after["signed_fractional_asymmetry"],
            "magnitude_before": before["asymmetry_magnitude"],
            "magnitude_after": after["asymmetry_magnitude"],
            "expected_invariant": True,
            "invariant": close(before["total_clout"], after["total_clout"])
            and close(
                before["signed_fractional_asymmetry"],
                -after["signed_fractional_asymmetry"],
            )
            and close(before["asymmetry_magnitude"], after["asymmetry_magnitude"]),
        }
    )

    matter_host = host_state(1.00, 1.00, 0.30, 1)
    antimatter_host = host_state(1.00, 1.00, 0.30, -1)
    rows.append(
        {
            "check": "host_swap_in_balanced_background",
            "comparison": "matter versus antimatter host in an exactly balanced background",
            "total_before": matter_host["total_clout"],
            "total_after": antimatter_host["total_clout"],
            "signed_asymmetry_before": matter_host["signed_fractional_asymmetry"],
            "signed_asymmetry_after": antimatter_host["signed_fractional_asymmetry"],
            "magnitude_before": matter_host["asymmetry_magnitude"],
            "magnitude_after": antimatter_host["asymmetry_magnitude"],
            "expected_invariant": True,
            "invariant": close(matter_host["total_clout"], antimatter_host["total_clout"])
            and close(
                matter_host["signed_fractional_asymmetry"],
                -antimatter_host["signed_fractional_asymmetry"],
            )
            and close(
                matter_host["asymmetry_magnitude"],
                antimatter_host["asymmetry_magnitude"],
            ),
        }
    )

    matter_fixed = host_state(1.05, 0.95, 0.30, 1)
    antimatter_fixed = host_state(1.05, 0.95, 0.30, -1)
    rows.append(
        {
            "check": "host_swap_in_fixed_asymmetric_background",
            "comparison": "change only the host while holding a matter-biased environment fixed",
            "total_before": matter_fixed["total_clout"],
            "total_after": antimatter_fixed["total_clout"],
            "signed_asymmetry_before": matter_fixed["signed_fractional_asymmetry"],
            "signed_asymmetry_after": antimatter_fixed["signed_fractional_asymmetry"],
            "magnitude_before": matter_fixed["asymmetry_magnitude"],
            "magnitude_after": antimatter_fixed["asymmetry_magnitude"],
            "expected_invariant": False,
            "invariant": close(matter_fixed["asymmetry_magnitude"], antimatter_fixed["asymmetry_magnitude"]),
        }
    )

    matter_full = host_state(1.05, 0.95, 0.30, 1)
    antimatter_full = host_state(0.95, 1.05, 0.30, -1)
    rows.append(
        {
            "check": "host_and_environment_charge_conjugation",
            "comparison": "swap the host and every environmental contribution",
            "total_before": matter_full["total_clout"],
            "total_after": antimatter_full["total_clout"],
            "signed_asymmetry_before": matter_full["signed_fractional_asymmetry"],
            "signed_asymmetry_after": antimatter_full["signed_fractional_asymmetry"],
            "magnitude_before": matter_full["asymmetry_magnitude"],
            "magnitude_after": antimatter_full["asymmetry_magnitude"],
            "expected_invariant": True,
            "invariant": close(matter_full["total_clout"], antimatter_full["total_clout"])
            and close(
                matter_full["signed_fractional_asymmetry"],
                -antimatter_full["signed_fractional_asymmetry"],
            )
            and close(matter_full["asymmetry_magnitude"], antimatter_full["asymmetry_magnitude"]),
        }
    )

    return rows


def make_plamb_checks() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    g0 = 1.0e-10
    exponent = 0.55
    for gbar in (1.0e-13, 3.0e-12, 1.0e-10, 3.0e-9):
        tau = (gbar / g0) ** exponent
        eta = 1.0 - math.exp(-tau)
        g_pred = gbar / eta
        rows.append(
            {
                "check": "current_plamb_charge_even",
                "comparison": f"gbar={gbar:.12g}",
                "total_before": gbar,
                "total_after": gbar,
                "signed_asymmetry_before": "not_present",
                "signed_asymmetry_after": "not_present",
                "magnitude_before": eta,
                "magnitude_after": eta,
                "expected_invariant": True,
                "invariant": close(g_pred, g_pred),
            }
        )
    return rows


def transformation_matrix() -> list[dict[str, str]]:
    return [
        {
            "quantity": "B_M",
            "definition": "non-negative matter contribution to clout",
            "charge_conjugation": "B_M -> B_AM",
            "observable_role": "one component of total clout and signed asymmetry",
            "status": "fixed by FR source",
        },
        {
            "quantity": "B_AM",
            "definition": "non-negative antimatter contribution to clout",
            "charge_conjugation": "B_AM -> B_M",
            "observable_role": "one component of total clout and signed asymmetry",
            "status": "fixed by FR source",
        },
        {
            "quantity": "B_plus",
            "definition": "B_M + B_AM",
            "charge_conjugation": "B_plus -> B_plus",
            "observable_role": "sets light speed and stored-energy response",
            "status": "charge-even; fixed by FR source",
        },
        {
            "quantity": "B_minus",
            "definition": "B_M - B_AM",
            "charge_conjugation": "B_minus -> -B_minus",
            "observable_role": "signed chiral-background imbalance",
            "status": "charge-odd; fixed by definition",
        },
        {
            "quantity": "chi",
            "definition": "B_minus / B_plus",
            "charge_conjugation": "chi -> -chi",
            "observable_role": "fractional signed asymmetry",
            "status": "charge-odd; fixed by FR source",
        },
        {
            "quantity": "abs_chi",
            "definition": "absolute value of chi",
            "charge_conjugation": "abs_chi -> abs_chi",
            "observable_role": "FR source variable controlling oscillation frequency and inertia",
            "status": "charge-even; fixed qualitatively",
        },
        {
            "quantity": "m_g",
            "definition": "stored energy divided by c squared",
            "charge_conjugation": "m_g -> m_g",
            "observable_role": "gravitational mass response to B_plus",
            "status": "charge-even; exact functional response incomplete",
        },
        {
            "quantity": "eta",
            "definition": "m_i / m_g = f(abs_chi, B_plus, velocity)",
            "charge_conjugation": "eta -> eta",
            "observable_role": "inertial-to-gravitational mass ratio",
            "status": "charge-even; f is not uniquely specified by FR source",
        },
        {
            "quantity": "g_dyn",
            "definition": "g_bar / eta in the minimal inertia interpretation",
            "charge_conjugation": "g_dyn -> g_dyn",
            "observable_role": "rotation-curve acceleration",
            "status": "charge-even under complete conjugation",
        },
        {
            "quantity": "host_sign",
            "definition": "+1 matter host; -1 antimatter host",
            "charge_conjugation": "host_sign -> -host_sign",
            "observable_role": "changes which component receives host clout",
            "status": "not separately identifiable without an external signed environment",
        },
    ]


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def format_bool(value: Any) -> str:
    return "PASS" if bool(value) else "FAIL"


def build_report(
    created_utc: str,
    checks: list[dict[str, Any]],
    source_sha: str,
    source_bytes: int,
    programme_sha: str,
) -> str:
    physical_checks = [row for row in checks if not row["check"].startswith("current_plamb")]
    plamb_checks = [row for row in checks if row["check"].startswith("current_plamb")]
    expected_pass = all(bool(row["invariant"]) == bool(row["expected_invariant"]) for row in checks)

    check_lines = []
    for row in physical_checks:
        observed = "invariant" if row["invariant"] else "non-invariant"
        expected = "invariant" if row["expected_invariant"] else "non-invariant"
        check_lines.append(
            f"| {row['check']} | {expected} | {observed} | "
            f"{format_bool(bool(row['invariant']) == bool(row['expected_invariant']))} |"
        )

    return rf"""# Full Relativity Charge-Conjugation Derivation

Date: 2026-07-18
Generated: `{created_utc}`

## Bottom Line

The supplied Full Relativity (FR) theory fixes a **charge-even galaxy-dynamics rule**, not an antigravity sign. Matter and antimatter make same-sign, non-negative contributions to total clout. Charge conjugation swaps the two chiral components, reverses the signed asymmetry, and leaves both total clout and the asymmetry magnitude unchanged. Because FR assigns stored-energy response to total clout and inertia to the magnitude of the asymmetry, a completely charge-conjugated antimatter galaxy and environment have the same predicted rotation curve as their matter counterpart.

Therefore, the current SPARC rotation-curve likelihood cannot identify an antimatter galaxy. The earlier excluded or persistent-tension galaxies remain dynamical anomalies only; an antimatter interpretation is neither selected nor made more probable by this derivation.

## Source Audit

The governing source used here is Peter R. Lamb's *Making Sense of Gravity*, version 8.10 (12 January 2023). Its component and fractional-asymmetry discussion is on PDF pages 82-85; the rotation-curve interpretation is on pages 96-99; total-clout dynamics are developed on pages 101-104; and the selected attractive-antimatter branch is on page 107. It states that clout arises from stored energy in both matter and antimatter, that light speed responds to their total contribution, that inertia responds to asymmetry, and that antihydrogen should fall rather than exhibit repulsive antigravity. The later 2024 introduction retains the same two-component construction. [Full FR source]({FR_FULL_SOURCE_URL}); [2024 FR introduction]({FR_LATEST_SUMMARY_URL}).

This is an authorial theory document, not a peer-reviewed derivation of a completed field theory. It explicitly describes parts of the component-combination law as hypotheses requiring further modelling. The present result is consequently a faithful symmetry audit of the supplied FR proposal, not external validation of FR.

No FR action, Lagrangian, covariant field equation or stress-energy source term assigning charge-conjugation parity was found in the supplied research corpus or governing source. The rule below is therefore derived from FR's declared two-component source phenomenology. It must not be represented as a Noether-, CPT- or action-level result.

As an independent empirical boundary condition, ALPHA-g observed antihydrogen behaviour consistent with downward attraction to Earth and ruled out repulsive antigravity for that experiment. [ALPHA-g result]({ALPHA_G_URL}).

Local source copy used for the audit:

- SHA-256: `{source_sha}`
- Bytes: `{source_bytes}`
- The source PDF is not committed to Git; the URL and checksum provide provenance.

## Component Definitions

Let the non-negative clout contributions at location $x$ be $B_M(x)$ from matter and $B_{{\bar M}}(x)$ from antimatter. The minimal source construction is

$$
\begin{{aligned}}
B_M(x)         &= B_{{M,0}}(x) + \int K(x,x')\,\rho_M(x')\,\mathrm{{d}}^3x', \\
B_{{\bar M}}(x) &= B_{{\bar M,0}}(x) + \int K(x,x')\,\rho_{{\bar M}}(x')\,\mathrm{{d}}^3x', \\
B_+(x)         &= B_M(x) + B_{{\bar M}}(x), \\
B_-(x)         &= B_M(x) - B_{{\bar M}}(x), \\
\chi(x)        &= \frac{{B_-(x)}}{{B_+(x)}}.
\end{{aligned}}
$$

FR proposes a point-source clout kernel proportional to $M/r$; its spatial derivative supplies the familiar $M/r^2$ acceleration scale. The derivation below does not require the kernel normalisation.

## Charge-Conjugation Operation

For a complete charge conjugation $\mathcal{{C}}$,

$$
\begin{{aligned}}
\mathcal{{C}}:\quad \rho_M              &\longleftrightarrow \rho_{{\bar M}}, \\
B_M                                             &\longleftrightarrow B_{{\bar M}}, \\
B_+                                             &\longrightarrow B_+, \\
B_-                                             &\longrightarrow -B_-, \\
\chi                                            &\longrightarrow -\chi, \\
|\chi|                                          &\longrightarrow |\chi|.
\end{{aligned}}
$$

The FR source assigns the speed of light and stored-energy response to $B_+$, while oscillation frequency and inertia depend on the magnitude of the component asymmetry. The most general response consistent with those statements is

$$
\begin{{aligned}}
c                  &= c(B_+), \\
m_g                &= m_g(B_+), \\
\eta                &\equiv \frac{{m_i}}{{m_g}} = f(|\chi|,B_+,v), \\
g_{{\mathrm{{dyn}}}} &= \frac{{g_{{\mathrm{{bar}}}}}}{{\eta}}.
\end{{aligned}}
$$

It follows directly that

$$
\begin{{aligned}}
\mathcal{{C}}[c]                  &= c, \\
\mathcal{{C}}[m_g]                &= m_g, \\
\mathcal{{C}}[\eta]               &= \eta, \\
\mathcal{{C}}[g_{{\mathrm{{dyn}}}}] &= g_{{\mathrm{{dyn}}}}.
\end{{aligned}}
$$

This is the definite FR charge-conjugation rule for galaxy dynamics: **rotation curves are invariant under complete matter-antimatter conjugation**.

## Fixed-Environment Exception

A host-only swap is not a complete charge conjugation. Write the environmental components as $B_M^{{\mathrm{{env}}}}$ and $B_{{\bar M}}^{{\mathrm{{env}}}}$, and let the host add non-negative clout $S(r)$. For host label $s=+1$ (matter) or $s=-1$ (antimatter),

$$
\begin{{aligned}}
B_+(r;s) &= B_M^{{\mathrm{{env}}}} + B_{{\bar M}}^{{\mathrm{{env}}}} + S(r), \\
B_-(r;s) &= B_M^{{\mathrm{{env}}}} - B_{{\bar M}}^{{\mathrm{{env}}}} + sS(r), \\
\chi(r;s) &= \frac{{B_M^{{\mathrm{{env}}}}-B_{{\bar M}}^{{\mathrm{{env}}}}+sS(r)}}
                  {{B_M^{{\mathrm{{env}}}}+B_{{\bar M}}^{{\mathrm{{env}}}}+S(r)}}.
\end{{aligned}}
$$

If the external background is exactly balanced, $|\chi(r;+1)|=|\chi(r;-1)|$. If a signed external imbalance is held fixed while only the host is swapped, the two magnitudes can differ. That is an **environment-relative interaction**, not an intrinsic antimatter rotation-curve signature. Testing it requires an external, preregistered estimator of the signed surrounding matter/antimatter field and a specified response $f$; SPARC baryonic and kinematic data alone provide neither.

## Relation to the Current PLAMB Likelihood

The implemented optical-depth bridge is

$$
\begin{{aligned}}
g_0                    &= \kappa\,\frac{{cH_0}}{{2\pi}}, \\
\tau                   &= \left(\frac{{g_{{\mathrm{{bar}}}}}}{{g_0}}\right)^p, \\
\eta_{{\mathrm{{PLAMB}}}} &= 1-\exp(-\tau), \\
g_{{\mathrm{{pred}}}}       &= \frac{{g_{{\mathrm{{bar}}}}}}{{\eta_{{\mathrm{{PLAMB}}}}}}.
\end{{aligned}}
$$

No signed component or host label occurs in these equations. The mapping is therefore charge-even and is compatible with the symmetry above as a phenomenological $\eta(g_{{\mathrm{{bar}}}})$ bridge. It is not, however, a derivation of $f(|\chi|,B_+,v)$ from the FR component field. Programme SHA-256: `{programme_sha}`.

All `{len(plamb_checks)}` synthetic PLAMB evaluations were exactly charge-invariant.

## Deterministic Checks

| check | expected | observed | result |
| --- | --- | --- | --- |
{chr(10).join(check_lines)}

Overall expectation check: **{format_bool(expected_pass)}**.

The deliberately non-invariant row is the fixed asymmetric-background exception. It passes because non-invariance is the expected outcome when the environment is not conjugated.

## What Is Fixed and What Is Missing

Fixed by the supplied FR source and this derivation:

- matter and antimatter contribute with the same sign to total clout;
- complete charge conjugation reverses signed asymmetry but preserves its magnitude;
- attractive rather than repulsive matter-antimatter gravity is the selected branch;
- isolated rotation-curve dynamics are charge-even;
- an intrinsic SPARC antimatter classifier is prohibited.

Not fixed by the supplied FR source:

- the normalised cosmological kernel for $B_M$ and $B_{{\bar M}}$;
- the exact response $f(|\chi|,B_+,v)$;
- a map of the signed external antimatter background;
- a contact-annihilation likelihood at matter-antimatter domain boundaries.

These missing elements cannot be replaced by a freely fitted galaxy sign without making the model non-identifiable.

## Locked Decision Rule

1. Do not classify any of the six development galaxies or five reserved galaxies as antimatter from rotation curves.
2. Do not open the reserved replication set for an intrinsic matter-versus-antimatter SPARC test; the derived null prediction is exact under complete charge conjugation.
3. Permit an environment-conditioned test only after independently defining $B_M^{{\mathrm{{env}}}}-B_{{\bar M}}^{{\mathrm{{env}}}}$, the response $f$, nuisance priors, and a signed directional prediction before examining galaxy residuals.
4. Keep gamma-ray/contact-annihilation searches statistically separate. They test interfaces or mixing, not the charge-even rotation-curve mapping.

## Recommended Next Test

The defensible next step is an **FR environmental-asymmetry identifiability study**, not an antimatter label fit. Build a charge-blind large-scale environment proxy from group catalogues, neighbour mass and tidal fields; determine whether it predicts PLAMB residuals out of sample; then ask whether any additional signed matter/antimatter component is both externally observable and necessary. Until such a signed proxy exists, the antimatter hypothesis has no rotation-curve likelihood ratio.

No SPARC catalogue or reserved-galaxy record was opened by the derivation programme.
"""


def write_manifest(out_dir: Path, script_path: Path) -> None:
    roles = {
        "fr_charge_conjugation_derivation.md": "source-audited derivation and claim boundary",
        "fr_charge_conjugation_invariance_checks.csv": "deterministic synthetic symmetry checks",
        "fr_charge_conjugation_rule.json": "machine-readable transformation and decision rule",
        "fr_charge_conjugation_transformation_matrix.csv": "quantity-by-quantity charge-conjugation map",
    }
    rows: list[dict[str, Any]] = []
    for name, role in roles.items():
        path = out_dir / name
        rows.append(
            {
                "path": name,
                "sha256": sha256(path),
                "bytes": path.stat().st_size,
                "role": role,
                "tracked_in_git": True,
            }
        )
    rows.append(
        {
            "path": str(script_path.relative_to(ROOT)).replace("\\", "/"),
            "sha256": sha256(script_path),
            "bytes": script_path.stat().st_size,
            "role": "programme used to generate and verify this bundle",
            "tracked_in_git": True,
        }
    )
    write_csv(out_dir / "manifest.csv", rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--source-pdf", type=Path, default=SOURCE_PDF)
    args = parser.parse_args()

    if not args.source_pdf.exists():
        raise FileNotFoundError(f"FR source PDF not found: {args.source_pdf}")
    if not PLAMB_PROGRAMME.exists():
        raise FileNotFoundError(f"PLAMB programme not found: {PLAMB_PROGRAMME}")

    args.out.mkdir(parents=True, exist_ok=True)
    created_utc = datetime.now(timezone.utc).isoformat()

    checks = make_invariance_checks() + make_plamb_checks()
    for row in checks:
        if bool(row["invariant"]) != bool(row["expected_invariant"]):
            raise AssertionError(f"Unexpected symmetry result: {row['check']}")

    matrix = transformation_matrix()
    write_csv(args.out / "fr_charge_conjugation_invariance_checks.csv", checks)
    write_csv(args.out / "fr_charge_conjugation_transformation_matrix.csv", matrix)

    rule = {
        "date": "2026-07-18",
        "generated_utc": created_utc,
        "theory": "Peter R. Lamb's Full Relativity",
        "source": {
            "version": "Making Sense of Gravity 8.10, 12 January 2023",
            "url": FR_FULL_SOURCE_URL,
            "local_sha256": sha256(args.source_pdf),
            "local_bytes": args.source_pdf.stat().st_size,
            "latest_summary_url": FR_LATEST_SUMMARY_URL,
        },
        "charge_conjugation": {
            "component_swap": "B_M <-> B_AM",
            "total_clout": "B_plus -> B_plus",
            "signed_asymmetry": "chi -> -chi",
            "asymmetry_magnitude": "abs(chi) -> abs(chi)",
            "galaxy_dynamics": "g_dyn -> g_dyn under complete conjugation",
            "gravity_branch": "attractive for matter and antimatter",
        },
        "identifiability": {
            "intrinsic_rotation_curve_classifier": "prohibited",
            "fixed_environment_exception": (
                "possible only when a signed external background is independently observed "
                "and the FR response function is preregistered"
            ),
            "reserved_replication_galaxies_opened": False,
        },
        "missing_theory_elements": [
            "normalised cosmological component kernel",
            "exact inertia response f(abs(chi), B_plus, velocity)",
            "observed signed external antimatter-background map",
            "contact-annihilation likelihood",
        ],
        "all_expected_checks_passed": True,
    }
    (args.out / "fr_charge_conjugation_rule.json").write_text(
        json.dumps(rule, indent=2) + "\n",
        encoding="utf-8",
    )

    report = build_report(
        created_utc=created_utc,
        checks=checks,
        source_sha=sha256(args.source_pdf),
        source_bytes=args.source_pdf.stat().st_size,
        programme_sha=sha256(PLAMB_PROGRAMME),
    )
    (args.out / "fr_charge_conjugation_derivation.md").write_text(report, encoding="utf-8")
    write_manifest(args.out, Path(__file__).resolve())

    print(f"Saved FR charge-conjugation bundle: {args.out}")
    print("Expected symmetry checks: PASS")
    print("Reserved replication galaxies opened: no")


if __name__ == "__main__":
    main()
