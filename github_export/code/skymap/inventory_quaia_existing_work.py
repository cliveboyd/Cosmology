#!/usr/bin/env python3
"""
Inventory earlier Quaia / skymap / HHT work for the current SU2/PLAMB programme.

The older Quaia work lives partly in the project root and partly under the
legacy .spyder-py3 folder. This helper records the important source scripts,
data products, plots, and reports so they can be referenced by current reports
without rediscovering them by hand.

Outputs:
    plamb_runs/diagnostics/quaia_existing_work/quaia_existing_work_registry.csv
    plamb_runs/diagnostics/quaia_existing_work/quaia_existing_work_summary.md
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUTDIR = ROOT / "plamb_runs" / "diagnostics" / "quaia_existing_work"


@dataclass
class RegistryRow:
    category: str
    role: str
    path: Path
    relevance: str


KEY_ROWS: list[RegistryRow] = [
    RegistryRow(
        "source_link",
        "Zenodo Quaia downloader",
        ROOT / ".spyder-py3" / "download_zenodo_quaia.py",
        "Records the Quaia DOI/source link: 10.5281/zenodo.8060755.",
    ),
    RegistryRow(
        "source_link",
        "IRSA TAP table listing",
        ROOT / "irsa_tap_table_list.txt",
        "Contains an IRSA TAP entry named 'quaia'.",
    ),
    RegistryRow(
        "catalog",
        "Quaia G<20 FITS catalog",
        ROOT / "quaia_G20.0.fits",
        "Local quasar catalog used for sky-map and redshift-dipole work.",
    ),
    RegistryRow(
        "catalog",
        "Legacy Quaia G<20 FITS catalog",
        ROOT / ".spyder-py3" / "quaia_G20.0.fits",
        "Legacy copy referenced by the earlier Spyder scripts.",
    ),
    RegistryRow(
        "catalog",
        "Legacy Quaia G<20.5 FITS catalog",
        ROOT / ".spyder-py3" / "quaia_G20.5.fits",
        "Deeper legacy Quaia catalog copy, useful for later completeness checks.",
    ),
    RegistryRow(
        "conversion",
        "FITS to RA/DEC/Z CSV converter",
        ROOT / ".spyder-py3" / "quaia_to_csv.py",
        "Extracts RA, DEC, and redshift_quaia into a compact CSV.",
    ),
    RegistryRow(
        "catalog",
        "Basic RA/DEC/Z table",
        ROOT / "quaia_G20p0_basic.csv",
        "Root-level compact Quaia table for direct analysis.",
    ),
    RegistryRow(
        "skymap",
        "General FR/KdS/HHT sky-map overlay",
        ROOT / "SM_SkyMap_FR_KdS.py",
        "Builds an anisotropy sky map and overlays FR/KdS vectors and HHT points.",
    ),
    RegistryRow(
        "dipole",
        "Redshift dipole fitter",
        ROOT / ".spyder-py3" / "quaia_redshift_dipole.py",
        "Fits z_mean(n) = z0 + D.n on HEALPix pixels.",
    ),
    RegistryRow(
        "dipole",
        "Selection-aware redshift dipole fitter",
        ROOT / ".spyder-py3" / "quaia_redshift_dipole_sel.py",
        "Extends the redshift-dipole work with selection handling.",
    ),
    RegistryRow(
        "dipole",
        "B-cut dipole grid",
        ROOT / ".spyder-py3" / "skymap2" / "quaia_outputs" / "quaia_dipole_bcut_grid.txt",
        "Scans Galactic latitude cuts and redshift slices for dipole amplitude/projection.",
    ),
    RegistryRow(
        "dipole",
        "Dipole summary",
        ROOT / ".spyder-py3" / "skymap2" / "quaia_outputs" / "quaia_dipole_summary.txt",
        "Summary of full-sample and redshift-sliced Quaia redshift dipoles.",
    ),
    RegistryRow(
        "random_controls",
        "Random mock dipole statistics",
        ROOT / "quaia_G20p0_random_mocks_stats.csv",
        "Random mock control table for the Quaia sky dipole tests.",
    ),
    RegistryRow(
        "random_controls",
        "Random mock dipole/quadrupole statistics",
        ROOT / "quaia_G20p0_random_mocks_zdipole_quad.csv",
        "Random mock control table including dipole and quadrupole residual terms.",
    ),
    RegistryRow(
        "hht",
        "Build Quaia HHT input",
        ROOT / ".spyder-py3" / "quaia_build_hht_input.py",
        "Builds 1D redshift series from Quaia dipole summary for HHT.",
    ),
    RegistryRow(
        "hht",
        "Build signed A_parallel Quaia HHT input",
        ROOT / ".spyder-py3" / "quaia_build_hht_apar_signed.py",
        "Builds signed dipole-projection series suitable for chirality/asymmetry checks.",
    ),
    RegistryRow(
        "hht",
        "Quaia HHT residual analysis",
        ROOT / ".spyder-py3" / "hht_quaia_residuals.py",
        "Runs EMD/HHT, Hilbert frequencies, and surrogate controls on Quaia series.",
    ),
    RegistryRow(
        "hht",
        "B-cut 25 A_parallel HHT IMF stats",
        ROOT / ".spyder-py3" / "plamb_runs" / "hht_quaia" / "bcut25_allz_apar" / "hht_quaia_imf_stats.csv",
        "Existing HHT statistics for the bcut25 A_parallel series.",
    ),
    RegistryRow(
        "hht",
        "B-cut 25 A_parallel HHT plot",
        ROOT / ".spyder-py3" / "plamb_runs" / "hht_quaia" / "bcut25_allz_apar" / "hht_quaia.png",
        "Existing HHT diagnostic figure for the bcut25 A_parallel series.",
    ),
]


def file_status(path: Path) -> dict[str, str]:
    exists = path.exists()
    rel = str(path.relative_to(ROOT)) if path.is_absolute() and path.is_relative_to(ROOT) else str(path)
    if not exists:
        return {
            "path": rel,
            "exists": "no",
            "bytes": "",
            "last_modified": "",
        }
    stat = path.stat()
    return {
        "path": rel,
        "exists": "yes",
        "bytes": str(stat.st_size),
        "last_modified": datetime.fromtimestamp(stat.st_mtime).isoformat(timespec="seconds"),
    }


def read_head(path: Path, limit: int = 5) -> list[str]:
    if not path.exists():
        return []
    lines: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped:
                lines.append(stripped)
            if len(lines) >= limit:
                break
    return lines


def parse_full_dipole_summary() -> dict[str, str]:
    path = ROOT / ".spyder-py3" / "skymap2" / "quaia_outputs" / "quaia_dipole_summary.txt"
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if parts and parts[0] == "full" and len(parts) >= 7:
                # The label can contain spaces, for example "full sample".
                # The numeric fields are the final five columns.
                return {
                    "N_good": parts[-5],
                    "amp": parts[-4],
                    "l_deg": parts[-3],
                    "b_deg": parts[-2],
                    "a0": parts[-1],
                }
    return {}


def parse_bcut_full_rows() -> list[dict[str, str]]:
    path = ROOT / ".spyder-py3" / "skymap2" / "quaia_outputs" / "quaia_dipole_bcut_grid.txt"
    rows: list[dict[str, str]] = []
    if not path.exists():
        return rows
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 9 and parts[1] == "full":
                rows.append(
                    {
                        "bcut": parts[0],
                        "N_good": parts[2],
                        "amp": parts[3],
                        "l_deg": parts[4],
                        "b_deg": parts[5],
                        "sep_cmb": parts[6],
                        "amp_par": parts[7],
                        "frac_par": parts[8],
                    }
                )
    return rows


def parse_hht_top_energy() -> list[dict[str, str]]:
    path = ROOT / ".spyder-py3" / "plamb_runs" / "hht_quaia" / "bcut25_allz_apar" / "hht_quaia_imf_stats.csv"
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    rows.sort(key=lambda row: float(row.get("energy", "nan")), reverse=True)
    return rows[:3]


def count_matching_files(pattern: str) -> int:
    return sum(1 for _ in ROOT.glob(pattern))


def write_registry_csv(path: Path) -> None:
    fields = [
        "category",
        "role",
        "path",
        "exists",
        "bytes",
        "last_modified",
        "relevance",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for item in KEY_ROWS:
            row = {
                "category": item.category,
                "role": item.role,
                "relevance": item.relevance,
            }
            row.update(file_status(item.path))
            writer.writerow(row)


def markdown_table(rows: list[dict[str, str]], headers: list[str]) -> list[str]:
    if not rows:
        return ["No rows found."]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(row.get(h, "") for h in headers) + " |")
    return lines


def write_summary_md(path: Path) -> None:
    full = parse_full_dipole_summary()
    bcut_rows = parse_bcut_full_rows()
    hht_top = parse_hht_top_energy()
    hht_dirs = count_matching_files(".spyder-py3/plamb_runs/hht_quaia/*")
    quaia_py = count_matching_files(".spyder-py3/*quaia*.py")
    quaia_reports = count_matching_files("*Quaia*.docx") + count_matching_files("*quaia*.docx")

    lines: list[str] = []
    lines.append("# Quaia existing-work registry")
    lines.append("")
    lines.append("This note records the older Quaia / skymap analysis already present in the Cosmology workspace.")
    lines.append("")
    lines.append("## Headline")
    lines.append("")
    lines.append("- Quaia was previously used in a sky-map/redshift-dipole workflow, not only as a raw catalog.")
    lines.append("- Local data include `quaia_G20.0.fits`, compact RA/DEC/Z CSV tables, HEALPix count/mean-z/residual maps, random mocks, and HHT products.")
    lines.append("- Legacy scripts live mainly under `.spyder-py3`; some still contain Mac-style `/Users/boyde/...` paths and should be wrapped before rerunning on Windows.")
    lines.append("- The most relevant current use is as an independent angular/quasar-redshift probe for anisotropy, parity, and HHT checks, rather than as a direct replacement for scalar SN/BAO/Planck likelihoods.")
    lines.append("")
    lines.append("## Existing source/data links")
    lines.append("")
    lines.append("- Zenodo source recorded by `download_zenodo_quaia.py`: `10.5281/zenodo.8060755`.")
    lines.append("- IRSA TAP local table list includes a `quaia` table in `irsa_tap_table_list.txt`.")
    lines.append("- Root catalog files include `quaia_G20.0.fits`, `quaia_G20p0_basic.csv`, and HEALPix map products.")
    lines.append("")
    lines.append("## External references to keep with the programme")
    lines.append("")
    lines.append("- Original local downloader/source DOI: https://doi.org/10.5281/zenodo.8060755")
    lines.append("- Later Quaia catalog Zenodo record: https://zenodo.org/records/10403370")
    lines.append("- Quaia cosmology-analysis data products: https://zenodo.org/records/8098636")
    lines.append("- Quaia catalog paper: https://arxiv.org/abs/2306.17749")
    lines.append("- Quaia + CMB lensing cosmology paper: https://arxiv.org/abs/2306.17748")
    lines.append("- Recent Quaia mock-catalog paper to consider for null tests: https://arxiv.org/abs/2509.15890")
    lines.append("")
    lines.append("## Local inventory counts")
    lines.append("")
    lines.append(f"- Quaia-specific Python scripts in `.spyder-py3`: {quaia_py}")
    lines.append(f"- Quaia HHT result directories in `.spyder-py3/plamb_runs/hht_quaia`: {hht_dirs}")
    lines.append(f"- Quaia report documents in the project root: {quaia_reports}")
    lines.append("")
    lines.append("## Full-sample dipole snapshot")
    lines.append("")
    if full:
        lines.extend(
            markdown_table(
                [full],
                ["N_good", "amp", "l_deg", "b_deg", "a0"],
            )
        )
    else:
        lines.append("No full-sample dipole row found.")
    lines.append("")
    lines.append("## Full-sample b-cut scan")
    lines.append("")
    lines.extend(
        markdown_table(
            bcut_rows,
            ["bcut", "N_good", "amp", "l_deg", "b_deg", "sep_cmb", "amp_par", "frac_par"],
        )
    )
    lines.append("")
    lines.append("## Existing HHT snapshot")
    lines.append("")
    lines.append("Top energy IMFs from `bcut25_allz_apar/hht_quaia_imf_stats.csv`:")
    lines.append("")
    lines.extend(markdown_table(hht_top, ["IMF", "energy", "median_freq", "p_val"]))
    lines.append("")
    lines.append("## Relevance to current SU2/PLAMB work")
    lines.append("")
    lines.append("- SU2: Quaia can test whether the preferred SU2/SU2R branch has any sky-direction or parity-correlated signature. Current SU2 fits are scalar distance fits; Quaia would add angular structure.")
    lines.append("- PLAMB: Quaia can probe the claim that inertia/asymmetry varies by location, because a quasar redshift dipole is explicitly spatial.")
    lines.append("- HHT: Existing Quaia A_parallel and signed A_parallel series can be rerun through the current `diagnose_hht_resonance.py --mode csv` path after a small converter normalises column names.")
    lines.append("- Dataset independence: Quaia is observationally independent of SN/BAO/Planck, but it still carries catalog selection, Galactic mask, Gaia/unWISE selection, redshift-estimation, and sky-coverage assumptions.")
    lines.append("")
    lines.append("## Recommended next programme update")
    lines.append("")
    lines.append("1. Create a Windows-native Quaia wrapper that points to the local FITS/CSV files and avoids editing legacy scripts directly.")
    lines.append("2. Promote the stable outputs into `plamb_runs/diagnostics/quaia_existing_work` or a new `quaia_skymap` diagnostic folder.")
    lines.append("3. Export the existing `A_par`, signed `A_par`, and `f_par` series into CSV format for `diagnose_hht_resonance.py --mode csv`.")
    lines.append("4. Add a report section: Quaia as independent angular/quasar-redshift anisotropy probe for SU2 chirality and PLAMB location-dependent inertia.")
    lines.append("5. Do not rerun SU2/PLAMB core chains solely because of Quaia until a formal likelihood hook exists; use Quaia first as a diagnostic/validation branch.")
    lines.append("")
    lines.append("## Key files")
    lines.append("")
    lines.append("See `quaia_existing_work_registry.csv` for file existence, sizes, and timestamps.")
    lines.append("")

    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    csv_path = OUTDIR / "quaia_existing_work_registry.csv"
    md_path = OUTDIR / "quaia_existing_work_summary.md"
    write_registry_csv(csv_path)
    write_summary_md(md_path)
    print(f"Saved registry: {csv_path}")
    print(f"Saved summary: {md_path}")


if __name__ == "__main__":
    main()
