from __future__ import annotations

import csv
import datetime as dt
import json
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(r"C:\Users\clive\Documents\Cosmology")
THREAD_ROOT = Path(r"C:\Users\clive\Documents\Codex\2026-07-14\ok-continu")
EXPORT = ROOT / "github_export"
DATE = "2026-07-14"

SOFFICE = Path(r"C:\Program Files\LibreOffice\program\soffice.com")
PDF_ROOT = EXPORT / "docs" / "rendered_pdfs"
MANIFEST_DIR = EXPORT / "manifests"
OUTPUT_DIR = THREAD_ROOT / "outputs"
LO_PROFILE = THREAD_ROOT / "work" / "lo_profile_research_pdf"
TIMEOUT_SECONDS = 180


SKYMAP_RE = re.compile(r"(skymap|sky.?map|quaia|zdipole|rbao|dipole|quadrupole)", re.I)

EXTRA_STREAMS: dict[str, list[str]] = {
    "plamb_core_sn_fr": [
        "PLamb_results_summary_report.docx",
        "PLamb_Master_Summary.docx",
        "PLamb_Update_Summary.docx",
        "PLamb_Cosmology_Summary.docx",
        "PLamb_Cosmology_Summary_with_Acronyms_and_Script.docx",
        "PLamb_Test_10V6c_plus_Program_Guide_v1.1.docx",
        "PLamb_Cosmology_Report_V10V3_Updated.docx",
        "PLamb_Cosmology_Report_ΛCDM V10V3_Updated.docx",
        "PLamb_FR_NoDE_Summary_Report.docx",
        "PLamb_FR_Results_V10V3.docx",
        "PLamb_Report_FR_Master_Updated.docx",
        "PLamb_SN_Summary_KdS_LCDM_FR_V10V6c.docx",
        "PLamb_Consolidated_Results_20250826_1051.docx",
        "Cosmology_BIC_Report_Aug31_with_residuals.docx",
        "Pantheon_Fitting_Report.docx",
        "Pantheon_FR_Run_Report.docx",
        "DESI_DR2_Lya_Dataset_Status_Report.docx",
    ],
    "kds_bh_stf": [
        "KdS_Cosmology_Results_Summary.docx",
        "KdS_Fit_Expanded_RealData_Report_TextOnly.docx",
        "KdS_Fit_Results_Real_PantheonPlusSHOES.docx",
        "KdS_MCMC_Report.docx",
        "KdS_rSTF_Analysis_Report.docx",
        "KdS_rSTF_Fit_Report_v2.docx",
        "Kerr_deSitter_Model_CMB.docx",
        "Kerr_deSitter_Model_Expanded_CliveBoyd_Final_Updated.docx",
        "STF_BAO_Fit_Summary.docx",
        "STF_Model_Summary_Report_Updated.docx",
        "STF_vs_LCDM_Residuals_Report.docx",
        "Black Hole Cosmology - Mass Accretion Detection..docx",
        "Black Hole Cosmology - SHOES_Pantheon_BH_Evolution_Report.docx",
        "SHOES_Pantheon_BH_Evolution_Report.docx",
        "Bayesian_Hierarchical_Model_Results_2025-07-09_20-20.docx",
    ],
    "su2_void_gw_jwst": [
        "SU2_Chiral_Gauge_Sector_Dark_Energy_Proposal.docx",
        "Void Analysis.docx",
        "Void_Analysis_Report.docx",
        "PLamb_Void_Analysis_Project_Thread.docx",
        "JWST CEPHID Super Nova Search.docx",
        "GWTC3_Scalar_Curvature_Bayes_Report_v2.docx",
        "MexicanHat_Scientific_Context_Summary.docx",
        "Vzz_MexicanHat_Model_Description.docx",
    ],
}


def slugify(name: str) -> str:
    path = Path(name)
    stem = path.stem
    stem = re.sub(r"[^A-Za-z0-9._-]+", "_", stem).strip("_")
    stem = re.sub(r"_+", "_", stem)
    return f"{stem}.pdf"


def display_safe(value: object) -> str:
    text = str(value)
    return text.encode("ascii", "backslashreplace").decode("ascii")


def rel_to_root(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT)).replace("\\", "/")
    except ValueError:
        return str(path)


def file_uri(path: Path) -> str:
    return "file:///" + str(path).replace("\\", "/").replace(" ", "%20")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def stream_sources() -> list[tuple[str, Path, str]]:
    items: list[tuple[str, Path, str]] = []

    skymap_manifest = MANIFEST_DIR / f"skymap_report_manifest_{DATE}.csv"
    if skymap_manifest.exists():
        for row in read_csv(skymap_manifest):
            source = Path(row["source_path"])
            if source.suffix.lower() in {".docx", ".rtf"} and source.exists():
                items.append(("skymap_quaia", source, "skymap_report_manifest"))

    hht_manifest = MANIFEST_DIR / f"hht_cmb_sn_docx_report_manifest_{DATE}.csv"
    if hht_manifest.exists():
        for row in read_csv(hht_manifest):
            source = Path(row["source_docx"])
            if source.exists() and not SKYMAP_RE.search(source.name):
                items.append(("hht_cmb_sn_bao", source, "hht_cmb_sn_docx_report_manifest"))

    for stream, names in EXTRA_STREAMS.items():
        for name in names:
            source = ROOT / name
            if source.exists():
                items.append((stream, source, "curated_top_level_report_list"))
            else:
                items.append((stream, source, "missing_curated_top_level_report"))

    seen: set[tuple[str, str]] = set()
    deduped: list[tuple[str, Path, str]] = []
    for stream, source, basis in items:
        key = (stream, str(source).lower())
        if key in seen:
            continue
        seen.add(key)
        deduped.append((stream, source, basis))
    return deduped


def page_count(path: Path) -> int | str:
    try:
        from pypdf import PdfReader  # type: ignore

        return len(PdfReader(str(path)).pages)
    except Exception as exc:
        return f"unreadable: {exc}"


def convert_to_pdf(stream: str, source: Path, basis: str) -> dict[str, object]:
    row: dict[str, object] = {
        "stream": stream,
        "selection_basis": basis,
        "source_path": str(source),
        "source_relative": rel_to_root(source),
        "source_size_bytes": source.stat().st_size if source.exists() else "",
        "source_modified_time": dt.datetime.fromtimestamp(source.stat().st_mtime).isoformat(timespec="seconds")
        if source.exists()
        else "",
        "status": "",
        "pdf_path": "",
        "pdf_relative": "",
        "pdf_size_bytes": "",
        "page_count": "",
        "duration_seconds": "",
        "message": "",
    }
    if not source.exists():
        row["status"] = "missing_source"
        return row
    if not SOFFICE.exists():
        row["status"] = "missing_soffice"
        row["message"] = str(SOFFICE)
        return row

    outdir = PDF_ROOT / stream
    outdir.mkdir(parents=True, exist_ok=True)
    LO_PROFILE.mkdir(parents=True, exist_ok=True)

    existing_out = outdir / f"{source.stem}.pdf"
    final_out = outdir / slugify(source.name)
    if final_out.exists():
        row["status"] = "rendered_existing"
        row["pdf_path"] = str(final_out)
        row["pdf_relative"] = rel_to_root(final_out)
        row["pdf_size_bytes"] = final_out.stat().st_size
        row["page_count"] = page_count(final_out)
        row["duration_seconds"] = "0.0"
        row["message"] = "Reused existing rendered PDF from earlier interrupted run."
        return row
    if existing_out.exists():
        existing_out.unlink()
    if final_out.exists():
        final_out.unlink()

    cmd = [
        str(SOFFICE),
        "--headless",
        f"-env:UserInstallation={file_uri(LO_PROFILE)}",
        "--convert-to",
        "pdf",
        "--outdir",
        str(outdir),
        str(source),
    ]
    started = time.time()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT_SECONDS)
    except subprocess.TimeoutExpired:
        row["status"] = "timeout"
        row["duration_seconds"] = TIMEOUT_SECONDS
        return row

    duration = round(time.time() - started, 2)
    row["duration_seconds"] = duration
    output = (proc.stdout or "") + (proc.stderr or "")
    produced = existing_out if existing_out.exists() else final_out if final_out.exists() else None
    if proc.returncode != 0 or produced is None:
        row["status"] = "failed"
        row["message"] = output[-1000:].replace("\r", " ").replace("\n", " ").strip()
        return row

    if produced != final_out:
        produced.replace(final_out)
    row["status"] = "rendered"
    row["pdf_path"] = str(final_out)
    row["pdf_relative"] = rel_to_root(final_out)
    row["pdf_size_bytes"] = final_out.stat().st_size
    row["page_count"] = page_count(final_out)
    row["message"] = output[-500:].replace("\r", " ").replace("\n", " ").strip()
    return row


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "stream",
        "selection_basis",
        "source_path",
        "source_relative",
        "source_size_bytes",
        "source_modified_time",
        "status",
        "pdf_path",
        "pdf_relative",
        "pdf_size_bytes",
        "page_count",
        "duration_seconds",
        "message",
    ]
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_readmes(rows: list[dict[str, object]]) -> None:
    PDF_ROOT.mkdir(parents=True, exist_ok=True)
    rendered = [r for r in rows if str(r["status"]).startswith("rendered")]
    failed = [r for r in rows if not str(r["status"]).startswith("rendered")]
    by_stream: dict[str, list[dict[str, object]]] = {}
    for row in rendered:
        by_stream.setdefault(str(row["stream"]), []).append(row)

    lines = [
        "# Rendered Research PDFs",
        "",
        f"Generated on {DATE} from selected local DOCX/RTF reports using LibreOffice headless conversion.",
        "",
        "These PDFs are archival reading snapshots. The original Office files remain outside Git; source paths and render details are recorded in `github_export/manifests/rendered_pdf_manifest_2026-07-14.csv`.",
        "",
        "## Streams",
        "",
    ]
    for stream in sorted(by_stream):
        total = sum(int(r["pdf_size_bytes"]) for r in by_stream[stream])
        lines.append(f"- `{stream}/`: {len(by_stream[stream])} PDFs, {total / 1048576:.1f} MiB")
    if failed:
        lines.extend(["", "## Render Notes", "", f"- Non-rendered sources: {len(failed)}. See manifest `status` and `message` columns."])
    lines.append("")
    (PDF_ROOT / "README.md").write_text("\n".join(lines), encoding="utf-8")

    readout = OUTPUT_DIR / f"rendered_research_pdfs_readout_{DATE}.md"
    readout_lines = [
        "# Rendered Research PDF Readout",
        "",
        f"Date: {DATE}",
        "",
        f"- Rendered PDFs: {len(rendered)}",
        f"- Non-rendered sources: {len(failed)}",
        "",
        "## Stream Counts",
        "",
    ]
    for stream in sorted(by_stream):
        total = sum(int(r["pdf_size_bytes"]) for r in by_stream[stream])
        readout_lines.append(f"- `{stream}`: {len(by_stream[stream])} PDFs, {total / 1048576:.1f} MiB")
    if failed:
        readout_lines.extend(["", "## Non-Rendered Sources", ""])
        for row in failed:
            readout_lines.append(f"- `{row['source_relative']}`: {row['status']} {row['message']}")
    readout_lines.extend(
        [
            "",
            "## Git Locations",
            "",
            "- `github_export/docs/rendered_pdfs/`",
            "- `github_export/manifests/rendered_pdf_manifest_2026-07-14.csv`",
            "",
        ]
    )
    readout.write_text("\n".join(readout_lines), encoding="utf-8")


def main() -> int:
    PDF_ROOT.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    sources = stream_sources()
    print(json.dumps({"sources": len(sources), "soffice": str(SOFFICE)}, indent=2), flush=True)

    rows: list[dict[str, object]] = []
    for idx, (stream, source, basis) in enumerate(sources, start=1):
        print(f"[{idx}/{len(sources)}] {stream}: {display_safe(source.name)}", flush=True)
        row = convert_to_pdf(stream, source, basis)
        rows.append(row)
        print(
            f"    -> {display_safe(row['status'])} "
            f"pages={display_safe(row['page_count'])} "
            f"size={display_safe(row['pdf_size_bytes'])}",
            flush=True,
        )

    manifest = MANIFEST_DIR / f"rendered_pdf_manifest_{DATE}.csv"
    write_csv(manifest, rows)
    write_readmes(rows)
    ok_rows = [r for r in rows if str(r["status"]).startswith("rendered")]
    summary = {
        "rendered": len(ok_rows),
        "not_rendered": len(rows) - len(ok_rows),
        "manifest": str(manifest),
        "readout": str(OUTPUT_DIR / f"rendered_research_pdfs_readout_{DATE}.md"),
    }
    print(json.dumps(summary, indent=2), flush=True)
    return 0 if summary["not_rendered"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
