from __future__ import annotations

import csv
import datetime as dt
import html
import json
import re
import shutil
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET


ROOT = Path(r"C:\Users\clive\Documents\Cosmology")
THREAD_ROOT = Path(r"C:\Users\clive\Documents\Codex\2026-07-14\ok-continu")
EXPORT = ROOT / "github_export"
DATE = "2026-07-14"

CODE_DIR = EXPORT / "code" / "skymap"
DOC_DIR = EXPORT / "docs" / "skymap_reports" / "extracted_text"
MEDIA_DIR = EXPORT / "docs" / "skymap_reports" / "extracted_media"
RESULT_DIR = EXPORT / "results" / DATE / "skymap"
DATA_DIR = EXPORT / "data" / "skymap"
MANIFEST_DIR = EXPORT / "manifests"
OUTPUT_DIR = THREAD_ROOT / "outputs"

DOC_RE = re.compile(
    r"(skymap|sky.?map|rbao|zdipole|quaia|dipole|quadrupole|anisotropy)",
    re.IGNORECASE,
)
CODE_RE = re.compile(r"(skymap|quaia|qso_dipole)", re.IGNORECASE)
TEXT_RESULT_EXTS = {".csv", ".txt", ".json"}
FIGURE_RE = re.compile(
    r"^(SM_RBAO_demo|quaia_G20p0_nside64_zdipole_|quaia_G20p0_random_mocks_dipole_hist)",
    re.IGNORECASE,
)

NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}


def slugify(name: str, *, keep_ext: bool = True) -> str:
    path = Path(name)
    stem = path.stem if keep_ext else name
    ext = path.suffix if keep_ext else ""
    stem = re.sub(r"[^A-Za-z0-9._-]+", "_", stem).strip("_")
    stem = re.sub(r"_+", "_", stem)
    return f"{stem}{ext.lower()}"


def ensure_dirs() -> None:
    for path in [CODE_DIR, DOC_DIR, MEDIA_DIR, RESULT_DIR, DATA_DIR, MANIFEST_DIR, OUTPUT_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def rel_to_root(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT)).replace("\\", "/")
    except ValueError:
        return str(path)


def file_meta(path: Path) -> dict[str, str | int]:
    stat = path.stat()
    return {
        "source_path": str(path),
        "source_relative": rel_to_root(path),
        "size_bytes": stat.st_size,
        "modified_time": dt.datetime.fromtimestamp(stat.st_mtime).isoformat(timespec="seconds"),
    }


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def para_text(elem: ET.Element) -> str:
    parts: list[str] = []
    for child in elem.iter():
        tag = child.tag.rsplit("}", 1)[-1]
        if tag == "t" and child.text:
            parts.append(child.text)
        elif tag == "tab":
            parts.append("\t")
        elif tag in {"br", "cr"}:
            parts.append("\n")
    return "".join(parts).strip()


def extract_docx_text(path: Path) -> tuple[str, dict[str, str]]:
    with zipfile.ZipFile(path) as zf:
        doc_xml = zf.read("word/document.xml")
        root = ET.fromstring(doc_xml)
        body = root.find("w:body", NS)
        lines: list[str] = []
        if body is not None:
            for elem in list(body):
                tag = elem.tag.rsplit("}", 1)[-1]
                if tag == "p":
                    text = para_text(elem)
                    if text:
                        lines.append(text)
                elif tag == "tbl":
                    for row in elem.findall(".//w:tr", NS):
                        cells = []
                        for cell in row.findall("./w:tc", NS):
                            cell_text = " ".join(
                                para_text(p) for p in cell.findall(".//w:p", NS) if para_text(p)
                            ).strip()
                            cells.append(cell_text)
                        if any(cells):
                            lines.append(" | ".join(cells))

        props: dict[str, str] = {}
        for member, prefix in [
            ("docProps/core.xml", "core"),
            ("docProps/app.xml", "app"),
        ]:
            if member not in zf.namelist():
                continue
            try:
                props_xml = ET.fromstring(zf.read(member))
            except ET.ParseError:
                continue
            for item in props_xml.iter():
                tag = item.tag.rsplit("}", 1)[-1]
                text = (item.text or "").strip()
                if text:
                    props[f"{prefix}_{tag}"] = text
    return "\n\n".join(lines).strip(), props


def extract_rtf_text(path: Path) -> str:
    raw = path.read_text(encoding="latin-1", errors="ignore")
    raw = re.sub(r"\\'[0-9a-fA-F]{2}", " ", raw)
    raw = re.sub(r"\\[a-zA-Z]+\d* ?", " ", raw)
    raw = raw.replace("{", " ").replace("}", " ").replace("\\", " ")
    raw = html.unescape(raw)
    lines = [re.sub(r"\s+", " ", line).strip() for line in raw.splitlines()]
    return "\n\n".join(line for line in lines if line)


def extract_pdf_text(path: Path) -> str:
    try:
        from pypdf import PdfReader  # type: ignore
    except Exception:
        return "[PDF text extraction unavailable: pypdf is not installed in this runtime.]"

    reader = PdfReader(str(path))
    pages = []
    for idx, page in enumerate(reader.pages, start=1):
        text = page.extract_text() or ""
        pages.append(f"## Page {idx}\n\n{text.strip()}")
    return "\n\n".join(pages).strip()


def doc_candidates() -> list[Path]:
    candidates = []
    for path in ROOT.iterdir():
        if path.is_file() and path.suffix.lower() in {".docx", ".pdf", ".rtf"} and DOC_RE.search(path.name):
            candidates.append(path)
    return sorted(candidates, key=lambda p: p.name.lower())


def code_candidates() -> list[Path]:
    candidates: list[Path] = []
    for path in [
        THREAD_ROOT / "work" / "capture_skymap_archive.py",
        ROOT / "SM_SkyMap_FR_KdS.py",
        ROOT / "inventory_quaia_existing_work.py",
        ROOT / "diagnose_su2_quaia_scan.py",
        ROOT / "diagnose_su2_quaia_mock_controls.py",
        ROOT / "diagnose_su2_quaia_hht_controls.py",
    ]:
        if path.exists():
            candidates.append(path)
    spyder = ROOT / ".spyder-py3"
    for path in spyder.glob("*.py"):
        if CODE_RE.search(path.name) or path.name in {
            "hht_quaia_residuals.py",
            "stack_quaia_imfs.py",
            "run_quaia_imf_robustness.py",
            "postprocess_quaia_shift_null.py",
            "overnight_quaia_hunt.py",
            "overnight_quaia_pipeline.py",
        }:
            candidates.append(path)
    unique: dict[str, Path] = {}
    for path in candidates:
        unique[str(path).lower()] = path
    return sorted(unique.values(), key=lambda p: rel_to_root(p).lower())


def module_summary(path: Path) -> str:
    text = path.read_text(encoding="utf-8", errors="replace")
    match = re.search(r'^\s*(?:(?:r|u|f|fr|rf|ur|ru)?("""|\'\'\'))(.*?)(?:\1)', text, re.S | re.I)
    if match:
        summary = re.sub(r"\s+", " ", match.group(2)).strip()
    else:
        lines = []
        for line in text.splitlines()[:30]:
            stripped = line.strip()
            if stripped.startswith("#") and not stripped.startswith("#!"):
                lines.append(stripped.lstrip("#").strip())
        summary = " ".join(lines)
    return summary[:500]


def copy_code() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    seen_names: set[str] = set()
    for path in code_candidates():
        base = slugify(path.name)
        if path.parent.name == ".spyder-py3":
            base = f"spyder_{base}"
        dest_name = base
        n = 2
        while dest_name.lower() in seen_names:
            dest_name = f"{Path(base).stem}_{n}{Path(base).suffix}"
            n += 1
        seen_names.add(dest_name.lower())
        dest = CODE_DIR / dest_name
        shutil.copy2(path, dest)
        meta = file_meta(path)
        rows.append(
            {
                **meta,
                "export_path": rel_to_root(dest),
                "module_summary": module_summary(path),
            }
        )
    return rows


def extract_docx_media(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    if path.suffix.lower() != ".docx":
        return rows
    doc_media_dir = MEDIA_DIR / slugify(path.stem, keep_ext=False)
    doc_media_dir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path) as zf:
        media_members = [
            item for item in zf.infolist() if item.filename.startswith("word/media/") and not item.is_dir()
        ]
        for idx, item in enumerate(media_members, start=1):
            ext = Path(item.filename).suffix.lower() or ".bin"
            dest = doc_media_dir / f"image_{idx:03d}{ext}"
            with zf.open(item) as src, dest.open("wb") as out:
                shutil.copyfileobj(src, out)
            rows.append(
                {
                    "source_path": str(path),
                    "source_relative": rel_to_root(path),
                    "source_media_member": item.filename,
                    "media_size_bytes": item.file_size,
                    "export_path": rel_to_root(dest),
                }
            )
    return rows


def extract_docs() -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    rows: list[dict[str, object]] = []
    media_rows: list[dict[str, object]] = []
    for path in doc_candidates():
        suffix = path.suffix.lower()
        props: dict[str, str] = {}
        this_media_rows: list[dict[str, object]] = []
        try:
            if suffix == ".docx":
                text, props = extract_docx_text(path)
                this_media_rows = extract_docx_media(path)
            elif suffix == ".rtf":
                text = extract_rtf_text(path)
            elif suffix == ".pdf":
                text = extract_pdf_text(path)
            else:
                continue
        except Exception as exc:
            text = f"[Extraction failed: {exc}]"

        dest = DOC_DIR / f"{slugify(path.name, keep_ext=True)}.md"
        meta = file_meta(path)
        word_count = len(re.findall(r"\b\w+\b", text))
        heading = path.stem.replace("_", " ")
        front = [
            f"# {heading}",
            "",
            f"- Source: `{meta['source_relative']}`",
            f"- Source size: {meta['size_bytes']} bytes",
            f"- Source modified: {meta['modified_time']}",
            f"- Extracted: {DATE}",
            f"- Word count estimate: {word_count}",
        ]
        title = props.get("core_title")
        if title:
            front.append(f"- Document title property: {title}")
        front.append("")
        front.append("## Extracted Text")
        front.append("")
        dest.write_text("\n".join(front) + (text or "[No text extracted.]") + "\n", encoding="utf-8")
        rows.append(
            {
                **meta,
                "export_path": rel_to_root(dest),
                "word_count_estimate": word_count,
                "doc_title": title or "",
                "embedded_media_files": len(this_media_rows),
                "embedded_media_bytes": sum(int(row["media_size_bytes"]) for row in this_media_rows),
            }
        )
        media_rows.extend(this_media_rows)
    return rows, media_rows


def result_candidates() -> list[Path]:
    paths: list[Path] = []
    skymap_run = ROOT / ".spyder-py3" / "plamb_runs" / "SkyMap"
    if skymap_run.exists():
        for path in skymap_run.rglob("*"):
            if path.is_file() and path.suffix.lower() in TEXT_RESULT_EXTS and "\\logs\\" not in str(path).lower():
                paths.append(path)
    quaia_outputs = ROOT / ".spyder-py3" / "skymap2" / "quaia_outputs"
    if quaia_outputs.exists():
        for path in quaia_outputs.rglob("*"):
            if path.is_file() and path.suffix.lower() in TEXT_RESULT_EXTS:
                paths.append(path)
    for name in [
        "quaia_G20p0_random_mocks_b30_stats.csv",
        "quaia_G20p0_random_mocks_stats.csv",
        "quaia_G20p0_random_mocks_zdipole_quad.csv",
    ]:
        path = ROOT / name
        if path.exists():
            paths.append(path)
    return sorted({str(p).lower(): p for p in paths}.values(), key=lambda p: rel_to_root(p).lower())


def copy_results() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for path in result_candidates():
        rel = rel_to_root(path)
        if ".spyder-py3/plamb_runs/SkyMap".lower() in rel.lower():
            sub = Path(rel.split(".spyder-py3/plamb_runs/SkyMap/", 1)[1])
            dest = RESULT_DIR / "plamb_runs_SkyMap" / sub
        elif ".spyder-py3/skymap2/quaia_outputs".lower() in rel.lower():
            sub = Path(rel.split(".spyder-py3/skymap2/quaia_outputs/", 1)[1])
            dest = RESULT_DIR / "skymap2_quaia_outputs" / sub
        else:
            dest = RESULT_DIR / "top_level_quaia_mock_summaries" / path.name
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, dest)
        meta = file_meta(path)
        rows.append({**meta, "export_path": rel_to_root(dest)})
    return rows


def copy_data_and_figures() -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    data_rows: list[dict[str, object]] = []
    figure_rows: list[dict[str, object]] = []

    for name in ["skymap_unified_catalog_sample.csv"]:
        path = ROOT / name
        if path.exists():
            dest = DATA_DIR / path.name
            shutil.copy2(path, dest)
            data_rows.append({**file_meta(path), "export_path": rel_to_root(dest), "note": "sample catalog"})

    for path in ROOT.iterdir():
        if path.is_file() and path.suffix.lower() == ".png" and FIGURE_RE.search(path.name):
            dest = RESULT_DIR / "figures" / path.name
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, dest)
            figure_rows.append({**file_meta(path), "export_path": rel_to_root(dest)})

    omitted = [
        ROOT / "skymap_unified_catalog_sample (1).csv",
        ROOT / "quaia_G20p0_basic.csv",
        ROOT / "quaia_G20p0_basic_gal.csv",
    ]
    for path in omitted:
        if path.exists():
            data_rows.append(
                {
                    **file_meta(path),
                    "export_path": "",
                    "note": "referenced only; omitted to avoid duplicate or large raw-data Git payload",
                }
            )
    return data_rows, figure_rows


def write_indexes(
    doc_rows: list[dict[str, object]],
    media_rows: list[dict[str, object]],
    code_rows: list[dict[str, object]],
    result_rows: list[dict[str, object]],
    data_rows: list[dict[str, object]],
    figure_rows: list[dict[str, object]],
) -> None:
    readme = CODE_DIR / "README.md"
    readme.write_text(
        "\n".join(
            [
                "# SkyMap Code Capture",
                "",
                f"Captured on {DATE}. This directory contains curated Python scripts for the SkyMap, R_BAO anisotropy, and Quaia redshift-dipole work.",
                "",
                "The copied scripts are source snapshots. The source paths, modified times, and short module summaries are recorded in `github_export/manifests/skymap_code_manifest_2026-07-14.csv`.",
                "",
                "Large raw catalogs and FITS/NPZ map products are not committed here; the manifests record the local source locations where relevant.",
                "",
            ]
        ),
        encoding="utf-8",
    )

    report_index = OUTPUT_DIR / f"skymap_capture_readout_{DATE}.md"
    doc_links = "\n".join(
        f"- `{row['source_relative']}` -> `{row['export_path']}` ({row['word_count_estimate']} words)"
        for row in doc_rows
    )
    code_links = "\n".join(
        f"- `{row['source_relative']}` -> `{row['export_path']}`"
        for row in code_rows[:20]
    )
    if len(code_rows) > 20:
        code_links += f"\n- ... plus {len(code_rows) - 20} more scripts in the code manifest."

    report_index.write_text(
        "\n".join(
            [
                "# SkyMap Git Capture Readout",
                "",
                f"Date: {DATE}",
                "",
                "## What Was Captured",
                "",
                f"- Extracted narrative reports: {len(doc_rows)} DOCX/RTF/PDF files to Markdown.",
                f"- Extracted embedded DOCX media files: {len(media_rows)}.",
                f"- Copied Python scripts: {len(code_rows)} curated SkyMap/Quaia/dipole scripts.",
                f"- Copied small text/CSV/JSON result files: {len(result_rows)}.",
                f"- Copied sample data files: {sum(1 for row in data_rows if row.get('export_path'))}.",
                f"- Copied representative figures: {len(figure_rows)}.",
                "",
                "No matching PDF reports were found in the top-level SkyMap/Quaia report scan at capture time.",
                "",
                "## Claim Boundary",
                "",
                "This archive is a provenance capture, not a new statistical endorsement. The SkyMap/Quaia material should be read as historical and methodological support for anisotropy follow-up. The later SU2/Quaia global look-elsewhere and template-regression gates remain the default claim boundary for current interpretation.",
                "",
                "## Extracted Reports",
                "",
                doc_links or "- None.",
                "",
                "## Code Snapshot",
                "",
                code_links or "- None.",
                "",
                "## Git Locations",
                "",
                "- `github_export/code/skymap/`",
                "- `github_export/docs/skymap_reports/extracted_text/`",
                "- `github_export/results/2026-07-14/skymap/`",
                "- `github_export/data/skymap/`",
                "- `github_export/manifests/skymap_*_manifest_2026-07-14.csv`",
                "",
            ]
        ),
        encoding="utf-8",
    )

    docs_index = DOC_DIR.parent / "README.md"
    docs_index.write_text(
        "\n".join(
            [
                "# SkyMap Report Text Capture",
                "",
                f"Captured on {DATE}. This folder stores Markdown text extracted from SkyMap, R_BAO, Quaia dipole, and related DOCX/RTF reports.",
                "",
                "Original DOCX/RTF binaries are not committed in this Git export. Extracted text and embedded DOCX media are captured here; use the manifests to locate the local originals if full Word layout is needed.",
                "",
                "See `github_export/manifests/skymap_report_manifest_2026-07-14.csv` and `github_export/manifests/skymap_report_media_manifest_2026-07-14.csv` for source paths and extraction metadata.",
                "",
            ]
        ),
        encoding="utf-8",
    )

    result_readme = RESULT_DIR / "README.md"
    result_readme.write_text(
        "\n".join(
            [
                "# SkyMap Result Capture",
                "",
                f"Captured on {DATE}. This directory stores small text, CSV, and JSON result summaries from the legacy SkyMap and Quaia dipole-map work.",
                "",
                "Heavy binary products such as FITS, NPZ, and most rendered PNG map grids are intentionally omitted. Representative small figures from the top-level working folder are included under `figures/`.",
                "",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    ensure_dirs()
    code_rows = copy_code()
    doc_rows, media_rows = extract_docs()
    result_rows = copy_results()
    data_rows, figure_rows = copy_data_and_figures()

    write_csv(
        MANIFEST_DIR / f"skymap_code_manifest_{DATE}.csv",
        code_rows,
        ["source_path", "source_relative", "size_bytes", "modified_time", "export_path", "module_summary"],
    )
    write_csv(
        MANIFEST_DIR / f"skymap_report_manifest_{DATE}.csv",
        doc_rows,
        [
            "source_path",
            "source_relative",
            "size_bytes",
            "modified_time",
            "export_path",
            "word_count_estimate",
            "doc_title",
            "embedded_media_files",
            "embedded_media_bytes",
        ],
    )
    write_csv(
        MANIFEST_DIR / f"skymap_report_media_manifest_{DATE}.csv",
        media_rows,
        ["source_path", "source_relative", "source_media_member", "media_size_bytes", "export_path"],
    )
    write_csv(
        MANIFEST_DIR / f"skymap_result_manifest_{DATE}.csv",
        result_rows,
        ["source_path", "source_relative", "size_bytes", "modified_time", "export_path"],
    )
    write_csv(
        MANIFEST_DIR / f"skymap_data_manifest_{DATE}.csv",
        data_rows,
        ["source_path", "source_relative", "size_bytes", "modified_time", "export_path", "note"],
    )
    write_csv(
        MANIFEST_DIR / f"skymap_figure_manifest_{DATE}.csv",
        figure_rows,
        ["source_path", "source_relative", "size_bytes", "modified_time", "export_path"],
    )

    write_indexes(doc_rows, media_rows, code_rows, result_rows, data_rows, figure_rows)
    print(
        json.dumps(
            {
                "code_files": len(code_rows),
                "reports": len(doc_rows),
                "report_media": len(media_rows),
                "results": len(result_rows),
                "data_rows": len(data_rows),
                "figures": len(figure_rows),
                "readout": str(OUTPUT_DIR / f"skymap_capture_readout_{DATE}.md"),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
