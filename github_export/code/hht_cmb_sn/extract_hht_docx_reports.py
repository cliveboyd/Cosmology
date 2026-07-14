from __future__ import annotations

import argparse
import csv
import re
import zipfile
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from xml.etree import ElementTree as ET


NS = {
    "w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main",
    "cp": "http://schemas.openxmlformats.org/package/2006/metadata/core-properties",
    "dc": "http://purl.org/dc/elements/1.1/",
    "dcterms": "http://purl.org/dc/terms/",
}


KEYWORDS = [
    "hht",
    "hilbert",
    "huang",
    "locking",
    "dechirp",
    "imf",
    "cmb",
    "pantheon",
    "shoes",
    "snia",
    "bao",
    "quaia",
    "skymap",
]


@dataclass
class Report:
    source: Path
    slug: str
    title: str
    created: str
    modified: str
    paragraphs: list[str]
    tables: list[list[list[str]]]


def slugify(path: Path) -> str:
    text = path.stem.lower()
    text = re.sub(r"[^a-z0-9]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text[:130] or "report"


def xml_text(node: ET.Element) -> str:
    chunks: list[str] = []
    for text in node.iter(f"{{{NS['w']}}}t"):
        if text.text:
            chunks.append(text.text)
    return "".join(chunks).strip()


def read_core_props(zf: zipfile.ZipFile) -> dict[str, str]:
    out: dict[str, str] = {}
    try:
        root = ET.fromstring(zf.read("docProps/core.xml"))
    except Exception:
        return out
    for key, xpath in {
        "title": "dc:title",
        "subject": "dc:subject",
        "creator": "dc:creator",
        "created": "dcterms:created",
        "modified": "dcterms:modified",
    }.items():
        node = root.find(xpath, NS)
        if node is not None and node.text:
            out[key] = node.text.strip()
    return out


def extract_docx(path: Path, slug_counts: dict[str, int]) -> Report | None:
    try:
        with zipfile.ZipFile(path) as zf:
            props = read_core_props(zf)
            root = ET.fromstring(zf.read("word/document.xml"))
    except Exception:
        return None

    body = root.find("w:body", NS)
    paragraphs: list[str] = []
    tables: list[list[list[str]]] = []
    if body is not None:
        for child in list(body):
            tag = child.tag.rsplit("}", 1)[-1]
            if tag == "p":
                text = xml_text(child)
                if text:
                    paragraphs.append(text)
            elif tag == "tbl":
                rows: list[list[str]] = []
                for tr in child.findall("w:tr", NS):
                    cells = []
                    for tc in tr.findall("w:tc", NS):
                        cell_text = " ".join(
                            xml_text(p) for p in tc.findall(".//w:p", NS) if xml_text(p)
                        ).strip()
                        cells.append(cell_text)
                    if cells:
                        rows.append(cells)
                if rows:
                    tables.append(rows)

    slug = slugify(path)
    count = slug_counts.get(slug, 0)
    slug_counts[slug] = count + 1
    if count:
        slug = f"{slug}_{count + 1}"

    title = props.get("title") or path.stem
    return Report(
        source=path,
        slug=slug,
        title=title,
        created=props.get("created", ""),
        modified=props.get("modified", ""),
        paragraphs=paragraphs,
        tables=tables,
    )


def matches(path: Path) -> bool:
    if path.name.startswith("._"):
        return False
    haystack = str(path).lower()
    if not any(keyword in haystack for keyword in KEYWORDS):
        return False
    return any(
        token in haystack
        for token in [
            "hht",
            "locking",
            "dechirp",
            "imf",
            "snia",
            "bao_hht",
            "hht_bao",
            "fr_hht",
            "skymap_hht",
            "quaia hht",
        ]
    )


def markdown_for(report: Report) -> str:
    lines: list[str] = [
        f"# {report.title}",
        "",
        "## Source",
        "",
        f"- Local DOCX: `{report.source}`",
        f"- Extracted: `{datetime.now().isoformat(timespec='seconds')}`",
    ]
    if report.created:
        lines.append(f"- DOCX created: `{report.created}`")
    if report.modified:
        lines.append(f"- DOCX modified: `{report.modified}`")
    lines.extend(
        [
            "",
            "## Extracted Text",
            "",
        ]
    )
    for paragraph in report.paragraphs:
        lines.append(paragraph)
        lines.append("")
    if report.tables:
        lines.extend(["## Extracted Tables", ""])
        for i, table in enumerate(report.tables, start=1):
            lines.append(f"### Table {i}")
            lines.append("")
            for row in table:
                lines.append(" | ".join(cell.replace("\n", " ") for cell in row))
            lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()

    docs = sorted(
        [p for p in args.root.rglob("*.docx") if matches(p)],
        key=lambda p: (p.stat().st_mtime, str(p)),
        reverse=True,
    )
    if args.limit:
        docs = docs[: args.limit]

    args.out.mkdir(parents=True, exist_ok=True)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    slug_counts: dict[str, int] = {}
    rows: list[dict[str, object]] = []
    for path in docs:
        report = extract_docx(path, slug_counts)
        if report is None:
            rows.append(
                {
                    "source_docx": str(path),
                    "status": "extract_failed",
                    "markdown": "",
                    "title": "",
                    "bytes": path.stat().st_size,
                    "mtime": datetime.fromtimestamp(path.stat().st_mtime).isoformat(timespec="seconds"),
                    "paragraphs": 0,
                    "tables": 0,
                }
            )
            continue
        md_name = f"{report.slug}.md"
        md_path = args.out / md_name
        md_path.write_text(markdown_for(report), encoding="utf-8", newline="\n")
        rows.append(
            {
                "source_docx": str(path),
                "status": "extracted",
                "markdown": str(md_path),
                "title": report.title,
                "bytes": path.stat().st_size,
                "mtime": datetime.fromtimestamp(path.stat().st_mtime).isoformat(timespec="seconds"),
                "paragraphs": len(report.paragraphs),
                "tables": len(report.tables),
            }
        )

    with args.manifest.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "source_docx",
                "status",
                "markdown",
                "title",
                "bytes",
                "mtime",
                "paragraphs",
                "tables",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Scanned DOCX reports: {len(docs)}")
    print(f"Manifest: {args.manifest}")
    print(f"Markdown output: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
