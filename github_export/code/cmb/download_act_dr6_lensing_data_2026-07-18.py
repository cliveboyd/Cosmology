#!/usr/bin/env python3
"""Download, safely extract and hash the official ACT DR6 lensing data."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import tarfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
URL = (
    "https://lambda.gsfc.nasa.gov/data/suborbital/ACT/ACT_dr6/"
    "likelihood/data/ACT_dr6_likelihood_v1.2.tgz"
)
BASE = ROOT / "external_datasets" / "current_cosmology_sources" / "ACT_DR6_lensing"
DEFAULT_ARCHIVE = BASE / "ACT_dr6_likelihood_v1.2.tgz"
DEFAULT_DATA_ROOT = BASE / "data"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    partial = output.with_suffix(output.suffix + ".part")
    offset = partial.stat().st_size if partial.exists() else 0
    request = urllib.request.Request(url, headers={"User-Agent": "Cosmology-CMB-audit/1.0"})
    if offset:
        request.add_header("Range", f"bytes={offset}-")
    with urllib.request.urlopen(request, timeout=120) as response:
        append = offset > 0 and getattr(response, "status", None) == 206
        with partial.open("ab" if append else "wb") as handle:
            shutil.copyfileobj(response, handle, length=1024 * 1024)
    partial.replace(output)


def validate_members(archive: Path) -> list[dict[str, int | str]]:
    rows: list[dict[str, int | str]] = []
    with tarfile.open(archive, "r:gz") as bundle:
        for member in bundle.getmembers():
            if member.issym() or member.islnk():
                raise RuntimeError(f"Links are forbidden in the data archive: {member.name}")
            target = Path(member.name)
            if target.is_absolute() or ".." in target.parts:
                raise RuntimeError(f"Unsafe archive member: {member.name}")
            if member.isfile():
                rows.append({"name": member.name, "bytes": member.size})
    required = {
        "v1.2/clkk_bandpowers_act.txt",
        "v1.2/binning_matrix_act.txt",
        "v1.2/covmat_act.txt",
        "v1.2/like_corrs/cosmo2017_10K_acc3_lensedCls.dat",
    }
    names = {str(row["name"]).replace("\\", "/") for row in rows}
    missing = required - names
    if missing:
        raise RuntimeError(f"ACT data archive is missing registered files: {sorted(missing)}")
    return rows


def extract(archive: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "r:gz") as bundle:
        bundle.extractall(destination, filter="data")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", default=URL)
    parser.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--force-download", action="store_true")
    parser.add_argument("--force-extract", action="store_true")
    args = parser.parse_args()
    if args.force_download or not args.archive.exists():
        download(args.url, args.archive)
    members = validate_members(args.archive)
    data_dir = args.data_root / "v1.2"
    if args.force_extract or not data_dir.exists():
        extract(args.archive, args.data_root)
    if not (data_dir / "covmat_act.txt").exists():
        raise RuntimeError(f"Extraction did not create the expected data directory: {data_dir}")
    manifest = {
        "downloaded_utc": datetime.now(timezone.utc).isoformat(),
        "source_url": args.url,
        "archive": str(args.archive),
        "archive_bytes": args.archive.stat().st_size,
        "archive_sha256": sha256_file(args.archive),
        "data_dir": str(data_dir),
        "file_count": len(members),
        "members": members,
    }
    manifest_path = args.archive.with_suffix(".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Archive SHA-256: {manifest['archive_sha256']}")
    print(f"Data directory: {data_dir}")
    print(f"Files: {len(members)}")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
