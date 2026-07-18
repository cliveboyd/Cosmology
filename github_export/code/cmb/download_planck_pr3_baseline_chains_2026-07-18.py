#!/usr/bin/env python3
"""Download, validate and hash the registered Planck PR3 chain bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
URL = (
    "https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/"
    "cosmoparams/COM_CosmoParams_base-plikHM-TTTEEE-lowl-lowE_R3.00.zip"
)
DEFAULT_OUTPUT = (
    ROOT
    / "external_datasets"
    / "current_cosmology_sources"
    / "Planck_PR3_chains"
    / "COM_CosmoParams_base-plikHM-TTTEEE-lowl-lowE_R3.00.zip"
)


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
        mode = "ab" if append else "wb"
        if offset and not append:
            offset = 0
        with partial.open(mode) as handle:
            shutil.copyfileobj(response, handle, length=1024 * 1024)
    partial.replace(output)


def validate_zip(path: Path) -> list[dict[str, int | str]]:
    with zipfile.ZipFile(path) as archive:
        bad_member = archive.testzip()
        if bad_member is not None:
            raise RuntimeError(f"Corrupt ZIP member: {bad_member}")
        members = [
            {
                "name": info.filename,
                "compressed_bytes": info.compress_size,
                "uncompressed_bytes": info.file_size,
                "crc32": f"{info.CRC:08x}",
            }
            for info in archive.infolist()
            if not info.is_dir()
        ]
    if not any(str(item["name"]).endswith(".paramnames") for item in members):
        raise RuntimeError("No GetDist parameter-name file found in the chain bundle")
    if not any(str(item["name"]).endswith("_1.txt") for item in members):
        raise RuntimeError("No first GetDist chain file found in the bundle")
    return members


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", default=URL)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    if args.force or not args.output.exists():
        download(args.url, args.output)
    members = validate_zip(args.output)
    manifest = {
        "downloaded_utc": datetime.now(timezone.utc).isoformat(),
        "source_url": args.url,
        "path": str(args.output),
        "bytes": args.output.stat().st_size,
        "sha256": sha256_file(args.output),
        "member_count": len(members),
        "members": members,
    }
    manifest_path = args.output.with_suffix(".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Saved: {args.output}")
    print(f"SHA-256: {manifest['sha256']}")
    print(f"Members: {len(members)}")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
