#!/usr/bin/env python3
"""Download and hash the official COBE/FIRAS monopole table."""

from __future__ import annotations

import argparse
import hashlib
import json
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
URL = "https://lambda.gsfc.nasa.gov/data/cobe/firas/monopole_spec/firas_monopole_spec_v1.txt"
DEFAULT_OUTPUT = (
    ROOT
    / "external_datasets"
    / "current_cosmology_sources"
    / "COBE_FIRAS"
    / "firas_monopole_spec_v1.txt"
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", default=URL)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(args.url, headers={"User-Agent": "Cosmology-CMB-audit/1.0"})
    with urllib.request.urlopen(request, timeout=60) as response:
        payload = response.read()
    args.output.write_bytes(payload)
    numeric_rows = [
        line
        for line in args.output.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if len(numeric_rows) != 43:
        raise RuntimeError(f"Expected 43 FIRAS channels, found {len(numeric_rows)}")
    manifest = {
        "downloaded_utc": datetime.now(timezone.utc).isoformat(),
        "source_url": args.url,
        "path": str(args.output),
        "bytes": args.output.stat().st_size,
        "sha256": sha256_file(args.output),
        "numeric_rows": len(numeric_rows),
    }
    manifest_path = args.output.with_suffix(".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Saved: {args.output}")
    print(f"SHA-256: {manifest['sha256']}")


if __name__ == "__main__":
    main()
