from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
import urllib.request
from pathlib import Path
from typing import Any

from astropy.io import fits


ROOT = Path(r"C:\Users\clive\Documents\Cosmology")
DEFAULT_PREREG = (
    ROOT
    / "github_export"
    / "results"
    / "2026-07-18"
    / "su2_quaia"
    / "sdss_dr16q_v4_cross_catalogue_preregistration_2026-07-18.json"
)
DEFAULT_DESTINATION = ROOT / "external_datasets" / "sdss_dr16q_v4" / "DR16Q_v4.fits"
CHUNK_BYTES = 8 * 1024 * 1024


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK_BYTES), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_preregistration(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        prereg = json.load(handle)
    if prereg.get("outcome_data_accessed_before_lock") is not False:
        raise ValueError("Preregistration does not certify a pre-outcome lock")
    return prereg


def acquire_lock(path: Path) -> int:
    try:
        descriptor = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError as exc:
        owner = path.read_text(encoding="ascii", errors="replace").strip()
        raise RuntimeError(f"Download lock already exists: {path} (owner {owner or 'unknown'})") from exc
    os.write(descriptor, f"pid={os.getpid()}\nstarted={time.strftime('%Y-%m-%dT%H:%M:%S%z')}\n".encode("ascii"))
    return descriptor


def download_resumable(url: str, destination: Path, expected_bytes: int) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    if destination.exists():
        size = destination.stat().st_size
        if size == expected_bytes:
            print(f"Official file already present with expected size: {destination}", flush=True)
            return destination
        raise RuntimeError(f"Existing destination has unexpected size {size}: {destination}")

    offset = partial.stat().st_size if partial.exists() else 0
    if offset > expected_bytes:
        raise RuntimeError(f"Partial file is larger than expected ({offset} > {expected_bytes}): {partial}")
    if offset == expected_bytes:
        os.replace(partial, destination)
        return destination

    headers = {"User-Agent": "Cosmology-DR16Q-validation/1.0"}
    if offset:
        headers["Range"] = f"bytes={offset}-"
        print(f"Resuming at byte {offset:,} of {expected_bytes:,}", flush=True)
    else:
        print(f"Downloading {expected_bytes:,} bytes from {url}", flush=True)

    request = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(request, timeout=120) as response:
        status = int(getattr(response, "status", response.getcode()))
        if offset and status != 206:
            print("Server did not honour Range; restarting the partial download", flush=True)
            offset = 0
            mode = "wb"
        else:
            mode = "ab" if offset else "wb"

        next_report = offset + 64 * 1024 * 1024
        started = time.monotonic()
        written = offset
        with partial.open(mode) as handle:
            while True:
                chunk = response.read(CHUNK_BYTES)
                if not chunk:
                    break
                handle.write(chunk)
                written += len(chunk)
                if written >= next_report:
                    elapsed = max(time.monotonic() - started, 1e-6)
                    rate = (written - offset) / elapsed / (1024 * 1024)
                    print(
                        f"downloaded={written:,}/{expected_bytes:,} ({100.0 * written / expected_bytes:.1f}%) "
                        f"rate={rate:.1f} MiB/s",
                        flush=True,
                    )
                    next_report = written + 64 * 1024 * 1024

    actual = partial.stat().st_size
    if actual != expected_bytes:
        raise RuntimeError(
            f"Download stopped at {actual:,} of {expected_bytes:,} bytes; rerun to resume from the partial file"
        )
    os.replace(partial, destination)
    print(f"Download complete: {destination}", flush=True)
    return destination


def validate_fits(path: Path, catalogue: dict[str, Any]) -> dict[str, Any]:
    expected_bytes = int(catalogue["expected_bytes"])
    actual_bytes = path.stat().st_size
    if actual_bytes != expected_bytes:
        raise ValueError(f"File size mismatch: expected {expected_bytes}, found {actual_bytes}")

    with fits.open(path, memmap=True, mode="readonly") as hdul:
        if len(hdul) < 2 or not isinstance(hdul[1], fits.BinTableHDU):
            raise ValueError("Expected a binary table in FITS extension 1")
        table = hdul[1]
        actual_rows = int(table.header.get("NAXIS2", -1))
        expected_rows = int(catalogue["expected_rows"])
        if actual_rows != expected_rows:
            raise ValueError(f"Row-count mismatch: expected {expected_rows}, found {actual_rows}")
        names = {name.upper() for name in (table.columns.names or [])}
        missing = [name for name in catalogue["required_columns"] if name.upper() not in names]
        if missing:
            raise ValueError(f"Required FITS columns are missing: {missing}")

    print("FITS structure validated; computing SHA-256", flush=True)
    return {
        "path": str(path),
        "bytes": actual_bytes,
        "rows": actual_rows,
        "sha256": sha256(path),
        "validated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "required_columns_present": True,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Resumably download and structurally validate official SDSS DR16Q v4.")
    parser.add_argument("--preregistration", type=Path, default=DEFAULT_PREREG)
    parser.add_argument("--destination", type=Path, default=DEFAULT_DESTINATION)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    prereg = load_preregistration(args.preregistration)
    catalogue = prereg["catalogue"]
    lock_path = args.destination.with_suffix(args.destination.suffix + ".download.lock")
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = acquire_lock(lock_path)
    try:
        path = download_resumable(str(catalogue["url"]), args.destination, int(catalogue["expected_bytes"]))
        manifest = validate_fits(path, catalogue)
        manifest.update(
            {
                "source_url": catalogue["url"],
                "preregistration": str(args.preregistration),
                "preregistration_sha256": sha256(args.preregistration),
            }
        )
        manifest_path = path.with_suffix(".download_manifest.json")
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        print(f"Saved download manifest: {manifest_path}", flush=True)
    finally:
        os.close(descriptor)
        lock_path.unlink(missing_ok=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr, flush=True)
        raise
