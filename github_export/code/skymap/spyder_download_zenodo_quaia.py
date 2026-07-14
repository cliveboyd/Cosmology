#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 16 20:57:00 2025

@author: boyde
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
download_zenodo_quaia.py — Download all files from a Zenodo record (0V2)

Default DOI:
    10.5281/zenodo.8060755   (Quaia: The Gaia–unWISE Quasar Catalog)

Usage examples
--------------
Download all files from the default Quaia DOI into the current directory:

    python download_zenodo_quaia.py

Specify a different DOI or record ID and output directory:

    python download_zenodo_quaia.py \
        --doi    10.5281/zenodo.8098636 \
        --outdir "/Users/boyde/Downloads/quaia_data"

Requires:
    requests  (pip install requests)
"""

import argparse
import os
import sys
import textwrap
from typing import Dict, Any

import requests


API_BASE     = "https://zenodo.org/api/records"
CHUNK_SIZE   = 1024 * 1024    # 1 MB
MAX_RETRIES  = 3


def parse_record_id(doi_or_id: str) -> str:
    """
    Extract the Zenodo record ID from a DOI, URL, or bare numeric ID.

    Accepts inputs like:
        "10.5281/zenodo.8060755"
        "https://zenodo.org/records/8060755"
        "8060755"
    """
    s             = doi_or_id.strip()

    # Case 1: bare numeric ID
    if s.isdigit():
        return s

    # Case 2: full DOI including "zenodo.<id>"
    if "zenodo." in s:
        return s.split("zenodo.")[-1]

    # Case 3: record URL
    if "zenodo.org/records/" in s:
        return s.rsplit("/", 1)[-1]

    raise ValueError(f"Could not parse Zenodo record ID from: {doi_or_id!r}")


def fetch_record_metadata(record_id: str) -> Dict[str, Any]:
    """
    Query the Zenodo API for a record and return the JSON metadata.
    """
    url           = f"{API_BASE}/{record_id}"
    resp          = requests.get(url, timeout = 30)

    if resp.status_code != 200:
        raise RuntimeError(
            f"Zenodo API request failed ({resp.status_code}) for {url}"
        )

    return resp.json()


def download_file(file_info: Dict[str, Any],
                  outdir: str,
                  overwrite: bool = False) -> None:
    """
    Download a single file from its Zenodo file metadata entry.

    Retries up to MAX_RETRIES times if a network error occurs.
    If a partial file exists and overwrite is False, we check its size
    against the expected size and force re-download if incomplete.
    """
    filename      = file_info.get("key") or file_info.get("filename")
    if not filename:
        raise RuntimeError("Could not determine filename from file_info")

    url           = file_info["links"].get("self") or file_info["links"].get("download")
    if not url:
        raise RuntimeError(f"No download link found for file {filename!r}")

    os.makedirs(outdir, exist_ok = True)

    out_path      = os.path.join(outdir, filename)
    size_bytes    = file_info.get("size", None)

    # If file exists and overwrite is False, check if it looks complete.
    if os.path.exists(out_path) and not overwrite:
        if isinstance(size_bytes, (int, float)):
            local_size  = os.path.getsize(out_path)
            if local_size == size_bytes:
                print(f"[skip] {filename} (already exists with correct size)")
                return
            else:
                print(
                    f"[warn] {filename} exists but size mismatch "
                    f"(local {local_size} vs expected {size_bytes}); "
                    f"will re-download."
                )
                overwrite  = True
        else:
            print(f"[skip] {filename} (already exists; no size info to check)")
            return

    size_str      = (
        f"{size_bytes / (1024 ** 2):.2f} MB"
        if isinstance(size_bytes, (int, float)) else "unknown size"
    )

    for attempt in range(1, MAX_RETRIES + 1):
        if overwrite and os.path.exists(out_path):
            try:
                os.remove(out_path)
            except OSError:
                pass

        print(
            f"[download] {filename}  ({size_str})  "
            f"[attempt {attempt}/{MAX_RETRIES}]"
        )
        print(f"           -> {out_path}")

        try:
            with requests.get(url, stream = True, timeout = 60) as r:
                r.raise_for_status()
                downloaded  = 0
                with open(out_path, "wb") as f:
                    for chunk in r.iter_content(chunk_size = CHUNK_SIZE):
                        if not chunk:
                            continue
                        f.write(chunk)
                        downloaded += len(chunk)

            # If we know the expected size, check it.
            if isinstance(size_bytes, (int, float)) and downloaded != size_bytes:
                print(
                    f"[warn] {filename}: downloaded size {downloaded} bytes "
                    f"!= expected {size_bytes} bytes; retrying."
                )
                continue

            print(
                f"[done] {filename} "
                f"(downloaded {downloaded / (1024 ** 2):.2f} MB)"
            )
            return

        except Exception as e:
            print(f"[error] {filename}: attempt {attempt}/{MAX_RETRIES} failed: {e}")
            try:
                if os.path.exists(out_path):
                    os.remove(out_path)
            except OSError:
                pass

            if attempt == MAX_RETRIES:
                print(f"[fail] {filename}: giving up after {MAX_RETRIES} attempts.")
                raise


def build_parser() -> argparse.ArgumentParser:
    parser        = argparse.ArgumentParser(
        description = "Download all files from a Zenodo record "
                      "(default: Quaia catalog 10.5281/zenodo.8060755)",
        formatter_class = argparse.RawDescriptionHelpFormatter,
        epilog      = textwrap.dedent(
            """\
            Examples
            --------
            Download default Quaia record into current directory:
                python download_zenodo_quaia.py

            Download another record into ~/Downloads/quaia_data:
                python download_zenodo_quaia.py \\
                    --doi    10.5281/zenodo.8098636 \\
                    --outdir ~/Downloads/quaia_data
            """
        ),
    )

    parser.add_argument(
        "--doi",
        default  = "10.5281/zenodo.8060755",
        help     = (
            "Zenodo DOI, record URL, or numeric record ID "
            "(default: 10.5281/zenodo.8060755)."
        ),
    )

    parser.add_argument(
        "--outdir",
        default  = ".",
        help     = (
            "Output directory for downloaded files "
            "(default: current directory)."
        ),
    )

    parser.add_argument(
        "--overwrite",
        action   = "store_true",
        help     = "Overwrite existing files instead of skipping them.",
    )

    return parser


def main(argv = None) -> None:
    if argv is None:
        argv       = sys.argv[1:]

    parser        = build_parser()
    args          = parser.parse_args(argv)

    try:
        record_id  = parse_record_id(args.doi)
    except ValueError as e:
        parser.error(str(e))
        return

    print(f"[info] Zenodo record ID: {record_id}")
    print(f"[info] Output directory: {os.path.abspath(args.outdir)}")

    meta          = fetch_record_metadata(record_id)

    title         = meta.get("metadata", {}).get("title", "Unknown title")
    version       = meta.get("metadata", {}).get("version", "")
    if version:
        print(f"[info] Title   : {title}  (version {version})")
    else:
        print(f"[info] Title   : {title}")

    files         = meta.get("files", [])
    if not files:
        print("[warn] No files listed for this record.")
        return

    print(f"[info] Found {len(files)} file(s) to download.\n")

    for file_info in files:
        download_file(
            file_info,
            outdir    = args.outdir,
            overwrite = args.overwrite,
        )


if __name__ == "__main__":
    main()
