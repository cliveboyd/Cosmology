from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import requests


DOI = "10.7910/DVN/EWCNL5"
DATAVERSE = "https://dataverse.harvard.edu"
FILENAMES = ("SFD_dust_4096_ngp.fits", "SFD_dust_4096_sgp.fits")


def md5sum(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def get_metadata() -> dict:
    url = f"{DATAVERSE}/api/datasets/:persistentId"
    response = requests.get(url, params={"persistentId": f"doi:{DOI}"}, timeout=60)
    response.raise_for_status()
    return response.json()


def download_file(file_id: int, md5: str, out_path: Path) -> None:
    if out_path.exists():
        existing = md5sum(out_path)
        if existing == md5:
            print(f"Existing file OK: {out_path}")
            return
        print(f"Existing file has wrong MD5, replacing: {out_path}")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_suffix(out_path.suffix + ".part")
    url = f"{DATAVERSE}/api/access/datafile/{file_id}"
    with requests.get(url, stream=True, timeout=120) as response:
        response.raise_for_status()
        with tmp_path.open("wb") as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)

    got = md5sum(tmp_path)
    if got != md5:
        tmp_path.unlink(missing_ok=True)
        raise RuntimeError(f"MD5 mismatch for {out_path.name}: expected {md5}, got {got}")
    tmp_path.replace(out_path)
    print(f"Downloaded: {out_path}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Fetch SFD dust FITS maps from Harvard Dataverse.")
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path(r"C:\Users\clive\Documents\Cosmology\external_datasets\sfd"),
    )
    args = parser.parse_args()

    metadata = get_metadata()
    files = metadata["data"]["latestVersion"]["files"]
    by_name = {row["dataFile"]["filename"]: row["dataFile"] for row in files}

    manifest = {
        "doi": DOI,
        "dataverse_dataset": metadata["data"].get("persistentUrl", ""),
        "files": [],
    }
    for filename in FILENAMES:
        if filename not in by_name:
            raise RuntimeError(f"{filename} not found in Dataverse metadata")
        data_file = by_name[filename]
        out_path = args.out_dir / filename
        download_file(int(data_file["id"]), str(data_file["md5"]), out_path)
        manifest["files"].append(
            {
                "filename": filename,
                "dataverse_file_id": int(data_file["id"]),
                "md5": str(data_file["md5"]),
                "local_path": str(out_path),
                "bytes": out_path.stat().st_size,
            }
        )

    manifest_path = args.out_dir / "sfd_dataverse_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Saved manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
