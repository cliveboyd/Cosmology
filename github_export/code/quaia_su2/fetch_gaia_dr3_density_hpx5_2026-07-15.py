from __future__ import annotations

import argparse
import json
from pathlib import Path

import pyvo


DEFAULT_OUT = Path(r"C:\Users\clive\Documents\Cosmology\external_datasets\gaia_dr3_density")
GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap"


def main() -> int:
    parser = argparse.ArgumentParser(description="Fetch a compact Gaia DR3 HEALPix source-density map.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--hpx-level", type=int, default=4)
    parser.add_argument("--g-max", type=float, default=20.5)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = args.out_dir / f"gaia_dr3_source_density_hpx{args.hpx_level}_glt{str(args.g_max).replace('.', 'p')}.csv"
    out_json = args.out_dir / f"gaia_dr3_source_density_hpx{args.hpx_level}_manifest.json"
    if out_csv.exists() and out_csv.stat().st_size > 1000:
        print(f"Existing density CSV: {out_csv}")
        return 0

    query = f"""
SELECT hpx, COUNT(*) AS n_source
FROM (
  SELECT gaia_healpix_index({args.hpx_level}, source_id) AS hpx
  FROM gaiadr3.gaia_source_lite
  WHERE phot_g_mean_mag < {args.g_max}
) AS sub
GROUP BY hpx
ORDER BY hpx
"""
    print(query.strip(), flush=True)
    service = pyvo.dal.TAPService(GAIA_TAP)
    table = service.search(query, maxrec=200000).to_table()
    table.write(out_csv, format="ascii.csv", overwrite=True)
    manifest = {
        "source": "ESA Gaia Archive TAP",
        "tap_url": GAIA_TAP,
        "table": "gaiadr3.gaia_source_lite",
        "query": query.strip(),
        "hpx_level": args.hpx_level,
        "g_max": args.g_max,
        "rows": len(table),
        "local_csv": str(out_csv),
    }
    out_json.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Saved density CSV: {out_csv}")
    print(f"Saved manifest: {out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
