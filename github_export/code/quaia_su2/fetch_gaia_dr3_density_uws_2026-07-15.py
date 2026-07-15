from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import requests


DEFAULT_OUT = Path(r"C:\Users\clive\Documents\Cosmology\external_datasets\gaia_dr3_density")
DEFAULT_TAP = "https://gaia.ari.uni-heidelberg.de/tap"


def main() -> int:
    parser = argparse.ArgumentParser(description="Fetch Gaia DR3 source density via direct TAP async/UWS.")
    parser.add_argument("--tap-url", default=DEFAULT_TAP)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--hpx-level", type=int, default=4)
    parser.add_argument("--g-max", type=float, default=20.5)
    parser.add_argument("--poll-seconds", type=float, default=15.0)
    parser.add_argument("--max-wait-seconds", type=float, default=1800.0)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = args.out_dir / f"gaia_dr3_source_density_hpx{args.hpx_level}_glt{str(args.g_max).replace('.', 'p')}.csv"
    out_manifest = args.out_dir / f"gaia_dr3_source_density_hpx{args.hpx_level}_manifest.json"

    query = f"""
SELECT hpx, COUNT(*) AS n_source
FROM (
  SELECT gaia_healpix_index({args.hpx_level}, source_id) AS hpx
  FROM gaiadr3.gaia_source_lite
  WHERE phot_g_mean_mag < {args.g_max}
) AS sub
GROUP BY hpx
ORDER BY hpx
""".strip()

    session = requests.Session()
    async_url = args.tap_url.rstrip("/") + "/async"
    submit = session.post(
        async_url,
        data={"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": query, "MAXREC": "200000"},
        allow_redirects=False,
        timeout=120,
    )
    if submit.status_code in {301, 302, 303}:
        job_url = submit.headers["Location"]
        if job_url.startswith("/"):
            job_url = args.tap_url.rstrip("/") + job_url
    else:
        submit.raise_for_status()
        job_url = submit.url
    print(f"Created TAP async job: {job_url}", flush=True)

    run = session.post(job_url.rstrip("/") + "/phase", data={"PHASE": "RUN"}, timeout=60)
    run.raise_for_status()

    start = time.time()
    while True:
        phase_resp = session.get(job_url.rstrip("/") + "/phase", timeout=60)
        phase_resp.raise_for_status()
        phase = phase_resp.text.strip()
        print(f"Gaia density job phase: {phase}", flush=True)
        if phase in {"COMPLETED", "ERROR", "ABORTED"}:
            break
        if time.time() - start > args.max_wait_seconds:
            raise TimeoutError(f"Timed out waiting for {job_url}")
        time.sleep(args.poll_seconds)

    if phase != "COMPLETED":
        error_text = session.get(job_url.rstrip("/") + "/error", timeout=60).text
        raise RuntimeError(f"TAP job ended in {phase}: {error_text[:1000]}")

    result_url = job_url.rstrip("/") + "/results/result"
    result = session.get(result_url, timeout=180)
    result.raise_for_status()
    out_csv.write_text(result.text, encoding="utf-8")
    rows = max(len(result.text.splitlines()) - 1, 0)
    manifest = {
        "source": "Gaia TAP async/UWS",
        "tap_url": args.tap_url,
        "job_url": job_url,
        "table": "gaiadr3.gaia_source_lite",
        "query": query,
        "hpx_level": args.hpx_level,
        "g_max": args.g_max,
        "rows": rows,
        "local_csv": str(out_csv),
    }
    out_manifest.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Saved density CSV: {out_csv}")
    print(f"Saved manifest: {out_manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
