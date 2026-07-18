from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(r"C:\Users\clive\Documents\Cosmology")
CODE_DIR = ROOT / "github_export" / "code" / "quaia_su2"
DOWNLOAD = CODE_DIR / "download_sdss_dr16q_v4_2026-07-18.py"
ANALYSE = CODE_DIR / "validate_sdss_dr16q_v4_cross_catalogue_2026-07-18.py"


def run(script: Path) -> None:
    command = [sys.executable, "-u", str(script)]
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Running: {' '.join(command)}", flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


def main() -> int:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Starting preregistered SDSS DR16Q v4 pipeline", flush=True)
    run(DOWNLOAD)
    run(ANALYSE)
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Pipeline complete", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
