

"""
quaia_plot_mc_shuffle.py  (0V2)

Visualisation & table helper for the pixel-shuffle Monte Carlo:

  - Reads quaia_mc_bcut_shuffle_summary.txt.
  - Extracts bcut, z_mid, A_obs, A_par_obs, f_par_obs, p(A), p(|A_par|).
  - Makes two plots: p(A) vs z_mid, p(|A_par|) vs z_mid, for each |b|-cut.
  - Prints a compact table suitable for notes / LaTeX.

Usage (from IPython / Spyder):

  %run /Users/boyde/.spyder-py3/quaia_plot_mc_shuffle.py
"""

from pathlib import Path
import re

import numpy             as np
import matplotlib.pyplot as plt

from   quaia_config      import OUT_DIR


# ---------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------

def parse_line_table_style(line: str):
    """
    Parse one table-style line, e.g.

    10.0  z0p10_0p30  0.20  218  4.873e-03  4.308e-03  +0.884  0.562  0.213  0.104
    """
    if not line or line.lstrip().startswith("#"):
        return None

    parts = line.split()
    
    # Expect at least: bcut, tag, z_mid, N_good, A_obs, A_par_obs, f_par_obs, pA, pAabs
    if len(parts) < 9:
        return None

    bcut      = float(parts[0])
    tag       = parts[1]
    z_mid     = float(parts[2])
    
    # N_good  = int(parts[3])   # not actually used by the plots, so we can ignore
    A_obs     = float(parts[4])
    A_par_obs = float(parts[5])
    f_par_obs = float(parts[6])
    pA        = float(parts[7])
    pAabs     = float(parts[8])
    # optional 10th column p(|f_par|) can be ignored or stored if you like

    return dict(
                bcut      = bcut,
                tag       = tag,
                z_mid     = z_mid,
                A_obs     = A_obs,
                A_par_obs = A_par_obs,
                f_par_obs = f_par_obs,
                pA        = pA,
                pAabs     = pAabs,
               )


_log_pattern = re.compile(
                          r"\[mc\]\s+\|b\|>(?P<bcut>[\d.]+),"
                          r"\s*(?P<z_lo>[\d.]+).*?z\s*<\s*(?P<z_hi>[\d.]+),"
                          r".*?A_obs=(?P<A_obs>[-+eE\d.]+),"
                          r"\s*A_par_obs=(?P<A_par_obs>[-+eE\d.]+),"
                          r"\s*f_par_obs=(?P<f_par_obs>[-+eE\d.]+),"
                          r"\s*p\(A\)=(?P<pA>[-+eE\d.]+),"
                          r"\s*p\(\|A_par\|\)=(?P<pAabs>[-+eE\d.]+)"
                         )


def parse_line_log_style(line: str):
    """
    Parse an old-style log line beginning with "[mc] |b|>...".
    """
    m = _log_pattern.search(line)
    if not m:
        return None
    
    d     = m.groupdict()
    z_lo  = float(d["z_lo"])
    z_hi  = float(d["z_hi"])
    z_mid = 0.5 * (z_lo + z_hi)

    return dict(
        bcut      = float(d["bcut"]),
        tag       = f"z{z_lo:.2f}_{z_hi:.2f}".replace(".", "p"),
        z_mid     = z_mid,
        A_obs     = float(d["A_obs"]),
        A_par_obs = float(d["A_par_obs"]),
        f_par_obs = float(d["f_par_obs"]),
        pA        = float(d["pA"]),
        pAabs     = float(d["pAabs"]),
    )


def parse_shuffle_summary(path: Path):
    """
    Try to parse the summary file in either the new *table* format or the
    old *log* format. Returns a list of dicts.
    """
    records = []
    with path.open("r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            rec = None
            # First try table-style (current format)
            rec = parse_line_table_style(line)
            if rec is None and line.startswith("[mc]"):
                # fall back to old log-style
                rec = parse_line_log_style(line)

            if rec is not None:
                records.append(rec)

    if not records:
        raise RuntimeError(f"No MC records parsed from {path}")

    return records


# ---------------------------------------------------------------------
# Plotting and tables
# ---------------------------------------------------------------------

def make_plots(records):
    """
    Make p(A) and p(|A_par|) vs z_mid plots for each bcut.
    """
    bcuts = sorted({r["bcut"] for r in records})

    # p(A) vs z
    plt.figure()
    for b in bcuts:
        sub = [r for r in records if np.isclose(r["bcut"], b)]
        sub = sorted(sub, key=lambda r: r["z_mid"])
        z   = np.array([r["z_mid"] for r in sub])
        pA  = np.array([r["pA"] for r in sub])
        plt.plot(z, pA, marker="o", label=f"|b|>{int(b)}°")
    plt.axhline(0.05, linestyle="--")
    plt.xlabel("z_mid")
    plt.ylabel("p(A)")
    plt.title("Pixel-shuffle MC: p(A) vs redshift")
    plt.legend()
    plt.grid(True, alpha=0.3)

    # p(|A_par|) vs z
    plt.figure()
    for b in bcuts:
        sub   = [r for r in records if np.isclose(r["bcut"], b)]
        sub   = sorted(sub, key=lambda r: r["z_mid"])
        z     = np.array([r["z_mid"] for r in sub])
        pAabs = np.array([r["pAabs"] for r in sub])
        plt.plot(z, pAabs, marker="o", label=f"|b|>{int(b)}°")
    plt.axhline(0.05, linestyle="--")
    plt.xlabel("z_mid")
    plt.ylabel("p(|A_par|)")
    plt.title("Pixel-shuffle MC: p(|A_par|) vs redshift")
    plt.legend()
    plt.grid(True, alpha=0.3)

    print("[plot-mc-shuffle] Plots generated (Agg backend -> check Plots pane).")


def print_compact_table(records):
    """
    Print a compact notes / LaTeX style table.
    """
    bcuts = sorted({r["bcut"] for r in records})

    print("# Pixel-shuffle MC summary")
    print("# Columns: z_mid  A_obs  A_par_obs  f_par_obs  p(A)  p(|A_par|)")

    for b in bcuts:
        print()
        print(f"=== |b|>{int(b)}° ===")
        sub = [r for r in records if np.isclose(r["bcut"], b)]
        sub = sorted(sub, key=lambda r: r["z_mid"])
        for r in sub:
            print(
                  f"{r['z_mid']:4.2f}  "
                  f"{r['A_obs']: .3e}  "
                  f"{r['A_par_obs']: .3e}  "
                  f"{r['f_par_obs']: .3f}  "
                  f"{r['pA']: .3f}  "
                  f"{r['pAabs']: .3f}"
                 )

def main():
    summary_path = OUT_DIR / "quaia_mc_bcut_shuffle_summary.txt"
    print(f"[plot-mc-shuffle] reading {summary_path}")
    records = parse_shuffle_summary(summary_path)
    print_compact_table(records)
    make_plots(records)


if __name__ == "__main__":
    main()
