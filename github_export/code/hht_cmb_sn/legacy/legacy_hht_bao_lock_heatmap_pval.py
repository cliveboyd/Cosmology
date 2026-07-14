#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 27 06:33:33 2025

@author: boyde
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Heatmaps of r(IMF, k_eff) and r(IMF, 1/k_eff) with permutation p-values and FDR masking.
"""

import argparse, numpy as np, pandas as pd, matplotlib.pyplot as plt
from pathlib import Path

# ------------------- args -------------------
parser                           = argparse.ArgumentParser()
parser.add_argument('--packs',   nargs='+',  required=True)
parser.add_argument('--labels',  nargs='*',  default=None)
parser.add_argument('--bao',     type=str,   required=True)
parser.add_argument('--rd',      type=float, default=147.09)
parser.add_argument('--imax',    type=int,   default=None)
parser.add_argument('--smooth',  type=int,   default=9)
parser.add_argument('--zmin',    type=float, default=0.02)
parser.add_argument('--zmax',    type=float, default=2.10)
parser.add_argument('--perms',   type=int,   default=2000, help='permutation count for 1/k test')
parser.add_argument('--alpha',   type=float, default=0.10, help='FDR q-level')
parser.add_argument('--out',     type=str,   default='plamb_runs/hht_bao_lock_real/heatmap_p')
parser.add_argument('--show',    action='store_true')
args                             = parser.parse_args()

OUTDIR                           = Path(args.out); OUTDIR.mkdir(parents=True, exist_ok=True)

# ------------------- BAO prep -------------------
bao                              = pd.read_csv(args.bao)
B_DM                             = bao[bao['kind'].str.contains('DM', case=False, na=False)]

def keff(dm_over_rd):
    Lc                           = dm_over_rd * args.rd
    return (2*np.pi) / Lc

# ------------------- helpers -------------------
def ridge_from_imf(imf_signal, smooth=0):
    r                            = np.asarray(imf_signal, float).copy()
    if smooth and (smooth % 2 == 1) and smooth >= 3:
        try:
            from scipy.signal import medfilt as _medfilt
            r                    = _medfilt(r, kernel_size=smooth)
        except Exception:
            pass
    return r

def pearson(a, b):
    m                            = np.isfinite(a) & np.isfinite(b)
    if np.count_nonzero(m) < 3:  return np.nan
    return np.corrcoef(a[m], b[m])[0, 1]

def perm_pval(a, b, nperm=2000, seed=12345):
    m                            = np.isfinite(a) & np.isfinite(b)
    a, b                         = np.asarray(a[m]), np.asarray(b[m])
    if a.size < 3:               return np.nan, np.nan
    rng                          = np.random.default_rng(seed)
    r_obs                        = pearson(a, b)
    gt                           = 0
    for _ in range(nperm):
        gt                      += (pearson(a, rng.permutation(b)) >= r_obs)
    return (gt + 1) / (nperm + 1), r_obs

def fdr_bh(pvals, alpha=0.10):
    p                            = np.asarray(pvals, float)
    idx                          = np.argsort(p)
    m                            = np.arange(1, p.size + 1)
    thresh                       = alpha * m / p.size
    ok                           = p[idx] <= thresh
    cut                          = np.max(np.where(ok)[0]) if np.any(ok) else -1
    qmask                        = np.zeros_like(p, dtype=bool)
    if cut >= 0:                 qmask[idx[:cut+1]] = True
    return qmask

# ------------------- iterate packs -------------------
heat_rk_list, heat_ri_list       = [], []
heat_p_list                      = []
pack_labels                      = []
rows                              = []

for i, npz_path in enumerate(args.packs):
    tag                          = (args.labels[i] if args.labels and i < len(args.labels)
                                    else Path(npz_path).parents[0].name)
    D                            = np.load(npz_path, allow_pickle=True)
    Ls                           = D['Ls']
    inst_freqs                   = D['inst_freqs']
    energies                     = D['energies']
    z                            = np.expm1(Ls)

    # BAO curves on HHT grid
    z_grid                       = np.linspace(max(args.zmin, z.min()), min(args.zmax, z.max()), 160)
    if len(B_DM):
        dm                       = np.interp(z_grid, B_DM['z'].to_numpy(), B_DM['value'].to_numpy(),
                                             left=np.nan, right=np.nan)
        k                        = keff(dm)
        invk                     = 1.0 / np.where(k != 0, k, np.nan)
    else:
        k                        = np.full_like(z_grid, np.nan)
        invk                     = np.full_like(z_grid, np.nan)

    k_on                         = np.interp(z, z_grid, k,    left=np.nan, right=np.nan)
    ik_on                        = np.interp(z, z_grid, invk, left=np.nan, right=np.nan)

    nimf                         = inst_freqs.shape[0]
    imax                         = nimf if args.imax is None else min(args.imax, nimf)

    rk_vals                      = np.full(imax, np.nan)
    ri_vals                      = np.full(imax, np.nan)
    p_vals                       = np.full(imax, np.nan)

    mwin                         = (z >= args.zmin) & (z <= args.zmax)

    for j in range(imax):
        ridge_j                  = ridge_from_imf(inst_freqs[j], smooth=args.smooth)
        a                        = ridge_j[mwin]
        rk_vals[j]               = pearson(a, k_on[mwin])
        p, r                     = perm_pval(a, ik_on[mwin], nperm=args.perms, seed=12345 + j)
        p_vals[j]                = p
        ri_vals[j]               = r

        rows.append({
            'pack'              : tag,
            'npz'               : npz_path,
            'IMF'               : j + 1,
            'energy'            : float(energies[j]) if j < len(energies) else np.nan,
            'r_with_k'          : rk_vals[j],
            'r_with_1_over_k'   : ri_vals[j],
            'p_perm_1_over_k'   : p_vals[j],
            'smooth'            : args.smooth,
            'zmin'              : args.zmin,
            'zmax'              : args.zmax
        })

    heat_rk_list.append(rk_vals)
    heat_ri_list.append(ri_vals)
    heat_p_list.append(p_vals)
    pack_labels.append(tag)

# ------------------- save table -------------------
df                               = pd.DataFrame(rows).sort_values(['pack', 'IMF'])
df.to_csv(OUTDIR/'bao_lock_heatmap_table_with_p.csv', index=False)
print(f"[ok] table+p -> {OUTDIR/'bao_lock_heatmap_table_with_p.csv'}")

# ------------------- FDR mask on 1/k correlations -------------------
p_concat                         = np.concatenate(heat_p_list)
mask_concat                      = fdr_bh(p_concat, alpha=args.alpha)
splits                           = np.cumsum([len(x) for x in heat_p_list])[:-1]
masks                            = np.split(mask_concat, splits)

def plot_heat(matrix_list, title, fname, mask=None):
    mat                          = np.vstack(matrix_list)
    fig, ax                      = plt.subplots(figsize=(1.2*mat.shape[1]+2.8, 1.1*mat.shape[0]+2.0))
    im                           = ax.imshow(mat, aspect='auto', cmap='coolwarm', vmin=-1.0, vmax=1.0)
    ax.set_title(title)
    ax.set_xlabel('IMF index'); ax.set_ylabel('Pack')
    ax.set_xticks(np.arange(mat.shape[1])); ax.set_xticklabels(np.arange(1, mat.shape[1]+1))
    ax.set_yticks(np.arange(mat.shape[0])); ax.set_yticklabels(pack_labels)
    cbar                         = plt.colorbar(im, ax=ax, shrink=0.85); cbar.set_label('Pearson r')
    if mask is not None:
        mcat                     = np.vstack(mask)
        yy, xx                   = np.where(mcat)
        ax.scatter(xx, yy, s=40, facecolors='none', edgecolors='k', linewidths=1.2, label='FDR pass')
        ax.legend(loc='upper right', fontsize=8)
    plt.tight_layout()
    out                          = OUTDIR/fname
    plt.savefig(out, dpi=150)
    if args.show:                plt.show()
    plt.close(fig)
    print(f"[ok] heatmap -> {out}")

plot_heat(heat_rk_list,
          f"r(IMF, k_eff) | smooth={args.smooth}, z∈[{args.zmin},{args.zmax}]",
          'heat_r_with_k.png')

plot_heat(heat_ri_list,
          f"r(IMF, 1/k_eff) | smooth={args.smooth}, z∈[{args.zmin},{args.zmax}]",
          'heat_r_with_1_over_k.png',
          mask=masks)
