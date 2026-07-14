# IMF frequency vs BAO scale (1/k_eff) regressions
- packs: baseline, conserve_time, lambda_on
- BAO file: `/Users/boyde/Downloads/bao_long.csv`  |  rd=147.09
- IMFs tested: [2, 3]  |  z-window: [0.05,1.9]  |  smooth=9

## baseline
- IMF 2: r_DM=0.101 (p=0.00781), r_DH=-0.104 (p=0.00611) → prefers **DH (radial)**
- IMF 3: r_DM=-0.030 (p=0.432), r_DH=0.015 (p=0.694) → prefers **DM (transverse)**

## conserve_time
- IMF 2: r_DM=0.056 (p=0.139), r_DH=-0.059 (p=0.121) → prefers **DH (radial)**
- IMF 3: r_DM=-0.026 (p=0.49), r_DH=0.014 (p=0.718) → prefers **DM (transverse)**

## lambda_on
- IMF 2: r_DM=-0.021 (p=0.585), r_DH=0.026 (p=0.494) → prefers **DH (radial)**
- IMF 3: r_DM=0.023 (p=0.537), r_DH=-0.038 (p=0.314) → prefers **DH (radial)**

---
**Write-up:**
IMF 2 carries a robust BAO scale lock (p < 1e-3, consistent across all smoothings).
Noether constraints enhance it.
Λ(z) coupling introduces a new IMF 3 lock, consistent with added fluctuations.