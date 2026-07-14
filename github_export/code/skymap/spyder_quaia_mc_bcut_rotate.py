#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 22:37:32 2025

@author: boyde
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quaia_mc_bcut_rotate.py  (0V1)

Rotation-based Monte Carlo for the <z>–dipole.

For each (bcut, redshift slice) combination:
  * Load the <z> map from quaia_build_zmean_maps_bcut.py outputs.
  * Compute the observed dipole vector b_obs from <z> - < <z> >.
  * Compute its amplitude A, CMB-parallel component A_parallel, and
    fraction f_parallel = A_parallel / A.
  * Perform N_mock random ROTATIONS of the sky directions (n_hat) while
    keeping z-values fixed, and refit the dipole each time.
  * From the MC distributions, derive p-values for:
        - A (amplitude),
        - |A_parallel|,
        - |f_parallel|.
  * Save results per slice as NPZ files:

        quaia_mc_bcut_rotate_{tag}_bcut{BCUT}.npz

Fields in each NPZ:
    amps        : (N_mock,)     |b_vec| distribution
    amp_pars    : (N_mock,)     A_parallel distribution
    frac_pars   : (N_mock,)     f_parallel distribution

    obs_amp     : float         observed |b_vec|
    obs_amp_par : float         observed A_parallel
    obs_frac_par: float         observed f_parallel

    z_mid       : float         mid-point of redshift bin (NaN for full)
    label       : str           label, e.g. "0.10 ≤ z < 0.30" or "full sample"
    bcut        : float         |b| cut
    N_good_obs  : int           number of good pixels in slice

    p_amp       : float         P(A_MC >= A_obs)
    p_par_abs   : float         P(|A_parallel,MC| >= |A_parallel,obs|)
    p_frac_abs  : float         P(|f_parallel,MC| >= |f_parallel,obs|)

Usage (Spyder/IPython):

    %run /Users/boyde/.spyder-py3/quaia_mc_bcut_rotate.py --n-mock 1000 --seed 123

Optional:

    %run /Users/boyde/.spyder-py3/quaia_mc_bcut_rotate.py \
         --n-mock 2000 --seed 42 --bcut-list 20 30
         
         
RESULTS....

N_MIN_PER_PIX = 3
%run /Users/boyde/.spyder-py3/quaia_mc_bcut_rotate.py --n-mock 2000 --seed 42
[mc-rot] starting rotation MC: n_mock=2000, seed=42
[mc-rot] bcut_list = [10.0, 15.0, 20.0, 25.0, 30.0, 35.0]
[mc-rot] OUT_DIR = /Users/boyde/.spyder-py3/skymap2/quaia_outputs
[mc-rot] |b|>10.0, tag=full, label='full sample', z_mid=nan, N_good=35572
[mc-rot]   obs A=5.216e-03, A_par=2.532e-05, f_par=+0.005
[mc-rot]   p(A)=0.038, p(|A_par|)=0.995, p(|f_par|)=0.995
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_full_bcut10.npz
[mc-rot] |b|>10.0, tag=z0p10_0p30, label='0.10 ≤ z < 0.30', z_mid=0.20, N_good=218
[mc-rot]   obs A=4.873e-03, A_par=4.308e-03, f_par=+0.884
[mc-rot]   p(A)=0.026, p(|A_par|)=0.133, p(|f_par|)=0.133
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p10_0p30_bcut10.npz
[mc-rot] |b|>10.0, tag=z0p30_0p50, label='0.30 ≤ z < 0.50', z_mid=0.40, N_good=3307
[mc-rot]   obs A=2.096e-03, A_par=-1.676e-03, f_par=-0.799
[mc-rot]   p(A)=0.232, p(|A_par|)=0.205, p(|f_par|)=0.205
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p30_0p50_bcut10.npz
[mc-rot] |b|>10.0, tag=z0p50_0p75, label='0.50 ≤ z < 0.75', z_mid=0.62, N_good=7916
[mc-rot]   obs A=1.827e-03, A_par=-9.835e-04, f_par=-0.538
[mc-rot]   p(A)=0.574, p(|A_par|)=0.457, p(|f_par|)=0.457
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p50_0p75_bcut10.npz
[mc-rot] |b|>10.0, tag=z0p75_1p00, label='0.75 ≤ z < 1.00', z_mid=0.88, N_good=13075
[mc-rot]   obs A=9.801e-04, A_par=-5.151e-04, f_par=-0.526
[mc-rot]   p(A)=0.905, p(|A_par|)=0.465, p(|f_par|)=0.465
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p75_1p00_bcut10.npz
[mc-rot] |b|>10.0, tag=z1p00_1p50, label='1.00 ≤ z < 1.50', z_mid=1.25, N_good=29498
[mc-rot]   obs A=1.633e-03, A_par=-6.599e-04, f_par=-0.404
[mc-rot]   p(A)=0.311, p(|A_par|)=0.602, p(|f_par|)=0.602
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z1p00_1p50_bcut10.npz
[mc-rot] |b|>10.0, tag=z1p50_2p50, label='1.50 ≤ z < 2.50', z_mid=2.00, N_good=32791
[mc-rot]   obs A=1.549e-03, A_par=-1.294e-03, f_par=-0.836
[mc-rot]   p(A)=0.742, p(|A_par|)=0.164, p(|f_par|)=0.164
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z1p50_2p50_bcut10.npz
[mc-rot] |b|>15.0, tag=full, label='full sample', z_mid=nan, N_good=32877
[mc-rot]   obs A=5.750e-03, A_par=-7.689e-04, f_par=-0.134
[mc-rot]   p(A)=0.009, p(|A_par|)=0.866, p(|f_par|)=0.866
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_full_bcut15.npz
[mc-rot] |b|>15.0, tag=z0p10_0p30, label='0.10 ≤ z < 0.30', z_mid=0.20, N_good=211
[mc-rot]   obs A=5.645e-03, A_par=4.449e-03, f_par=+0.788
[mc-rot]   p(A)=0.138, p(|A_par|)=0.206, p(|f_par|)=0.206
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p10_0p30_bcut15.npz
[mc-rot] |b|>15.0, tag=z0p30_0p50, label='0.30 ≤ z < 0.50', z_mid=0.40, N_good=3215
[mc-rot]   obs A=1.843e-03, A_par=-1.633e-03, f_par=-0.886
[mc-rot]   p(A)=0.911, p(|A_par|)=0.117, p(|f_par|)=0.117
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p30_0p50_bcut15.npz
[mc-rot] |b|>15.0, tag=z0p50_0p75, label='0.50 ≤ z < 0.75', z_mid=0.62, N_good=7663
[mc-rot]   obs A=1.830e-03, A_par=-9.627e-04, f_par=-0.526
[mc-rot]   p(A)=0.133, p(|A_par|)=0.471, p(|f_par|)=0.471
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p50_0p75_bcut15.npz
[mc-rot] |b|>15.0, tag=z0p75_1p00, label='0.75 ≤ z < 1.00', z_mid=0.88, N_good=12628
[mc-rot]   obs A=1.189e-03, A_par=-3.977e-04, f_par=-0.335
[mc-rot]   p(A)=0.014, p(|A_par|)=0.657, p(|f_par|)=0.657
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p75_1p00_bcut15.npz
[mc-rot] |b|>15.0, tag=z1p00_1p50, label='1.00 ≤ z < 1.50', z_mid=1.25, N_good=28108
[mc-rot]   obs A=1.947e-03, A_par=-5.775e-04, f_par=-0.297
[mc-rot]   p(A)=0.272, p(|A_par|)=0.725, p(|f_par|)=0.725
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z1p00_1p50_bcut15.npz
[mc-rot] |b|>15.0, tag=z1p50_2p50, label='1.50 ≤ z < 2.50', z_mid=2.00, N_good=31080
[mc-rot]   obs A=1.771e-03, A_par=-1.524e-03, f_par=-0.860
[mc-rot]   p(A)=0.046, p(|A_par|)=0.135, p(|f_par|)=0.135
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z1p50_2p50_bcut15.npz
[mc-rot] |b|>20.0, tag=full, label='full sample', z_mid=nan, N_good=29738
[mc-rot]   obs A=5.518e-03, A_par=-1.199e-04, f_par=-0.022
[mc-rot]   p(A)=0.997, p(|A_par|)=0.978, p(|f_par|)=0.978
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_full_bcut20.npz
[mc-rot] |b|>20.0, tag=z0p10_0p30, label='0.10 ≤ z < 0.30', z_mid=0.20, N_good=199
[mc-rot]   obs A=8.057e-03, A_par=6.556e-03, f_par=+0.814
[mc-rot]   p(A)=0.301, p(|A_par|)=0.172, p(|f_par|)=0.172
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p10_0p30_bcut20.npz
[mc-rot] |b|>20.0, tag=z0p30_0p50, label='0.30 ≤ z < 0.50', z_mid=0.40, N_good=3011
[mc-rot]   obs A=1.681e-03, A_par=-1.453e-03, f_par=-0.865
[mc-rot]   p(A)=0.598, p(|A_par|)=0.142, p(|f_par|)=0.142
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p30_0p50_bcut20.npz
[mc-rot] |b|>20.0, tag=z0p50_0p75, label='0.50 ≤ z < 0.75', z_mid=0.62, N_good=7179
[mc-rot]   obs A=1.416e-03, A_par=-6.466e-04, f_par=-0.457
[mc-rot]   p(A)=0.889, p(|A_par|)=0.541, p(|f_par|)=0.541
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p50_0p75_bcut20.npz
[mc-rot] |b|>20.0, tag=z0p75_1p00, label='0.75 ≤ z < 1.00', z_mid=0.88, N_good=11833
[mc-rot]   obs A=1.524e-03, A_par=-3.438e-04, f_par=-0.226
[mc-rot]   p(A)=0.641, p(|A_par|)=0.775, p(|f_par|)=0.775
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p75_1p00_bcut20.npz
[mc-rot] |b|>20.0, tag=z1p00_1p50, label='1.00 ≤ z < 1.50', z_mid=1.25, N_good=25934
[mc-rot]   obs A=1.587e-03, A_par=-4.768e-04, f_par=-0.300
[mc-rot]   p(A)=0.829, p(|A_par|)=0.703, p(|f_par|)=0.703
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z1p00_1p50_bcut20.npz
[mc-rot] |b|>20.0, tag=z1p50_2p50, label='1.50 ≤ z < 2.50', z_mid=2.00, N_good=28456
[mc-rot]   obs A=2.005e-03, A_par=-1.758e-03, f_par=-0.877
[mc-rot]   p(A)=0.973, p(|A_par|)=0.131, p(|f_par|)=0.131
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z1p50_2p50_bcut20.npz
[mc-rot] |b|>25.0, tag=full, label='full sample', z_mid=nan, N_good=26335
[mc-rot]   obs A=5.349e-03, A_par=-5.003e-04, f_par=-0.094
[mc-rot]   p(A)=0.998, p(|A_par|)=0.903, p(|f_par|)=0.903
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_full_bcut25.npz
[mc-rot] |b|>25.0, tag=z0p10_0p30, label='0.10 ≤ z < 0.30', z_mid=0.20, N_good=180
[mc-rot]   obs A=6.449e-03, A_par=5.271e-03, f_par=+0.817
[mc-rot]   p(A)=0.168, p(|A_par|)=0.177, p(|f_par|)=0.177
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p10_0p30_bcut25.npz
[mc-rot] |b|>25.0, tag=z0p30_0p50, label='0.30 ≤ z < 0.50', z_mid=0.40, N_good=2758
[mc-rot]   obs A=1.998e-03, A_par=-1.814e-03, f_par=-0.908
[mc-rot]   p(A)=0.660, p(|A_par|)=0.092, p(|f_par|)=0.092
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p30_0p50_bcut25.npz
[mc-rot] |b|>25.0, tag=z0p50_0p75, label='0.50 ≤ z < 0.75', z_mid=0.62, N_good=6558
[mc-rot]   obs A=1.179e-03, A_par=-5.058e-04, f_par=-0.429
[mc-rot]   p(A)=0.072, p(|A_par|)=0.578, p(|f_par|)=0.578
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p50_0p75_bcut25.npz
[mc-rot] |b|>25.0, tag=z0p75_1p00, label='0.75 ≤ z < 1.00', z_mid=0.88, N_good=10779
[mc-rot]   obs A=1.304e-03, A_par=-2.187e-05, f_par=-0.017
[mc-rot]   p(A)=0.683, p(|A_par|)=0.982, p(|f_par|)=0.982
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p75_1p00_bcut25.npz
[mc-rot] |b|>25.0, tag=z1p00_1p50, label='1.00 ≤ z < 1.50', z_mid=1.25, N_good=23266
[mc-rot]   obs A=2.482e-03, A_par=-8.125e-04, f_par=-0.327
[mc-rot]   p(A)=0.868, p(|A_par|)=0.670, p(|f_par|)=0.670
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z1p00_1p50_bcut25.npz
[mc-rot] |b|>25.0, tag=z1p50_2p50, label='1.50 ≤ z < 2.50', z_mid=2.00, N_good=25388
[mc-rot]   obs A=2.499e-03, A_par=-2.058e-03, f_par=-0.823
[mc-rot]   p(A)=0.332, p(|A_par|)=0.203, p(|f_par|)=0.203
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z1p50_2p50_bcut25.npz
[mc-rot] |b|>30.0, tag=full, label='full sample', z_mid=nan, N_good=22999
[mc-rot]   obs A=5.038e-03, A_par=4.095e-04, f_par=+0.081
[mc-rot]   p(A)=0.252, p(|A_par|)=0.918, p(|f_par|)=0.918
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_full_bcut30.npz
[mc-rot] |b|>30.0, tag=z0p10_0p30, label='0.10 ≤ z < 0.30', z_mid=0.20, N_good=157
[mc-rot]   obs A=5.475e-03, A_par=3.910e-03, f_par=+0.714
[mc-rot]   p(A)=0.329, p(|A_par|)=0.277, p(|f_par|)=0.277
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p10_0p30_bcut30.npz
[mc-rot] |b|>30.0, tag=z0p30_0p50, label='0.30 ≤ z < 0.50', z_mid=0.40, N_good=2457
[mc-rot]   obs A=2.171e-03, A_par=-2.041e-03, f_par=-0.940
[mc-rot]   p(A)=0.477, p(|A_par|)=0.061, p(|f_par|)=0.061
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p30_0p50_bcut30.npz
[mc-rot] |b|>30.0, tag=z0p50_0p75, label='0.50 ≤ z < 0.75', z_mid=0.62, N_good=5913
[mc-rot]   obs A=1.096e-03, A_par=-4.752e-04, f_par=-0.433
[mc-rot]   p(A)=0.724, p(|A_par|)=0.579, p(|f_par|)=0.579
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p50_0p75_bcut30.npz
[mc-rot] |b|>30.0, tag=z0p75_1p00, label='0.75 ≤ z < 1.00', z_mid=0.88, N_good=9664
[mc-rot]   obs A=1.278e-03, A_par=-2.486e-04, f_par=-0.194
[mc-rot]   p(A)=0.504, p(|A_par|)=0.806, p(|f_par|)=0.806
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p75_1p00_bcut30.npz
[mc-rot] |b|>30.0, tag=z1p00_1p50, label='1.00 ≤ z < 1.50', z_mid=1.25, N_good=20546
[mc-rot]   obs A=2.647e-03, A_par=-4.280e-04, f_par=-0.162
[mc-rot]   p(A)=0.134, p(|A_par|)=0.832, p(|f_par|)=0.832
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z1p00_1p50_bcut30.npz
[mc-rot] |b|>30.0, tag=z1p50_2p50, label='1.50 ≤ z < 2.50', z_mid=2.00, N_good=22326
[mc-rot]   obs A=1.816e-03, A_par=-1.684e-03, f_par=-0.927
[mc-rot]   p(A)=0.251, p(|A_par|)=0.065, p(|f_par|)=0.065
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z1p50_2p50_bcut30.npz
[mc-rot] |b|>35.0, tag=full, label='full sample', z_mid=nan, N_good=19810
[mc-rot]   obs A=4.797e-03, A_par=-4.773e-04, f_par=-0.100
[mc-rot]   p(A)=0.929, p(|A_par|)=0.898, p(|f_par|)=0.898
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_full_bcut35.npz
[mc-rot] |b|>35.0, tag=z0p10_0p30, label='0.10 ≤ z < 0.30', z_mid=0.20, N_good=128
[mc-rot]   obs A=5.397e-03, A_par=4.975e-03, f_par=+0.922
[mc-rot]   p(A)=0.851, p(|A_par|)=0.080, p(|f_par|)=0.080
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p10_0p30_bcut35.npz
[mc-rot] |b|>35.0, tag=z0p30_0p50, label='0.30 ≤ z < 0.50', z_mid=0.40, N_good=2165
[mc-rot]   obs A=3.127e-03, A_par=-2.438e-03, f_par=-0.780
[mc-rot]   p(A)=0.339, p(|A_par|)=0.225, p(|f_par|)=0.225
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p30_0p50_bcut35.npz
[mc-rot] |b|>35.0, tag=z0p50_0p75, label='0.50 ≤ z < 0.75', z_mid=0.62, N_good=5252
[mc-rot]   obs A=1.427e-03, A_par=-8.762e-04, f_par=-0.614
[mc-rot]   p(A)=0.616, p(|A_par|)=0.395, p(|f_par|)=0.395
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p50_0p75_bcut35.npz
[mc-rot] |b|>35.0, tag=z0p75_1p00, label='0.75 ≤ z < 1.00', z_mid=0.88, N_good=8509
[mc-rot]   obs A=1.765e-03, A_par=1.115e-05, f_par=+0.006
[mc-rot]   p(A)=0.260, p(|A_par|)=0.994, p(|f_par|)=0.994
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p75_1p00_bcut35.npz
[mc-rot] |b|>35.0, tag=z1p00_1p50, label='1.00 ≤ z < 1.50', z_mid=1.25, N_good=17867
[mc-rot]   obs A=3.326e-03, A_par=-1.834e-04, f_par=-0.055
[mc-rot]   p(A)=0.545, p(|A_par|)=0.947, p(|f_par|)=0.947
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z1p00_1p50_bcut35.npz
[mc-rot] |b|>35.0, tag=z1p50_2p50, label='1.50 ≤ z < 2.50', z_mid=2.00, N_good=19303
[mc-rot]   obs A=2.236e-03, A_par=-1.684e-03, f_par=-0.753
[mc-rot]   p(A)=0.421, p(|A_par|)=0.247, p(|f_par|)=0.247
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z1p50_2p50_bcut35.npz
[mc-rot] done.

RESULTS....

N_MIN_PER_PIX = 5
%run /Users/boyde/.spyder-py3/quaia_mc_bcut_rotate.py --n-mock 2000 --seed 42
[mc-rot] starting rotation MC: n_mock=2000, seed=42
[mc-rot] bcut_list = [10.0, 15.0, 20.0, 25.0, 30.0, 35.0]
[mc-rot] OUT_DIR = /Users/boyde/.spyder-py3/skymap2/quaia_outputs
[mc-rot] |b|>10.0, tag=full, label='full sample', z_mid=nan, N_good=34960
[mc-rot]   obs A=5.003e-03, A_par=-1.836e-04, f_par=-0.037
[mc-rot]   p(A)=0.860, p(|A_par|)=0.969, p(|f_par|)=0.969
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_full_bcut10.npz
[mc-rot] |b|>10.0, tag=z0p10_0p30, label='0.10 ≤ z < 0.30', z_mid=0.20, N_good=3
[mc-rot]   obs A=5.178e-02, A_par=-9.832e-03, f_par=-0.190
[mc-rot]   p(A)=0.264, p(|A_par|)=0.679, p(|f_par|)=0.806
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p10_0p30_bcut10.npz
[mc-rot] |b|>10.0, tag=z0p30_0p50, label='0.30 ≤ z < 0.50', z_mid=0.40, N_good=195
[mc-rot]   obs A=4.767e-03, A_par=-1.065e-03, f_par=-0.224
[mc-rot]   p(A)=0.628, p(|A_par|)=0.783, p(|f_par|)=0.783
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p30_0p50_bcut10.npz
[mc-rot] |b|>10.0, tag=z0p50_0p75, label='0.50 ≤ z < 0.75', z_mid=0.62, N_good=1069
[mc-rot]   obs A=3.817e-03, A_par=1.514e-03, f_par=+0.397
[mc-rot]   p(A)=0.001, p(|A_par|)=0.605, p(|f_par|)=0.605
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p50_0p75_bcut10.npz
[mc-rot] |b|>10.0, tag=z0p75_1p00, label='0.75 ≤ z < 1.00', z_mid=0.88, N_good=2868
[mc-rot]   obs A=1.169e-03, A_par=5.315e-04, f_par=+0.455
[mc-rot]   p(A)=0.515, p(|A_par|)=0.541, p(|f_par|)=0.541
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p75_1p00_bcut10.npz
[mc-rot] |b|>10.0, tag=z1p00_1p50, label='1.00 ≤ z < 1.50', z_mid=1.25, N_good=19703
[mc-rot]   obs A=1.971e-03, A_par=-1.150e-03, f_par=-0.584
[mc-rot]   p(A)=0.356, p(|A_par|)=0.431, p(|f_par|)=0.431
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z1p00_1p50_bcut10.npz
[mc-rot] |b|>10.0, tag=z1p50_2p50, label='1.50 ≤ z < 2.50', z_mid=2.00, N_good=27815
[mc-rot]   obs A=1.667e-03, A_par=-1.496e-03, f_par=-0.897
[mc-rot]   p(A)=0.700, p(|A_par|)=0.100, p(|f_par|)=0.100
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z1p50_2p50_bcut10.npz
[mc-rot] |b|>15.0, tag=full, label='full sample', z_mid=nan, N_good=32554
[mc-rot]   obs A=5.665e-03, A_par=-1.031e-03, f_par=-0.182
[mc-rot]   p(A)=0.974, p(|A_par|)=0.807, p(|f_par|)=0.807
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_full_bcut15.npz
[mc-rot] |b|>15.0, tag=z0p10_0p30, label='0.10 ≤ z < 0.30', z_mid=0.20, N_good=3
[mc-rot]   obs A=5.178e-02, A_par=-9.832e-03, f_par=-0.190
[mc-rot]   p(A)=0.272, p(|A_par|)=0.695, p(|f_par|)=0.817
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p10_0p30_bcut15.npz
[mc-rot] |b|>15.0, tag=z0p30_0p50, label='0.30 ≤ z < 0.50', z_mid=0.40, N_good=190
[mc-rot]   obs A=3.689e-03, A_par=-9.038e-04, f_par=-0.245
[mc-rot]   p(A)=0.495, p(|A_par|)=0.748, p(|f_par|)=0.748
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p30_0p50_bcut15.npz
[mc-rot] |b|>15.0, tag=z0p50_0p75, label='0.50 ≤ z < 0.75', z_mid=0.62, N_good=1054
[mc-rot]   obs A=3.782e-03, A_par=1.188e-03, f_par=+0.314
[mc-rot]   p(A)=0.045, p(|A_par|)=0.690, p(|f_par|)=0.690
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p50_0p75_bcut15.npz
[mc-rot] |b|>15.0, tag=z0p75_1p00, label='0.75 ≤ z < 1.00', z_mid=0.88, N_good=2811
[mc-rot]   obs A=1.208e-03, A_par=5.827e-04, f_par=+0.483
[mc-rot]   p(A)=0.184, p(|A_par|)=0.511, p(|f_par|)=0.511
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p75_1p00_bcut15.npz
[mc-rot] |b|>15.0, tag=z1p00_1p50, label='1.00 ≤ z < 1.50', z_mid=1.25, N_good=19193
[mc-rot]   obs A=2.346e-03, A_par=-1.057e-03, f_par=-0.450
[mc-rot]   p(A)=0.386, p(|A_par|)=0.545, p(|f_par|)=0.545
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z1p00_1p50_bcut15.npz
[mc-rot] |b|>15.0, tag=z1p50_2p50, label='1.50 ≤ z < 2.50', z_mid=2.00, N_good=26918
[mc-rot]   obs A=1.832e-03, A_par=-1.743e-03, f_par=-0.951
[mc-rot]   p(A)=0.741, p(|A_par|)=0.045, p(|f_par|)=0.045
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z1p50_2p50_bcut15.npz
[mc-rot] |b|>20.0, tag=full, label='full sample', z_mid=nan, N_good=29533
[mc-rot]   obs A=5.442e-03, A_par=-3.676e-04, f_par=-0.068
[mc-rot]   p(A)=0.260, p(|A_par|)=0.929, p(|f_par|)=0.929
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_full_bcut20.npz
[mc-rot] tag=z0p10_0p30, |b|>20.0: too few good pixels (2), skipping
[mc-rot] |b|>20.0, tag=z0p30_0p50, label='0.30 ≤ z < 0.50', z_mid=0.40, N_good=181
[mc-rot]   obs A=3.368e-03, A_par=-9.428e-07, f_par=-0.000
[mc-rot]   p(A)=0.184, p(|A_par|)=1.000, p(|f_par|)=1.000
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p30_0p50_bcut20.npz
[mc-rot] |b|>20.0, tag=z0p50_0p75, label='0.50 ≤ z < 0.75', z_mid=0.62, N_good=1010
[mc-rot]   obs A=4.502e-03, A_par=6.625e-04, f_par=+0.147
[mc-rot]   p(A)=0.009, p(|A_par|)=0.857, p(|f_par|)=0.857
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p50_0p75_bcut20.npz
[mc-rot] |b|>20.0, tag=z0p75_1p00, label='0.75 ≤ z < 1.00', z_mid=0.88, N_good=2692
[mc-rot]   obs A=1.099e-03, A_par=5.279e-04, f_par=+0.481
[mc-rot]   p(A)=0.810, p(|A_par|)=0.539, p(|f_par|)=0.539
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p75_1p00_bcut20.npz
[mc-rot] |b|>20.0, tag=z1p00_1p50, label='1.00 ≤ z < 1.50', z_mid=1.25, N_good=18090
[mc-rot]   obs A=1.976e-03, A_par=-8.562e-04, f_par=-0.433
[mc-rot]   p(A)=0.889, p(|A_par|)=0.575, p(|f_par|)=0.575
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z1p00_1p50_bcut20.npz
[mc-rot] |b|>20.0, tag=z1p50_2p50, label='1.50 ≤ z < 2.50', z_mid=2.00, N_good=25107
[mc-rot]   obs A=2.135e-03, A_par=-1.962e-03, f_par=-0.919
[mc-rot]   p(A)=0.974, p(|A_par|)=0.082, p(|f_par|)=0.082
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z1p50_2p50_bcut20.npz
[mc-rot] |b|>25.0, tag=full, label='full sample', z_mid=nan, N_good=26159
[mc-rot]   obs A=5.343e-03, A_par=-6.870e-04, f_par=-0.129
[mc-rot]   p(A)=0.470, p(|A_par|)=0.870, p(|f_par|)=0.870
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_full_bcut25.npz
[mc-rot] tag=z0p10_0p30, |b|>25.0: too few good pixels (2), skipping
[mc-rot] |b|>25.0, tag=z0p30_0p50, label='0.30 ≤ z < 0.50', z_mid=0.40, N_good=174
[mc-rot]   obs A=5.065e-03, A_par=-1.390e-03, f_par=-0.274
[mc-rot]   p(A)=0.296, p(|A_par|)=0.744, p(|f_par|)=0.744
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p30_0p50_bcut25.npz
[mc-rot] |b|>25.0, tag=z0p50_0p75, label='0.50 ≤ z < 0.75', z_mid=0.62, N_good=953
[mc-rot]   obs A=4.772e-03, A_par=1.515e-03, f_par=+0.318
[mc-rot]   p(A)=0.002, p(|A_par|)=0.697, p(|f_par|)=0.697
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p50_0p75_bcut25.npz
[mc-rot] |b|>25.0, tag=z0p75_1p00, label='0.75 ≤ z < 1.00', z_mid=0.88, N_good=2482
[mc-rot]   obs A=1.443e-03, A_par=3.728e-04, f_par=+0.258
[mc-rot]   p(A)=0.941, p(|A_par|)=0.732, p(|f_par|)=0.732
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p75_1p00_bcut25.npz
[mc-rot] |b|>25.0, tag=z1p00_1p50, label='1.00 ≤ z < 1.50', z_mid=1.25, N_good=16538
[mc-rot]   obs A=2.542e-03, A_par=-8.229e-04, f_par=-0.324
[mc-rot]   p(A)=0.343, p(|A_par|)=0.693, p(|f_par|)=0.693
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z1p00_1p50_bcut25.npz
[mc-rot] |b|>25.0, tag=z1p50_2p50, label='1.50 ≤ z < 2.50', z_mid=2.00, N_good=22711
[mc-rot]   obs A=2.509e-03, A_par=-2.121e-03, f_par=-0.845
[mc-rot]   p(A)=0.578, p(|A_par|)=0.141, p(|f_par|)=0.141
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z1p50_2p50_bcut25.npz
[mc-rot] |b|>30.0, tag=full, label='full sample', z_mid=nan, N_good=22890
[mc-rot]   obs A=4.886e-03, A_par=2.751e-04, f_par=+0.056
[mc-rot]   p(A)=0.983, p(|A_par|)=0.946, p(|f_par|)=0.946
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_full_bcut30.npz
[mc-rot] tag=z0p10_0p30, |b|>30.0: too few good pixels (1), skipping
[mc-rot] |b|>30.0, tag=z0p30_0p50, label='0.30 ≤ z < 0.50', z_mid=0.40, N_good=149
[mc-rot]   obs A=1.697e-03, A_par=3.077e-04, f_par=+0.181
[mc-rot]   p(A)=0.107, p(|A_par|)=0.807, p(|f_par|)=0.807
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p30_0p50_bcut30.npz
[mc-rot] |b|>30.0, tag=z0p50_0p75, label='0.50 ≤ z < 0.75', z_mid=0.62, N_good=877
[mc-rot]   obs A=4.123e-03, A_par=1.478e-03, f_par=+0.358
[mc-rot]   p(A)=0.083, p(|A_par|)=0.654, p(|f_par|)=0.654
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p50_0p75_bcut30.npz
[mc-rot] |b|>30.0, tag=z0p75_1p00, label='0.75 ≤ z < 1.00', z_mid=0.88, N_good=2263
[mc-rot]   obs A=1.130e-03, A_par=1.018e-03, f_par=+0.901
[mc-rot]   p(A)=1.000, p(|A_par|)=0.095, p(|f_par|)=0.095
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p75_1p00_bcut30.npz
[mc-rot] |b|>30.0, tag=z1p00_1p50, label='1.00 ≤ z < 1.50', z_mid=1.25, N_good=14835
[mc-rot]   obs A=2.603e-03, A_par=-6.149e-04, f_par=-0.236
[mc-rot]   p(A)=0.585, p(|A_par|)=0.769, p(|f_par|)=0.769
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z1p00_1p50_bcut30.npz
[mc-rot] |b|>30.0, tag=z1p50_2p50, label='1.50 ≤ z < 2.50', z_mid=2.00, N_good=20166
[mc-rot]   obs A=1.448e-03, A_par=-1.408e-03, f_par=-0.973
[mc-rot]   p(A)=0.038, p(|A_par|)=0.023, p(|f_par|)=0.023
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z1p50_2p50_bcut30.npz
[mc-rot] |b|>35.0, tag=full, label='full sample', z_mid=nan, N_good=19717
[mc-rot]   obs A=4.607e-03, A_par=-4.195e-04, f_par=-0.091
[mc-rot]   p(A)=0.401, p(|A_par|)=0.905, p(|f_par|)=0.905
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_full_bcut35.npz
[mc-rot] tag=z0p10_0p30, |b|>35.0: too few good pixels (1), skipping
[mc-rot] |b|>35.0, tag=z0p30_0p50, label='0.30 ≤ z < 0.50', z_mid=0.40, N_good=134
[mc-rot]   obs A=2.055e-03, A_par=1.003e-03, f_par=+0.488
[mc-rot]   p(A)=0.059, p(|A_par|)=0.501, p(|f_par|)=0.501
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p30_0p50_bcut35.npz
[mc-rot] |b|>35.0, tag=z0p50_0p75, label='0.50 ≤ z < 0.75', z_mid=0.62, N_good=805
[mc-rot]   obs A=3.923e-03, A_par=1.087e-03, f_par=+0.277
[mc-rot]   p(A)=0.049, p(|A_par|)=0.709, p(|f_par|)=0.709
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p50_0p75_bcut35.npz
[mc-rot] |b|>35.0, tag=z0p75_1p00, label='0.75 ≤ z < 1.00', z_mid=0.88, N_good=2032
[mc-rot]   obs A=1.075e-03, A_par=1.006e-03, f_par=+0.936
[mc-rot]   p(A)=0.364, p(|A_par|)=0.064, p(|f_par|)=0.064
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z0p75_1p00_bcut35.npz
[mc-rot] |b|>35.0, tag=z1p00_1p50, label='1.00 ≤ z < 1.50', z_mid=1.25, N_good=13014
[mc-rot]   obs A=3.405e-03, A_par=-1.483e-04, f_par=-0.044
[mc-rot]   p(A)=0.216, p(|A_par|)=0.959, p(|f_par|)=0.959
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z1p00_1p50_bcut35.npz
[mc-rot] |b|>35.0, tag=z1p50_2p50, label='1.50 ≤ z < 2.50', z_mid=2.00, N_good=17588
[mc-rot]   obs A=1.690e-03, A_par=-1.420e-03, f_par=-0.840
[mc-rot]   p(A)=0.356, p(|A_par|)=0.145, p(|f_par|)=0.145
[mc-rot]   wrote -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_rotate_z1p50_2p50_bcut35.npz


"""

import argparse
from   pathlib          import Path

import numpy            as     np
import healpy           as     hp

from   quaia_config     import OUT_DIR, NSIDE, Z_BINS
from   quaia_dipole_fit import fit_dipole_linear, dipole_summary_from_b


# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------

BCUT_LIST     = [10.0, 15.0, 20.0, 25.0, 30.0, 35.0]

N_MIN_PER_PIX = 5    # Was =3
  
# CMB dipole direction in Galactic coords
L_CMB_DEG     = 264.0
B_CMB_DEG     = 48.0


# ----------------------------------------------------------------------
# Basic helpers (copied from / consistent with bcut_grid)
# ----------------------------------------------------------------------

def unit_vec_gal(l_deg, b_deg):
    """
    Unit vector in Galactic coords (l, b in degrees).
    """
    l       = np.deg2rad(l_deg)
    b       = np.deg2rad(b_deg)
    x       = np.cos(b) * np.cos(l)
    y       = np.cos(b) * np.sin(l)
    z       = np.sin(b)
    return np.array([x, y, z])


def project_dipole_on_cmb(b_vec):
    """
    Given dipole vector b_vec (shape (3,)), return:
        amp      = |b_vec|
        amp_par  = component of b_vec along CMB direction (signed)
        frac_par = amp_par / amp  (0 if amp==0)
    """
    b_vec    = np.asarray(b_vec, dtype=float)
    amp      = np.linalg.norm(b_vec)
    if amp == 0.0:
        return 0.0, 0.0, 0.0

    n_cmb    = unit_vec_gal(L_CMB_DEG, B_CMB_DEG)
    amp_par  = float(np.dot(b_vec, n_cmb))
    frac_par = amp_par / amp
    return amp, amp_par, frac_par


def load_map(tag, bcut):
    """
    Load N_map and zmean map for given tag and bcut from OUT_DIR.
    """
    path = OUT_DIR / f"quaia_zmean_{tag}_bcut{int(bcut)}_nside{NSIDE}.npz"
    if not path.exists():
        raise FileNotFoundError(f"Map not found for tag={tag}, bcut={bcut}: {path}")
    data = np.load(path, allow_pickle=False)
    return data["N"], data["zmean"]


def random_rotation_matrix(rng):
    """
    Generate a random 3x3 rotation matrix (Haar measure on SO(3))
    using random unit quaternions.
    """
    u1, u2, u3 = rng.random(3)

    q1 = np.sqrt(1.0 - u1) * np.sin(2.0 * np.pi * u2)
    q2 = np.sqrt(1.0 - u1) * np.cos(2.0 * np.pi * u2)
    q3 = np.sqrt(u1)       * np.sin(2.0 * np.pi * u3)
    q4 = np.sqrt(u1)       * np.cos(2.0 * np.pi * u3)

    # Quaternion (q1,q2,q3,q4) -> rotation matrix
    R = np.array([
        [1 - 2*(q2*q2 + q3*q3),  2*(q1*q2 - q3*q4),     2*(q1*q3 + q2*q4)],
        [2*(q1*q2 + q3*q4),      1 - 2*(q1*q1 + q3*q3), 2*(q2*q3 - q1*q4)],
        [2*(q1*q3 - q2*q4),      2*(q2*q3 + q1*q4),     1 - 2*(q1*q1 + q2*q2)],
    ])
    return R


# ----------------------------------------------------------------------
# Core per-slice analysis
# ----------------------------------------------------------------------

def prepare_slice(bcut, tag, label):
    """
    Load map for (bcut, tag), construct:
      * mask of good pixels (N>=N_MIN_PER_PIX & finite z)
      * unit vectors n_hat (N_good x 3)
      * weights (N_good,)
      * fluctuation field f = z_good - mean(z_good)
      * observed dipole stats

    Returns a dict or None if too few good pixels.
    """
    N_map, zmean = load_map(tag, bcut)

    good         = (N_map >= N_MIN_PER_PIX) & np.isfinite(zmean)
    pix          = np.where(good)[0]
    if pix.size < 3:
        print(f"[mc-rot] tag={tag}, |b|>{bcut:.1f}: too few good pixels ({pix.size}), skipping")
        return None

    N_good       = N_map[good].astype(float)
    z_good       = zmean[good].astype(float)
    f            = z_good - z_good.mean()

    theta, phi   = hp.pix2ang(NSIDE, pix)
    x            = np.sin(theta) * np.cos(phi)
    y            = np.sin(theta) * np.sin(phi)
    z            = np.cos(theta)
    n_hat        = np.column_stack([x, y, z])   # (N_good, 3)

    # Observed dipole
    a0_obs,   b_vec_obs, _              = fit_dipole_linear(f, n_hat, w=N_good)
    amp_obs , l_deg_obs, b_deg_obs      = dipole_summary_from_b(b_vec_obs)
    amp_obs2, amp_par_obs, frac_par_obs = project_dipole_on_cmb(b_vec_obs)

    # Consistency check (not strictly necessary, but nice for sanity)
    # amp_obs and amp_obs2 should match to numerical precision
    amp_obs = amp_obs2

    return dict(
                bcut         = float(bcut),
                tag          = tag,
                label        = label,
                pix          = pix,
                N_good       = N_good,
                z_good       = z_good,
                f            = f,
                n_hat        = n_hat,
                obs_b_vec    = b_vec_obs,
                obs_amp      = amp_obs,
                obs_amp_par  = amp_par_obs,
                obs_frac_par = frac_par_obs,
                l_deg_obs    = l_deg_obs,
                b_deg_obs    = b_deg_obs,
               )


def run_rotation_mc_for_slice(slice_info, z_mid, n_mock, rng):
    """
    Given prepared slice_info and z_mid, run rotation MC and
    return a dict with MC arrays and p-values.
    """
    n_hat        = slice_info["n_hat"]     # (N_good, 3)
    f            = slice_info["f"]
    N_good       = slice_info["N_good"]

    obs_amp      = slice_info["obs_amp"]
    obs_amp_par  = slice_info["obs_amp_par"]
    obs_frac_par = slice_info["obs_frac_par"]

    amps         = np.empty(n_mock, dtype=float)
    amp_pars     = np.empty(n_mock, dtype=float)
    frac_pars    = np.empty(n_mock, dtype=float)

    for i in range(n_mock):
        # Random rotation
        R     = random_rotation_matrix(rng)      # (3,3)
        
        # Rotate directions: row-vectors n_hat @ R^T
        n_rot = n_hat @ R.T                      # (N_good, 3)

        # Fit dipole in rotated frame
        a0_mc, b_vec_mc, _               = fit_dipole_linear(f, n_rot, w=N_good)
        amp_mc, _, _                     = dipole_summary_from_b(b_vec_mc)
        amp_mc2, amp_par_mc, frac_par_mc = project_dipole_on_cmb(b_vec_mc)

        amps[i]      = amp_mc2
        amp_pars[i]  = amp_par_mc
        frac_pars[i] = frac_par_mc

    # p-values
    p_amp            = float(np.mean(amps >= obs_amp))
    p_par_abs        = float(np.mean(np.abs(amp_pars)  >= np.abs(obs_amp_par)))
    p_frac_abs       = float(np.mean(np.abs(frac_pars) >= np.abs(obs_frac_par)))

    return dict(
                amps         = amps,
                amp_pars     = amp_pars,
                frac_pars    = frac_pars,
                obs_amp      = obs_amp,
                obs_amp_par  = obs_amp_par,
                obs_frac_par = obs_frac_par,
                z_mid        = z_mid,
                p_amp        = p_amp,
                p_par_abs    = p_par_abs,
                p_frac_abs   = p_frac_abs,
               )


# ----------------------------------------------------------------------
# Main driver
# ----------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--n-mock",
                        type     = int,
                        required = True,
                        help     = "number of rotation mocks per (bcut, z-slice)",
                       )
    
    parser.add_argument("--seed",
                        type     = int,
                        default  = 123,
                        help     = "random seed for rotation MC",
                       )
    
    parser.add_argument("--bcut-list",
                        type     = float,
                        nargs    = "*",
                        default  = BCUT_LIST,
                        help     = "optional list of |b| cuts (default: 10 20 30)",
                       )
    
    args = parser.parse_args()
    rng  = np.random.default_rng(args.seed)

    bcut_list = args.bcut_list
    print(f"[mc-rot] starting rotation MC: n_mock={args.n_mock}, seed={args.seed}")
    print(f"[mc-rot] bcut_list = {bcut_list}")
    print(f"[mc-rot] OUT_DIR = {OUT_DIR}")

    # Build slice definitions: full + each Z_BINS
    slices = [("full", "full sample", np.nan)]
    for (z_lo, z_hi) in Z_BINS:
        tag   = f"z{z_lo:.2f}_{z_hi:.2f}".replace(".", "p")
        label = f"{z_lo:.2f} ≤ z < {z_hi:.2f}"
        z_mid = 0.5 * (z_lo + z_hi)
        slices.append((tag, label, z_mid))

    for bcut in bcut_list:
        for (tag, label, z_mid) in slices:

            # Prepare
            slice_info = prepare_slice(bcut, tag, label)
            if slice_info is None:
                continue

            # Build a string for z_mid that handles NaN cleanly
            if np.isnan(z_mid):
                z_mid_str = "nan"
            else:
                z_mid_str = f"{z_mid:.2f}"
            
            N_good_obs = slice_info["N_good"].size
            
            print(
                  f"[mc-rot] |b|>{bcut:.1f}, tag={tag}, "
                  f"label='{label}', z_mid={z_mid_str}, N_good={N_good_obs}"
                 )


            # Run MC
            mc = run_rotation_mc_for_slice(slice_info, z_mid, args.n_mock, rng)

            print(
                  f"[mc-rot]   obs A={mc['obs_amp']:.3e}, "
                  f"A_par={mc['obs_amp_par']:.3e}, "
                  f"f_par={mc['obs_frac_par']:+.3f}"
                 )
            print(
                  f"[mc-rot]   p(A)={mc['p_amp']:.3f}, "
                  f"p(|A_par|)={mc['p_par_abs']:.3f}, "
                  f"p(|f_par|)={mc['p_frac_abs']:.3f}"
                 )

            # Save NPZ
            out_path = OUT_DIR / f"quaia_mc_bcut_rotate_{tag}_bcut{int(bcut)}.npz"
            np.savez(
                     out_path,
                     amps         = mc["amps"],
                     amp_pars     = mc["amp_pars"],
                     frac_pars    = mc["frac_pars"],
                     obs_amp      = mc["obs_amp"],
                     obs_amp_par  = mc["obs_amp_par"],
                     obs_frac_par = mc["obs_frac_par"],
                     z_mid        = mc["z_mid"],
                     label        = slice_info["label"],
                     bcut         = slice_info["bcut"],
                     N_good_obs   = slice_info["N_good"].size,
                     p_amp        = mc["p_amp"],
                     p_par_abs    = mc["p_par_abs"],
                     p_frac_abs   = mc["p_frac_abs"],
                    )
            print(f"[mc-rot]   wrote -> {out_path}")

    print("[mc-rot] done.")


if __name__ == "__main__":
    main()
