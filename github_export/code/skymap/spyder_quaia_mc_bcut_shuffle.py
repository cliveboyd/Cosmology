#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 18:26:38 2025

@author: boyde
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quaia_mc_bcut_shuffle.py  (0V1)

Monte Carlo test of the <z> dipole vs redshift, using pixel-wise shuffling of
z_mean within each map as the isotropic null.

For each Galactic latitude cut |b| > b_cut and each redshift slice, we:
  * load the pre-built <z> maps (quaia_build_zmean_maps_bcut.py),
  * fit the observed dipole,
  * generate N_mock shuffled realisations,
  * refit the dipole for each realisation,
  * build null distributions for A, A_par and f_par = A_par / A,
  * estimate simple p-values for the observed values.

Usage (from Spyder/IPython):

    %run /Users/boyde/.spyder-py3/quaia_mc_bcut_shuffle.py --n-mock 1000 --seed 123

Outputs in OUT_DIR:
  - quaia_mc_bcut_shuffle_{tag}_bcut{BCUT}.npz  (raw MC distributions)
  - quaia_mc_bcut_shuffle_summary.txt           (per-slice summary & p-values)


RESULTS...
N_MIN_PER_PIX = 3

%run /Users/boyde/.spyder-py3/quaia_mc_bcut_shuffle.py --n-mock 2000 --seed 42
[mc] Quaia <z> dipole Monte Carlo (pixel shuffling)
[mc] N_mock per slice = 2000, seed = 42
[mc] OUT_DIR = /Users/boyde/.spyder-py3/skymap2/quaia_outputs
[mc] |b|>10.0, 0.10 ≤ z < 0.30, N_good=   218, A_obs=4.873e-03, A_par_obs=4.308e-03, f_par_obs=+0.884, p(A)=0.562, p(|A_par|)=0.213
[mc] |b|>10.0, 0.30 ≤ z < 0.50, N_good=  3307, A_obs=2.096e-03, A_par_obs=-1.676e-03, f_par_obs=-0.799, p(A)=0.277, p(|A_par|)=0.097
[mc] |b|>10.0, 0.50 ≤ z < 0.75, N_good=  7916, A_obs=1.827e-03, A_par_obs=-9.835e-04, f_par_obs=-0.538, p(A)=0.231, p(|A_par|)=0.244
[mc] |b|>10.0, 0.75 ≤ z < 1.00, N_good= 13075, A_obs=9.801e-04, A_par_obs=-5.151e-04, f_par_obs=-0.526, p(A)=0.486, p(|A_par|)=0.393
[mc] |b|>10.0, 1.00 ≤ z < 1.50, N_good= 29498, A_obs=1.633e-03, A_par_obs=-6.599e-04, f_par_obs=-0.404, p(A)=0.180, p(|A_par|)=0.359
[mc] |b|>10.0, 1.50 ≤ z < 2.50, N_good= 32791, A_obs=1.549e-03, A_par_obs=-1.294e-03, f_par_obs=-0.836, p(A)=0.637, p(|A_par|)=0.251
[mc] |b|>15.0, 0.10 ≤ z < 0.30, N_good=   211, A_obs=5.645e-03, A_par_obs=4.449e-03, f_par_obs=+0.788, p(A)=0.436, p(|A_par|)=0.190
[mc] |b|>15.0, 0.30 ≤ z < 0.50, N_good=  3215, A_obs=1.843e-03, A_par_obs=-1.633e-03, f_par_obs=-0.886, p(A)=0.381, p(|A_par|)=0.106
[mc] |b|>15.0, 0.50 ≤ z < 0.75, N_good=  7663, A_obs=1.830e-03, A_par_obs=-9.627e-04, f_par_obs=-0.526, p(A)=0.247, p(|A_par|)=0.257
[mc] |b|>15.0, 0.75 ≤ z < 1.00, N_good= 12628, A_obs=1.189e-03, A_par_obs=-3.977e-04, f_par_obs=-0.335, p(A)=0.368, p(|A_par|)=0.537
[mc] |b|>15.0, 1.00 ≤ z < 1.50, N_good= 28108, A_obs=1.947e-03, A_par_obs=-5.775e-04, f_par_obs=-0.297, p(A)=0.098, p(|A_par|)=0.439
[mc] |b|>15.0, 1.50 ≤ z < 2.50, N_good= 31080, A_obs=1.771e-03, A_par_obs=-1.524e-03, f_par_obs=-0.860, p(A)=0.551, p(|A_par|)=0.221
[mc] |b|>20.0, 0.10 ≤ z < 0.30, N_good=   199, A_obs=8.057e-03, A_par_obs=6.556e-03, f_par_obs=+0.814, p(A)=0.168, p(|A_par|)=0.053
[mc] |b|>20.0, 0.30 ≤ z < 0.50, N_good=  3011, A_obs=1.681e-03, A_par_obs=-1.453e-03, f_par_obs=-0.865, p(A)=0.495, p(|A_par|)=0.174
[mc] |b|>20.0, 0.50 ≤ z < 0.75, N_good=  7179, A_obs=1.416e-03, A_par_obs=-6.466e-04, f_par_obs=-0.457, p(A)=0.517, p(|A_par|)=0.469
[mc] |b|>20.0, 0.75 ≤ z < 1.00, N_good= 11833, A_obs=1.524e-03, A_par_obs=-3.438e-04, f_par_obs=-0.226, p(A)=0.209, p(|A_par|)=0.621
[mc] |b|>20.0, 1.00 ≤ z < 1.50, N_good= 25934, A_obs=1.587e-03, A_par_obs=-4.768e-04, f_par_obs=-0.300, p(A)=0.276, p(|A_par|)=0.544
[mc] |b|>20.0, 1.50 ≤ z < 2.50, N_good= 28456, A_obs=2.005e-03, A_par_obs=-1.758e-03, f_par_obs=-0.877, p(A)=0.470, p(|A_par|)=0.161
[mc] |b|>25.0, 0.10 ≤ z < 0.30, N_good=   180, A_obs=6.449e-03, A_par_obs=5.271e-03, f_par_obs=+0.817, p(A)=0.393, p(|A_par|)=0.161
[mc] |b|>25.0, 0.30 ≤ z < 0.50, N_good=  2758, A_obs=1.998e-03, A_par_obs=-1.814e-03, f_par_obs=-0.908, p(A)=0.434, p(|A_par|)=0.115
[mc] |b|>25.0, 0.50 ≤ z < 0.75, N_good=  6558, A_obs=1.179e-03, A_par_obs=-5.058e-04, f_par_obs=-0.429, p(A)=0.694, p(|A_par|)=0.601
[mc] |b|>25.0, 0.75 ≤ z < 1.00, N_good= 10779, A_obs=1.304e-03, A_par_obs=-2.187e-05, f_par_obs=-0.017, p(A)=0.407, p(|A_par|)=0.981
[mc] |b|>25.0, 1.00 ≤ z < 1.50, N_good= 23266, A_obs=2.482e-03, A_par_obs=-8.125e-04, f_par_obs=-0.327, p(A)=0.055, p(|A_par|)=0.344
[mc] |b|>25.0, 1.50 ≤ z < 2.50, N_good= 25388, A_obs=2.499e-03, A_par_obs=-2.058e-03, f_par_obs=-0.823, p(A)=0.363, p(|A_par|)=0.132
[mc] |b|>30.0, 0.10 ≤ z < 0.30, N_good=   157, A_obs=5.475e-03, A_par_obs=3.910e-03, f_par_obs=+0.714, p(A)=0.593, p(|A_par|)=0.343
[mc] |b|>30.0, 0.30 ≤ z < 0.50, N_good=  2457, A_obs=2.171e-03, A_par_obs=-2.041e-03, f_par_obs=-0.940, p(A)=0.445, p(|A_par|)=0.131
[mc] |b|>30.0, 0.50 ≤ z < 0.75, N_good=  5913, A_obs=1.096e-03, A_par_obs=-4.752e-04, f_par_obs=-0.433, p(A)=0.762, p(|A_par|)=0.648
[mc] |b|>30.0, 0.75 ≤ z < 1.00, N_good=  9664, A_obs=1.278e-03, A_par_obs=-2.486e-04, f_par_obs=-0.194, p(A)=0.474, p(|A_par|)=0.783
[mc] |b|>30.0, 1.00 ≤ z < 1.50, N_good= 20546, A_obs=2.647e-03, A_par_obs=-4.280e-04, f_par_obs=-0.162, p(A)=0.059, p(|A_par|)=0.625
[mc] |b|>30.0, 1.50 ≤ z < 2.50, N_good= 22326, A_obs=1.816e-03, A_par_obs=-1.684e-03, f_par_obs=-0.927, p(A)=0.674, p(|A_par|)=0.260
[mc] |b|>35.0, 0.10 ≤ z < 0.30, N_good=   128, A_obs=5.397e-03, A_par_obs=4.975e-03, f_par_obs=+0.922, p(A)=0.718, p(|A_par|)=0.288
[mc] |b|>35.0, 0.30 ≤ z < 0.50, N_good=  2165, A_obs=3.127e-03, A_par_obs=-2.438e-03, f_par_obs=-0.780, p(A)=0.218, p(|A_par|)=0.089
[mc] |b|>35.0, 0.50 ≤ z < 0.75, N_good=  5252, A_obs=1.427e-03, A_par_obs=-8.762e-04, f_par_obs=-0.614, p(A)=0.662, p(|A_par|)=0.461
[mc] |b|>35.0, 0.75 ≤ z < 1.00, N_good=  8509, A_obs=1.765e-03, A_par_obs=1.115e-05, f_par_obs=+0.006, p(A)=0.296, p(|A_par|)=0.989
[mc] |b|>35.0, 1.00 ≤ z < 1.50, N_good= 17867, A_obs=3.326e-03, A_par_obs=-1.834e-04, f_par_obs=-0.055, p(A)=0.025, p(|A_par|)=0.866
[mc] |b|>35.0, 1.50 ≤ z < 2.50, N_good= 19303, A_obs=2.236e-03, A_par_obs=-1.684e-03, f_par_obs=-0.753, p(A)=0.587, p(|A_par|)=0.307
[mc] wrote summary -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_shuffle_summary.txt
<Figure size 864x576 with 0 Axes>

%runfile /Users/boyde/.spyder-py3/quaia_plot_bcut_grid.py --wdir
Reloaded modules: quaia_config, quaia_dipole_fit
[plot-grid] reading /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_dipole_bcut_grid.txt
# bcut  z_mid   amp        amp_par    frac_par
10.0    0.20  4.873e-03  4 .308e-03   +0.884
10.0    0.40  2.096e-03  -1.675e-03   -0.799
10.0    0.62  1.827e-03  -9.835e-04   -0.538
10.0    0.88  9.801e-04  -5.152e-04   -0.526
10.0    1.25  1.633e-03  -6.600e-04   -0.404
10.0    2.00  1.549e-03  -1.295e-03   -0.836
15.0    0.20  5.645e-03   4.449e-03   +0.788
15.0    0.40  1.843e-03  -1.632e-03   -0.886
15.0    0.62  1.830e-03  -9.630e-04   -0.526
15.0    0.88  1.189e-03  -3.977e-04   -0.334
15.0    1.25  1.947e-03  -5.777e-04   -0.297
15.0    2.00  1.771e-03  -1.524e-03   -0.860
20.0    0.20  8.057e-03   6.555e-03   +0.814
20.0    0.40  1.681e-03  -1.453e-03   -0.865
20.0    0.62  1.416e-03  -6.466e-04   -0.457
20.0    0.88  1.524e-03  -3.439e-04   -0.226
20.0    1.25  1.587e-03  -4.770e-04   -0.301
20.0    2.00  2.005e-03  -1.758e-03   -0.877
25.0    0.20  6.449e-03   5.271e-03   +0.817
25.0    0.40  1.998e-03  -1.814e-03   -0.908
25.0    0.62  1.179e-03  -5.057e-04   -0.429
25.0    0.88  1.304e-03  -2.185e-05   -0.017
25.0    1.25  2.482e-03  -8.126e-04   -0.327
25.0    2.00  2.499e-03  -2.058e-03   -0.823
30.0    0.20  5.475e-03   3.910e-03   +0.714
30.0    0.40  2.171e-03  -2.041e-03   -0.940
30.0    0.62  1.096e-03  -4.751e-04   -0.434
30.0    0.88  1.278e-03  -2.485e-04   -0.194
30.0    1.25  2.647e-03  -4.278e-04   -0.162
30.0    2.00  1.816e-03  -1.684e-03   -0.927
35.0    0.20  5.397e-03   4.975e-03   +0.922
35.0    0.40  3.127e-03  -2.438e-03   -0.780
35.0    0.62  1.427e-03  -8.760e-04   -0.614
35.0    0.88  1.765e-03   1.109e-05   +0.006
35.0    1.25  3.326e-03  -1.833e-04   -0.055
35.0    2.00  2.236e-03  -1.684e-03   -0.753

RESULTS....
N_MIN_PER_PIX = 5
%run /Users/boyde/.spyder-py3/quaia_mc_bcut_shuffle.py --n-mock 2000 --seed 52
[mc] Quaia <z> dipole Monte Carlo (pixel shuffling)
[mc] N_mock per slice = 2000, seed = 52
[mc] OUT_DIR = /Users/boyde/.spyder-py3/skymap2/quaia_outputs
[mc] |b|>10.0, 0.10 ≤ z < 0.30, N_good=     3, A_obs=5.178e-02, A_par_obs=-9.832e-03, f_par_obs=-0.190, p(A)=0.332, p(|A_par|)=0.494
[mc] |b|>10.0, 0.30 ≤ z < 0.50, N_good=   195, A_obs=4.767e-03, A_par_obs=-1.065e-03, f_par_obs=-0.224, p(A)=0.641, p(|A_par|)=0.757
[mc] |b|>10.0, 0.50 ≤ z < 0.75, N_good=  1069, A_obs=3.817e-03, A_par_obs=1.514e-03, f_par_obs=+0.397, p(A)=0.261, p(|A_par|)=0.406
[mc] |b|>10.0, 0.75 ≤ z < 1.00, N_good=  2868, A_obs=1.169e-03, A_par_obs=5.315e-04, f_par_obs=+0.455, p(A)=0.775, p(|A_par|)=0.621
[mc] |b|>10.0, 1.00 ≤ z < 1.50, N_good= 19703, A_obs=1.971e-03, A_par_obs=-1.150e-03, f_par_obs=-0.584, p(A)=0.116, p(|A_par|)=0.137
[mc] |b|>10.0, 1.50 ≤ z < 2.50, N_good= 27815, A_obs=1.667e-03, A_par_obs=-1.496e-03, f_par_obs=-0.897, p(A)=0.563, p(|A_par|)=0.203
[mc] |b|>15.0, 0.10 ≤ z < 0.30, N_good=     3, A_obs=5.178e-02, A_par_obs=-9.832e-03, f_par_obs=-0.190, p(A)=0.339, p(|A_par|)=0.503
[mc] |b|>15.0, 0.30 ≤ z < 0.50, N_good=   190, A_obs=3.689e-03, A_par_obs=-9.038e-04, f_par_obs=-0.245, p(A)=0.804, p(|A_par|)=0.807
[mc] |b|>15.0, 0.50 ≤ z < 0.75, N_good=  1054, A_obs=3.782e-03, A_par_obs=1.188e-03, f_par_obs=+0.314, p(A)=0.274, p(|A_par|)=0.502
[mc] |b|>15.0, 0.75 ≤ z < 1.00, N_good=  2811, A_obs=1.208e-03, A_par_obs=5.827e-04, f_par_obs=+0.483, p(A)=0.768, p(|A_par|)=0.614
[mc] |b|>15.0, 1.00 ≤ z < 1.50, N_good= 19193, A_obs=2.346e-03, A_par_obs=-1.057e-03, f_par_obs=-0.450, p(A)=0.048, p(|A_par|)=0.168
[mc] |b|>15.0, 1.50 ≤ z < 2.50, N_good= 26918, A_obs=1.832e-03, A_par_obs=-1.743e-03, f_par_obs=-0.951, p(A)=0.506, p(|A_par|)=0.131
[mc] tag=z0p10_0p30, |b|>20.0: too few good pixels (2)
[mc] |b|>20.0, 0.30 ≤ z < 0.50, N_good=   181, A_obs=3.368e-03, A_par_obs=-9.428e-07, f_par_obs=-0.000, p(A)=0.865, p(|A_par|)=1.000
[mc] |b|>20.0, 0.50 ≤ z < 0.75, N_good=  1010, A_obs=4.502e-03, A_par_obs=6.625e-04, f_par_obs=+0.147, p(A)=0.170, p(|A_par|)=0.720
[mc] |b|>20.0, 0.75 ≤ z < 1.00, N_good=  2692, A_obs=1.099e-03, A_par_obs=5.279e-04, f_par_obs=+0.481, p(A)=0.827, p(|A_par|)=0.658
[mc] |b|>20.0, 1.00 ≤ z < 1.50, N_good= 18090, A_obs=1.976e-03, A_par_obs=-8.562e-04, f_par_obs=-0.433, p(A)=0.137, p(|A_par|)=0.300
[mc] |b|>20.0, 1.50 ≤ z < 2.50, N_good= 25107, A_obs=2.135e-03, A_par_obs=-1.962e-03, f_par_obs=-0.919, p(A)=0.406, p(|A_par|)=0.115
[mc] tag=z0p10_0p30, |b|>25.0: too few good pixels (2)
[mc] |b|>25.0, 0.30 ≤ z < 0.50, N_good=   174, A_obs=5.065e-03, A_par_obs=-1.390e-03, f_par_obs=-0.274, p(A)=0.660, p(|A_par|)=0.724
[mc] |b|>25.0, 0.50 ≤ z < 0.75, N_good=   953, A_obs=4.772e-03, A_par_obs=1.515e-03, f_par_obs=+0.318, p(A)=0.174, p(|A_par|)=0.456
[mc] |b|>25.0, 0.75 ≤ z < 1.00, N_good=  2482, A_obs=1.443e-03, A_par_obs=3.728e-04, f_par_obs=+0.258, p(A)=0.698, p(|A_par|)=0.766
[mc] |b|>25.0, 1.00 ≤ z < 1.50, N_good= 16538, A_obs=2.542e-03, A_par_obs=-8.229e-04, f_par_obs=-0.324, p(A)=0.053, p(|A_par|)=0.349
[mc] |b|>25.0, 1.50 ≤ z < 2.50, N_good= 22711, A_obs=2.509e-03, A_par_obs=-2.121e-03, f_par_obs=-0.845, p(A)=0.341, p(|A_par|)=0.107
[mc] tag=z0p10_0p30, |b|>30.0: too few good pixels (1)
[mc] |b|>30.0, 0.30 ≤ z < 0.50, N_good=   149, A_obs=1.697e-03, A_par_obs=3.077e-04, f_par_obs=+0.181, p(A)=0.983, p(|A_par|)=0.942
[mc] |b|>30.0, 0.50 ≤ z < 0.75, N_good=   877, A_obs=4.123e-03, A_par_obs=1.478e-03, f_par_obs=+0.358, p(A)=0.341, p(|A_par|)=0.494
[mc] |b|>30.0, 0.75 ≤ z < 1.00, N_good=  2263, A_obs=1.130e-03, A_par_obs=1.018e-03, f_par_obs=+0.901, p(A)=0.866, p(|A_par|)=0.446
[mc] |b|>30.0, 1.00 ≤ z < 1.50, N_good= 14835, A_obs=2.603e-03, A_par_obs=-6.149e-04, f_par_obs=-0.236, p(A)=0.083, p(|A_par|)=0.516
[mc] |b|>30.0, 1.50 ≤ z < 2.50, N_good= 20166, A_obs=1.448e-03, A_par_obs=-1.408e-03, f_par_obs=-0.973, p(A)=0.779, p(|A_par|)=0.316
[mc] tag=z0p10_0p30, |b|>35.0: too few good pixels (1)
[mc] |b|>35.0, 0.30 ≤ z < 0.50, N_good=   134, A_obs=2.055e-03, A_par_obs=1.003e-03, f_par_obs=+0.488, p(A)=0.980, p(|A_par|)=0.856
[mc] |b|>35.0, 0.50 ≤ z < 0.75, N_good=   805, A_obs=3.923e-03, A_par_obs=1.087e-03, f_par_obs=+0.277, p(A)=0.457, p(|A_par|)=0.682
[mc] |b|>35.0, 0.75 ≤ z < 1.00, N_good=  2032, A_obs=1.075e-03, A_par_obs=1.006e-03, f_par_obs=+0.936, p(A)=0.894, p(|A_par|)=0.495
[mc] |b|>35.0, 1.00 ≤ z < 1.50, N_good= 13014, A_obs=3.405e-03, A_par_obs=-1.483e-04, f_par_obs=-0.044, p(A)=0.030, p(|A_par|)=0.889
[mc] |b|>35.0, 1.50 ≤ z < 2.50, N_good= 17588, A_obs=1.690e-03, A_par_obs=-1.420e-03, f_par_obs=-0.840, p(A)=0.765, p(|A_par|)=0.382
[mc] wrote summary -> /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_mc_bcut_shuffle_summary.txt

%runfile /Users/boyde/.spyder-py3/quaia_plot_bcut_grid.py --wdir
Reloaded modules: quaia_config
[plot-grid] reading /Users/boyde/.spyder-py3/skymap2/quaia_outputs/quaia_dipole_bcut_grid.txt
# bcut  z_mid   amp        amp_par    frac_par
10.0   0.20  5.178e-02  -9.836e-03   -0.190
10.0   0.40  4.767e-03  -1.066e-03   -0.224
10.0   0.62  3.817e-03  1.514e-03   +0.397
10.0   0.88  1.169e-03  5.314e-04   +0.455
10.0   1.25  1.971e-03  -1.150e-03   -0.584
10.0   2.00  1.667e-03  -1.496e-03   -0.897
15.0   0.20  5.178e-02  -9.836e-03   -0.190
15.0   0.40  3.689e-03  -9.037e-04   -0.245
15.0   0.62  3.782e-03  1.188e-03   +0.314
15.0   0.88  1.208e-03  5.829e-04   +0.483
15.0   1.25  2.346e-03  -1.056e-03   -0.450
15.0   2.00  1.832e-03  -1.743e-03   -0.951
20.0   0.40  3.368e-03  -1.176e-06   -0.000
20.0   0.62  4.502e-03  6.623e-04   +0.147
20.0   0.88  1.099e-03  5.281e-04   +0.481
20.0   1.25  1.976e-03  -8.563e-04   -0.433
20.0   2.00  2.135e-03  -1.962e-03   -0.919
25.0   0.40  5.065e-03  -1.390e-03   -0.274
25.0   0.62  4.772e-03  1.515e-03   +0.317
25.0   0.88  1.443e-03  3.727e-04   +0.258
25.0   1.25  2.542e-03  -8.230e-04   -0.324
25.0   2.00  2.509e-03  -2.121e-03   -0.845
30.0   0.40  1.697e-03  3.078e-04   +0.181
30.0   0.62  4.123e-03  1.478e-03   +0.358
30.0   0.88  1.130e-03  1.018e-03   +0.901
30.0   1.25  2.603e-03  -6.152e-04   -0.236
30.0   2.00  1.448e-03  -1.409e-03   -0.973
35.0   0.40  2.055e-03  1.003e-03   +0.488
35.0   0.62  3.923e-03  1.087e-03   +0.277
35.0   0.88  1.075e-03  1.006e-03   +0.936
35.0   1.25  3.405e-03  -1.485e-04   -0.044
35.0   2.00  1.690e-03  -1.420e-03   -0.840
"""

import argparse
from   pathlib          import Path

import numpy            as     np
import healpy           as     hp

from   quaia_config     import (OUT_DIR,
                                NSIDE,
                                N_MIN_PER_PIX,
                                Z_SLICES
                               ) 

from   quaia_dipole_fit import fit_dipole_linear, dipole_summary_from_b

# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------

BCUT_LIST = [10.0, 15.0, 20.0, 25.0, 30.0, 35.0]

# Use N_MIN_PER_PIX from quaia_config (do NOT override it here)

# CMB dipole direction (Galactic)
L_CMB_DEG = 264.0
B_CMB_DEG =  48.0

# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

def unit_vec_gal(l_deg, b_deg):
    """Unit vector in Galactic coords (l, b in deg)."""
    l    = np.deg2rad(l_deg)
    b    = np.deg2rad(b_deg)
    x    = np.cos(b) * np.cos(l)
    y    = np.cos(b) * np.sin(l)
    z    = np.sin(b)
    return np.array([x, y, z])


def project_dipole_on_cmb(b_vec):
    """
    Given dipole vector b_vec (shape (3,)), return:
      amp      = |b_vec|
      amp_par  = component of b_vec along CMB direction (signed)
      frac_par = amp_par / amp  (0 if amp == 0)
    """
    b_vec = np.asarray(b_vec, dtype=float)
    amp   = np.linalg.norm(b_vec)
    if amp == 0.0:
        return 0.0, 0.0, 0.0

    n_cmb    = unit_vec_gal(L_CMB_DEG, B_CMB_DEG)
    amp_par  = float(np.dot(b_vec, n_cmb))
    frac_par = amp_par / amp
    return amp, amp_par, frac_par


def load_map(tag, bcut):
    path = OUT_DIR / f"quaia_zmean_{tag}_bcut{int(bcut)}_nside{NSIDE}.npz"
    data = np.load(path, allow_pickle=False)
    return data["N"], data["zmean"]


def mc_for_slice(tag, label, z_mid, bcut, n_mock, rng):
    """
    Run MC for one (bcut, tag) combination.

    Returns a dict with observed values, MC stats and p-values,
    plus writes an NPZ file with raw distributions.
    """
    N_map, zmean = load_map(tag, bcut)

    good = (N_map >= N_MIN_PER_PIX) & np.isfinite(zmean)
    pix  = np.where(good)[0]

    if pix.size < 3:
        print(f"[mc] tag={tag}, |b|>{bcut:.1f}: too few good pixels ({pix.size})")
        return None

    N_good                         = N_map[good].astype(float)
    z_good                         = zmean[good].astype(float)

    # HEALPix geometry for good pixels
    theta, phi                     = hp.pix2ang(NSIDE, pix)
    x                              = np.sin(theta) * np.cos(phi)
    y                              = np.sin(theta) * np.sin(phi)
    z                              = np.cos(theta)
    n_hat                          = np.column_stack([x, y, z])

    # Observed dipole
    f_obs                          = z_good - z_good.mean()
    a0_obs, b_vec_obs, _           = fit_dipole_linear(f_obs, n_hat, w=N_good)
    amp_obs, amp_par_obs, frac_obs = project_dipole_on_cmb(b_vec_obs)

    # Monte Carlo null
    amps      = np.empty(n_mock, dtype=float)
    amp_pars  = np.empty(n_mock, dtype=float)
    frac_pars = np.empty(n_mock, dtype=float)

    for i in range(n_mock):
        z_shuff                     = rng.permutation(z_good)
        f_mc                        = z_shuff - z_shuff.mean()
        a0_mc, b_vec_mc, _          = fit_dipole_linear(f_mc, n_hat, w=N_good)
        amp_mc, amp_par_mc, frac_mc = project_dipole_on_cmb(b_vec_mc)

        amps[i]                     = amp_mc
        amp_pars[i]                 = amp_par_mc
        frac_pars[i]                = frac_mc

    # One-sided p-values
    p_amp            = float(np.mean(amps >= amp_obs))
    p_par_abs        = float(np.mean(np.abs(amp_pars)  >= np.abs(amp_par_obs)))
    p_frac_abs       = float(np.mean(np.abs(frac_pars) >= np.abs(frac_obs)))

    # Save raw distributions
    out_npz = OUT_DIR / f"quaia_mc_bcut_shuffle_{tag}_bcut{int(bcut)}.npz"
    np.savez(
        out_npz,
        amps         = amps,
        amp_pars     = amp_pars,
        frac_pars    = frac_pars,
        obs_amp      = amp_obs,
        obs_amp_par  = amp_par_obs,
        obs_frac_par = frac_obs,
        z_mid        = z_mid,
        label        = label,
        bcut         = bcut,
        N_good       = int(pix.size),
        p_amp        = p_amp,
        p_par_abs    = p_par_abs,
        p_frac_abs   = p_frac_abs,
    )

    print(
        f"[mc] |b|>{bcut:4.1f}, {label:15s}, N_good={pix.size:6d}, "
        f"A_obs={amp_obs:.3e}, A_par_obs={amp_par_obs:.3e}, f_par_obs={frac_obs:+.3f}, "
        f"p(A)={p_amp:.3f}, p(|A_par|)={p_par_abs:.3f}"
    )

    return dict(
        bcut       = bcut,
        tag        = tag,
        label      = label,
        z_mid      = z_mid,
        N_good     = int(pix.size),
        A_obs      = amp_obs,
        A_par_obs  = amp_par_obs,
        f_par_obs  = frac_obs,
        p_amp      = p_amp,
        p_par_abs  = p_par_abs,
        p_frac_abs = p_frac_abs,
    )


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-mock", type=int, default=1000,
                        help="Number of Monte Carlo shuffles per slice")
    parser.add_argument("--seed",   type=int, default=12345,
                        help="Random seed")
    args   = parser.parse_args()

    rng    = np.random.default_rng(args.seed)

    print("[mc] Quaia <z> dipole Monte Carlo (pixel shuffling)")
    print(f"[mc] N_mock per slice = {args.n_mock}, seed = {args.seed}")
    print(f"[mc] OUT_DIR = {OUT_DIR}")

    results = []

    # Loop over b-cuts and z-slices, skip the 'full' slice
    for bcut in BCUT_LIST:
        for (z_lo, z_hi, label, tag) in Z_SLICES:
            if tag == "full":
                continue

            z_mid = 0.5 * (z_lo + z_hi)
            res   = mc_for_slice(tag, label, z_mid, bcut, args.n_mock, rng)
            if res is not None:
                results.append(res)

    # Write compact text summary
    out_txt = OUT_DIR / "quaia_mc_bcut_shuffle_summary.txt"
    with out_txt.open("w") as f:
        f.write("# Monte Carlo summary for <z> dipole vs z (pixel shuffling)\n")
        f.write("# bcut  tag            z_mid   N_good   "
                "A_obs      A_par_obs   f_par_obs   p_A    p_|A_par|  p_|f_par|\n")
        for r in results:
            f.write(
                f"{r['bcut']:4.1f}  {r['tag']:12s}  {r['z_mid']:5.2f}  "
                f"{r['N_good']:6d}  "
                f"{r['A_obs']:.3e}  {r['A_par_obs']:.3e}  {r['f_par_obs']:+.3f}  "
                f"{r['p_amp']:.3f}  {r['p_par_abs']:.3f}  {r['p_frac_abs']:.3f}\n"
            )

    print(f"[mc] wrote summary -> {out_txt}")


if __name__ == "__main__":
    main()
