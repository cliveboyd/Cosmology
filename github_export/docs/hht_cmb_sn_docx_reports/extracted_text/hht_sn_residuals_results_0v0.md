# HHT_SN_Residuals Results 0V0

## Source

- Local DOCX: `C:\Users\clive\Documents\Cosmology\HHT_SN_Residuals Results 0V0.docx`
- Extracted: `2026-07-14T21:53:50`
- DOCX created: `2025-09-26T14:33:00Z`
- DOCX modified: `2025-09-26T17:43:00Z`

## Extracted Text

HHT_SN_Residuals Results 0V0

%run '/Users/boyde/.spyder-py3/hht_bao_lock_check.py' \

'/Users/boyde/Downloads/bao_long.csv' \

--rd 147.09 --out 'plamb_runs/hht_bao_lock_check'

[ok] corr(ridge,k)=-0.918, corr(ridge,1/k)=0.832

[ok] files at plamb_runs/hht_bao_lock_check

<Figure size 864x576 with 0 Axes>

%run '/Users/boyde/.spyder-py3/hht_ridge_bao_lock.py' \

'/Users/boyde/Downloads/bao_long.csv' \

--rd 147.09 --out 'plamb_runs/hht_bao_lock'

[ok] corr = -0.918

[ok] files at plamb_runs/hht_bao_lock

<Figure size 864x576 with 0 Axes>

%run '/Users/boyde/.spyder-py3/hht_surr_significance.py' \

'/Users/boyde/.spyder-py3/PLamb_Test_10V6c_plus.py' \

'/Users/boyde/Downloads/Pantheon+SH0ES.dat' \

--nsurr 500 --method simpson --out 'plamb_runs/hht_surr'

usage: PLamb_Test_10V6c_plus.py [-h] [--out OUT] [--H0 H0] [--Om OM] [--Ol OL]

[--A_acc A_ACC] [--n_acc N_ACC]

[--gamma_c GAMMA_C] [--epsilon_M EPSILON_M]

[--zmax ZMAX] [--zmin ZMIN]

[--nwalkers NWALKERS] [--nsteps NSTEPS]

[--nburn NBURN] [--seed SEED]

[--grid-points GRID_POINTS] [--nint NINT]

[--thin THIN] [--fix_H0] [--fix_Om] [--fix_Ol]

[--fix_A_acc] [--fix_n_acc] [--fix_gamma_c]

[--fix_epsilon_M] [--show] [--verbose]

[--flat] [--save-chain] [--link-em-gc]

[--planck-c-variant]

[--model {LCDM,KdS,FR,PBH}]

[--dM-mode {linear_z,log1pz}]

[--sampler {emcee,grid,none}] [--data DATA]

[--cov COV] [--bao BAO] [--bao-cov BAO_COV]

[--cc CC] [--planck PLANCK] [--rd RD]

[--planck-zstar PLANCK_ZSTAR] [--rs RS]

[--planck-nint PLANCK_NINT]

[--kappa-em-gc KAPPA_EM_GC] [--tight-priors]

[--corner] [--diag-cov] [--preflight]

[--fluct-Lambda]

[--Lambda-model {white,OU,piecewise,spline}]

[--Lambda-amp LAMBDA_AMP]

[--Lambda-corrz LAMBDA_CORRZ]

[--Lambda-knots LAMBDA_KNOTS]

[--seed-Lambda SEED_LAMBDA]

[--couple-Lambda-Aacc]

[--lambda-aacc-gamma LAMBDA_AACC_GAMMA]

[--discrete-events] [--event-size EVENT_SIZE]

[--noether {none,conserve_time,conserve_scale,conserve_diffeo,break_time_soft,break_scale_soft,break_diffeo_soft,break_space_soft,conserve_time_soft,conserve_scale_soft,conserve_c_soft,conserve_lum_soft,conserve_space_soft,combo_soft}]

[--noether-strength NOETHER_STRENGTH]

[--bao-c-mode {c0,c_of_z}] [--fr-nonexpanding]

[--link-sigma LINK_SIGMA] [--link-soft]

[--c-model {linear_z,log1pz,piecewise,density_coupled,power1pz,expz,saturating,poly_log1pz,cheb_log1pz}]

[--gamma2_c GAMMA2_C] [--z_t Z_T]

[--alpha_rho ALPHA_RHO] [--c-power C_POWER]

[--c-z0 C_Z0] [--c-sat-kind {tanh,exp}]

[--c-coefs C_COEFS]

[--planck-c-mode {none,R_only,R_and_rs}]

[--rs-power RS_POWER] [--pbh-tests]

[--pbh-couple] [--pbh-M PBH_M]

[--pbh-a-dim PBH_A_DIM] [--pbh-G PBH_G]

[--pbh-c PBH_C] [--pbh-hbar PBH_HBAR]

[--pbh-kB PBH_KB] [--spec-N SPEC_N]

[--spec-beta SPEC_BETA]

[--spec-seed SPEC_SEED]

[--hht-slices HHT_SLICES] [--pbh-out PBH_OUT]

[--pbh-plots]

[--pbh-plot-slices PBH_PLOT_SLICES]

PLamb_Test_10V6c_plus.py: error: unrecognized arguments: --residuals-only

/Users/boyde/.spyder-py3/hht_surr_significance.py:63: FutureWarning: The 'delim_whitespace' keyword in pd.read_csv is deprecated and will be removed in a future version. Use ``sep='\s+'`` instead

df               = pd.read_csv(args.sn_table, delim_whitespace=True, comment='#')

/opt/anaconda3/lib/python3.12/site-packages/matplotlib/ticker.py:2119: RuntimeWarning: overflow encountered in multiply

steps = self._extended_steps * scale

[ok] CSV  -> plamb_runs/hht_surr/hht_imf_significance.csv

[ok] Plot -> plamb_runs/hht_surr/hht_imf_significance.png

<Figure size 864x576 with 0 Axes>

%run '/Users/boyde/.spyder-py3/hht_sn_residuals.py' \

'/Users/boyde/.spyder-py3/PLamb_Test_10V6c_plus.py' \

'/Users/boyde/Downloads/Pantheon+SH0ES.dat' \

'plamb_runs/hht_sn_noether/combo_soft6' \

--noether combo_soft --noether-strength 6 --zmin 0.02

[out] writing to: /Users/boyde/.spyder-py3/plamb_runs/hht_sn_noether/combo_soft6

[clip] kept 0.98 of points.

[grid] uniform L bins = 1111, L ∈ [0.0012,1.1821]

[save] IMF stats -> /Users/boyde/.spyder-py3/plamb_runs/hht_sn_noether/combo_soft6/hht_sn_imf_stats.csv, /Users/boyde/.spyder-py3/plamb_runs/hht_sn_noether/combo_soft6/hht_sn_imf_stats.npz

[files] /Users/boyde/.spyder-py3/plamb_runs/hht_sn_noether/combo_soft6/hht_sn.png

<Figure size 864x576 with 0 Axes>

%run '/Users/boyde/.spyder-py3/hht_sn_residuals.py' \

'/Users/boyde/.spyder-py3/PLamb_Test_10V6c_plus.py' \

'/Users/boyde/Downloads/Pantheon+SH0ES.dat' \

'plamb_runs/hht_sn_noether/conserve_time' \

--noether conserve_time --zmin 0.02

[out] writing to: /Users/boyde/.spyder-py3/plamb_runs/hht_sn_noether/conserve_time

[clip] kept 0.98 of points.

[grid] uniform L bins = 1111, L ∈ [0.0012,1.1821]

[save] IMF stats -> /Users/boyde/.spyder-py3/plamb_runs/hht_sn_noether/conserve_time/hht_sn_imf_stats.csv, /Users/boyde/.spyder-py3/plamb_runs/hht_sn_noether/conserve_time/hht_sn_imf_stats.npz

[files] /Users/boyde/.spyder-py3/plamb_runs/hht_sn_noether/conserve_time/hht_sn.png

<Figure size 864x576 with 0 Axes>

%run '/Users/boyde/.spyder-py3/hht_sn_residuals.py' \

'/Users/boyde/.spyder-py3/PLamb_Test_10V6c_plus.py' \

'/Users/boyde/Downloads/Pantheon+SH0ES.dat' \

'plamb_runs/hht_sn_zmin010' --zmin 0.10

[out] writing to: /Users/boyde/.spyder-py3/plamb_runs/hht_sn_zmin010

[clip] kept 0.98 of points.

[grid] uniform L bins = 1111, L ∈ [0.0012,1.1821]

[save] IMF stats -> /Users/boyde/.spyder-py3/plamb_runs/hht_sn_zmin010/hht_sn_imf_stats.csv, /Users/boyde/.spyder-py3/plamb_runs/hht_sn_zmin010/hht_sn_imf_stats.npz

[files] /Users/boyde/.spyder-py3/plamb_runs/hht_sn_zmin010/hht_sn.png

<Figure size 864x576 with 0 Axes>

%run '/Users/boyde/.spyder-py3/hht_sn_residuals.py' \

'/Users/boyde/.spyder-py3/PLamb_Test_10V6c_plus.py' \

'/Users/boyde/Downloads/Pantheon+SH0ES.dat' \

'plamb_runs/hht_sn_zmin010' --zmin 0.10

[out] writing to: /Users/boyde/.spyder-py3/plamb_runs/hht_sn_zmin010

[clip] kept 0.98 of points.

[grid] uniform L bins = 1111, L ∈ [0.0012,1.1821]

[save] IMF stats -> /Users/boyde/.spyder-py3/plamb_runs/hht_sn_zmin010/hht_sn_imf_stats.csv, /Users/boyde/.spyder-py3/plamb_runs/hht_sn_zmin010/hht_sn_imf_stats.npz

[files] /Users/boyde/.spyder-py3/plamb_runs/hht_sn_zmin010/hht_sn.png

<Figure size 864x576 with 0 Axes>

%run '/Users/boyde/.spyder-py3/hht_sn_residuals.py' \

'/Users/boyde/.spyder-py3/PLamb_Test_10V6c_plus.py' \

'/Users/boyde/Downloads/Pantheon+SH0ES.dat' \

'plamb_runs/hht_sn_zmin002' --zmin 0.02

[out] writing to: /Users/boyde/.spyder-py3/plamb_runs/hht_sn_zmin002

[clip] kept 0.98 of points.

[grid] uniform L bins = 1111, L ∈ [0.0012,1.1821]

[save] IMF stats -> /Users/boyde/.spyder-py3/plamb_runs/hht_sn_zmin002/hht_sn_imf_stats.csv, /Users/boyde/.spyder-py3/plamb_runs/hht_sn_zmin002/hht_sn_imf_stats.npz

[files] /Users/boyde/.spyder-py3/plamb_runs/hht_sn_zmin002/hht_sn.png

<Figure size 864x576 with 0 Axes>

%run '/Users/boyde/.spyder-py3/hht_integrator_ab.py' \

'/Users/boyde/.spyder-py3/PLamb_Test_10V6c_plus.py' \

'/Users/boyde/Downloads/Pantheon+SH0ES.dat'

[files] plamb_runs/hht_integrator/hht_integrator_compare.png

<Figure size 864x576 with 0 Axes>

%run '/Users/boyde/.spyder-py3/hht_sn_residuals.py' \

'/Users/boyde/.spyder-py3/PLamb_Test_10V6c_plus.py' \

'/Users/boyde/Downloads/Pantheon+SH0ES.dat' \

'plamb_runs/hht_couple/ou_seed777_ON'

[out] writing to: /Users/boyde/.spyder-py3/plamb_runs/hht_couple/ou_seed777_ON

[clip] kept 0.98 of points.

[grid] uniform L bins = 1111, L ∈ [0.0012,1.1821]

[save] IMF stats -> /Users/boyde/.spyder-py3/plamb_runs/hht_couple/ou_seed777_ON/hht_sn_imf_stats.csv, /Users/boyde/.spyder-py3/plamb_runs/hht_couple/ou_seed777_ON/hht_sn_imf_stats.npz

[files] /Users/boyde/.spyder-py3/plamb_runs/hht_couple/ou_seed777_ON/hht_sn.png

<Figure size 864x576 with 0 Axes>

%run '/Users/boyde/.spyder-py3/PLamb_Test_10V6c_plus.py' \

--model FR \

--data '/Users/boyde/Downloads/Pantheon+SH0ES.dat' \

--sampler none --discrete-events --event-size 0.005 \

--out 'plamb_runs/hht_compare/discrete_0005'

[status] parsing done.

[status] data loaded.

[status] active params: ['H0', 'Om', 'Ol', 'A_acc', 'n_acc', 'gamma_c', 'epsilon_M']  frozen: {}

[Integrator] discrete-events enabled  event-size Δz=0.005

[link] soft=False  kappa=1.0  sigma=0.02

Best-fit (model=FR, sampler=none): {'H0': 70.0, 'Om': 0.3, 'Ol': 0.7, 'Ok': 0.0, 'A_acc': 0.0, 'n_acc': 0.0, 'gamma_c': 0.0, 'epsilon_M': 0.0, 'lnpost': -615.6178439615401, 'chi2_SN': 1231.2356879230801, 'chi2_BAO': 0.0, 'chi2_CC': 0.0, 'chi2_Planck': 0.0, 'chi2': 1231.2356879230801, 'dof_total': 1694, 'dof_SN': 1694, 'chi2/dof_total': 0.7268215395059505, 'chi2_SN/dof_SN': 0.7268215395059505, 'AIC': 1245.2356879230801, 'BIC': 1283.3084890698512, 'N_SN': 1701, 'N_eff': 1701, 'active': ['H0', 'Om', 'Ol', 'A_acc', 'n_acc', 'gamma_c', 'epsilon_M'], 'frozen': {}, 'w_eff': -1.0, 'bao_c_mode': 'c_of_z', 'c_model': 'linear_z', 'planck_c_mode': 'none', 'rs_power': 1.0}

[fit]  chi2_total=1231.24  dof_total=1694  chi2/dof_total=0.727   AIC=1245.24  BIC=1283.31

[fit]  (SN-only) chi2_SN=1231.24  dof_SN=1694  chi2_SN/dof_SN=0.727

[timing] total 0.35 s

DONE 1) integrator fingerprint test (Simpson vs discrete-events)

# Simpson baseline

%run '/Users/boyde/.spyder-py3/PLamb_Test_10V6c_plus.py' \

--model FR --data '/Users/boyde/Downloads/Pantheon+SH0ES.dat' \

--sampler none --out 'plamb_runs/hht_compare/simpson'

[status] parsing done.

[status] data loaded.

[status] active params: ['H0', 'Om', 'Ol', 'A_acc', 'n_acc', 'gamma_c', 'epsilon_M']  frozen: {}

[link] soft=False  kappa=1.0  sigma=0.02

Best-fit (model=FR, sampler=none): {'H0': 70.0, 'Om': 0.3, 'Ol': 0.7, 'Ok': 0.0, 'A_acc': 0.0, 'n_acc': 0.0, 'gamma_c': 0.0, 'epsilon_M': 0.0, 'lnpost': -605.368998425507, 'chi2_SN': 1210.737996851014, 'chi2_BAO': 0.0, 'chi2_CC': 0.0, 'chi2_Planck': 0.0, 'chi2': 1210.737996851014, 'dof_total': 1694, 'dof_SN': 1694, 'chi2/dof_total': 0.7147213676806459, 'chi2_SN/dof_SN': 0.7147213676806459, 'AIC': 1224.737996851014, 'BIC': 1262.8107979977851, 'N_SN': 1701, 'N_eff': 1701, 'active': ['H0', 'Om', 'Ol', 'A_acc', 'n_acc', 'gamma_c', 'epsilon_M'], 'frozen': {}, 'w_eff': -1.0, 'bao_c_mode': 'c_of_z', 'c_model': 'linear_z', 'planck_c_mode': 'none', 'rs_power': 1.0}

[fit]  chi2_total=1210.74  dof_total=1694  chi2/dof_total=0.715   AIC=1224.74  BIC=1262.81

[fit]  (SN-only) chi2_SN=1210.74  dof_SN=1694  chi2_SN/dof_SN=0.715

[timing] total 0.04 s

Best-fit saved to plamb_runs/hht_compare/simpson/FR_none_bestfit.txt

DONE 1) integrator fingerprint test (Simpson vs discrete-events)

%run '/Users/boyde/.spyder-py3/hht_sn_residuals.py' \

'/Users/boyde/.spyder-py3/PLamb_Test_10V6c_plus.py' \

'/Users/boyde/Downloads/Pantheon+SH0ES.dat' \

'plamb_runs/hht_sn_FRPBH'

[out] writing to: /Users/boyde/.spyder-py3/plamb_runs/hht_sn_FRPBH

[clip] kept 0.98 of points.

[grid] uniform L bins = 1111, L ∈ [0.0012,1.1821]

[save] IMF stats -> /Users/boyde/.spyder-py3/plamb_runs/hht_sn_FRPBH/hht_sn_imf_stats.csv, /Users/boyde/.spyder-py3/plamb_runs/hht_sn_FRPBH/hht_sn_imf_stats.npz

[files] /Users/boyde/.spyder-py3/plamb_runs/hht_sn_FRPBH/hht_sn.png

<Figure size 864x576 with 0 Axes>

[status] parsing done.

[status] data loaded.

[status] active params: ['H0', 'Om', 'Ol', 'A_acc', 'n_acc', 'gamma_c', 'epsilon_M']  frozen: {}

[link] soft=False  kappa=1.0  sigma=0.02

Best-fit (model=FR, sampler=none): {'H0': 70.0, 'Om': 0.3, 'Ol': 0.7, 'Ok': 0.0, 'A_acc': 0.0, 'n_acc': 0.0, 'gamma_c': 0.0, 'epsilon_M': 0.0, 'lnpost': -605.368998425507, 'chi2_SN': 1210.737996851014, 'chi2_BAO': 0.0, 'chi2_CC': 0.0, 'chi2_Planck': 0.0, 'chi2': 1210.737996851014, 'dof_total': 1694, 'dof_SN': 1694, 'chi2/dof_total': 0.7147213676806459, 'chi2_SN/dof_SN': 0.7147213676806459, 'AIC': 1224.737996851014, 'BIC': 1262.8107979977851, 'N_SN': 1701, 'N_eff': 1701, 'active': ['H0', 'Om', 'Ol', 'A_acc', 'n_acc', 'gamma_c', 'epsilon_M'], 'frozen': {}, 'w_eff': -1.0, 'bao_c_mode': 'c_of_z', 'c_model': 'linear_z', 'planck_c_mode': 'none', 'rs_power': 1.0}

[fit]  chi2_total=1210.74  dof_total=1694  chi2/dof_total=0.715   AIC=1224.74  BIC=1262.81

[fit]  (SN-only) chi2_SN=1210.74  dof_SN=1694  chi2_SN/dof_SN=0.715

[timing] total 0.04 s

Best-fit saved to plamb_runs/hht_compare/simpson/FR_none_bestfit.txt

DONE # discrete-events (Δz = 0.01)

%run '/Users/boyde/.spyder-py3/PLamb_Test_10V6c_plus.py' \

--model FR --data '/Users/boyde/Downloads/Pantheon+SH0ES.dat' \

--sampler none --discrete-events --event-size 0.01 \

--out 'plamb_runs/hht_compare/discrete_001'

[status] parsing done.

[status] data loaded.

[status] active params: ['H0', 'Om', 'Ol', 'A_acc', 'n_acc', 'gamma_c', 'epsilon_M']  frozen: {}

[Integrator] discrete-events enabled  event-size Δz=0.01

[link] soft=False  kappa=1.0  sigma=0.02

Best-fit (model=FR, sampler=none): {'H0': 70.0, 'Om': 0.3, 'Ol': 0.7, 'Ok': 0.0, 'A_acc': 0.0, 'n_acc': 0.0, 'gamma_c': 0.0, 'epsilon_M': 0.0, 'lnpost': -625.5295460871644, 'chi2_SN': 1251.0590921743287, 'chi2_BAO': 0.0, 'chi2_CC': 0.0, 'chi2_Planck': 0.0, 'chi2': 1251.0590921743287, 'dof_total': 1694, 'dof_SN': 1694, 'chi2/dof_total': 0.738523667163122, 'chi2_SN/dof_SN': 0.738523667163122, 'AIC': 1265.0590921743287, 'BIC': 1303.1318933210998, 'N_SN': 1701, 'N_eff': 1701, 'active': ['H0', 'Om', 'Ol', 'A_acc', 'n_acc', 'gamma_c', 'epsilon_M'], 'frozen': {}, 'w_eff': -1.0, 'bao_c_mode': 'c_of_z', 'c_model': 'linear_z', 'planck_c_mode': 'none', 'rs_power': 1.0}

[fit]  chi2_total=1251.06  dof_total=1694  chi2/dof_total=0.739   AIC=1265.06  BIC=1303.13

[fit]  (SN-only) chi2_SN=1251.06  dof_SN=1694  chi2_SN/dof_SN=0.739

[timing] total 0.18 s

Best-fit saved to plamb_runs/hht_compare/discrete_001/FR_none_bestfit.txt

DONE     # realize Λ(z) once (OU) and keep params fixed (sampler none)

%run '/Users/boyde/.spyder-py3/PLamb_Test_10V6c_plus.py' \

--model FR --data '/Users/boyde/Downloads/Pantheon+SH0ES.dat' \

--sampler none --fluct-Lambda --Lambda-model OU --Lambda-amp 1e-120 \

--Lambda-corrz 0.3 --seed-Lambda 777 \

--out 'plamb_runs/hht_couple/ou_seed777'

[status] parsing done.

[status] data loaded.

[status] active params: ['H0', 'Om', 'Ol', 'A_acc', 'n_acc', 'gamma_c', 'epsilon_M']  frozen: {}

[Lambda] model=OU  amp=1.000e-120  corrz=0.3  knots=6  seed=777  coupling=False  gamma=0.0

[link] soft=False  kappa=1.0  sigma=0.02

Best-fit (model=FR, sampler=none): {'H0': 70.0, 'Om': 0.3, 'Ol': 0.7, 'Ok': 0.0, 'A_acc': 0.0, 'n_acc': 0.0, 'gamma_c': 0.0, 'epsilon_M': 0.0, 'lnpost': -605.368998425507, 'chi2_SN': 1210.737996851014, 'chi2_BAO': 0.0, 'chi2_CC': 0.0, 'chi2_Planck': 0.0, 'chi2': 1210.737996851014, 'dof_total': 1694, 'dof_SN': 1694, 'chi2/dof_total': 0.7147213676806459, 'chi2_SN/dof_SN': 0.7147213676806459, 'AIC': 1224.737996851014, 'BIC': 1262.8107979977851, 'N_SN': 1701, 'N_eff': 1701, 'active': ['H0', 'Om', 'Ol', 'A_acc', 'n_acc', 'gamma_c', 'epsilon_M'], 'frozen': {}, 'w_eff': -1.0, 'bao_c_mode': 'c_of_z', 'c_model': 'linear_z', 'planck_c_mode': 'none', 'rs_power': 1.0}

[fit]  chi2_total=1210.74  dof_total=1694  chi2/dof_total=0.715   AIC=1224.74  BIC=1262.81

[fit]  (SN-only) chi2_SN=1210.74  dof_SN=1694  chi2_SN/dof_SN=0.715

[timing] total 0.04 s

Best-fit saved to plamb_runs/hht_couple/ou_seed777/FR_none_bestfit.txt

DONE 2) Λ(z) coupling / phase-locking probe

%run '/Users/boyde/.spyder-py3/hht_sn_residuals.py' \

'/Users/boyde/.spyder-py3/PLamb_Test_10V6c_plus.py' \

'/Users/boyde/Downloads/Pantheon+SH0ES.dat' \

'plamb_runs/hht_couple/ou_seed777'

[out] writing to: /Users/boyde/.spyder-py3/plamb_runs/hht_couple/ou_seed777

[clip] kept 0.98 of points.

[grid] uniform L bins = 1111, L ∈ [0.0012,1.1821]

[save] IMF stats -> /Users/boyde/.spyder-py3/plamb_runs/hht_couple/ou_seed777/hht_sn_imf_stats.csv, /Users/boyde/.spyder-py3/plamb_runs/hht_couple/ou_seed777/hht_sn_imf_stats.npz

[files] /Users/boyde/.spyder-py3/plamb_runs/hht_couple/ou_seed777/hht_sn.png

<Figure size 864x576 with 0 Axes>
