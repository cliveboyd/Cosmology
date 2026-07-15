# SU2 / Quaia Scan-Law + Colour Over-Absorption Probe

Date: 2026-07-16

## Purpose

This post-processes the completed 120k-row injection-recovery table to separate three effects:

1. catalogue-quality matched recovery,
2. physical recovery after scan-law+colour controls,
3. final BIC-gated scan-law+colour acceptance.

Strict physical recovery uses the original non-zero injection rule: SNR >= 3, angular separation <= 30 deg, and recovered/injected amplitude ratio in [0.5, 1.5].

Loose physical recovery is a diagnostic only: SNR >= 3, angular separation <= 45 deg, and recovered/injected amplitude ratio in [0.25, 2.5]. It is not proposed as a promotion gate.

The BIC gate is the current likelihood requirement Delta BIC = BIC(full dipole + controls) - BIC(controls only) < -10.

## Bottom Line

The scan-law+colour gate is over-rejecting injected angular modes at the current 1x amplitude. The strict physical recovery rate is 0.2833, but the final scan-law+colour pass rate is only 0.0011; 0.2823 of all trials are physically recovered and then rejected only by the BIC gate. The locked observed-direction injections show a separate failure mode: they are usually high-SNR after controls, but are rotated or amplitude-inflated enough to fail the strict recovery cuts.

## Amplitude Summary

| amp_scale | matched_pass_rate | scanlaw_strict_physical_rate | scanlaw_loose_physical_rate | scanlaw_colour_pass_rate | strict_physical_bic_rejected_rate | matched_to_strict_physical_loss_rate | delta_bic_p50 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0.25 | 0 | 0 | 0 | 0 | 0 | 0 | 26.4669 |
| 0.5 | 0.02375 | 0.00515 | 0.1989 | 0 | 0.00515 | 0.02235 | 23.2216 |
| 1 | 0.5012 | 0.28335 | 0.73445 | 0.0011 | 0.28225 | 0.2491 | 9.74463 |
| 1.5 | 0.80865 | 0.6278 | 0.9279 | 0.2694 | 0.3584 | 0.20915 | -11.0205 |
| 2 | 0.9135 | 0.85955 | 0.98845 | 0.6538 | 0.20575 | 0.0844 | -38.7547 |

## Direction-Family Summary

| category | amp_scale | matched_pass_rate | scanlaw_strict_physical_rate | scanlaw_loose_physical_rate | scanlaw_colour_pass_rate | strict_physical_bic_rejected_rate | matched_to_strict_physical_loss_rate | delta_bic_p50 | scanlaw_snr_p50 | scanlaw_sep_deg_p50 | scanlaw_amp_ratio_p50 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| cmb_or_anti | 0.25 | 0 | 0 | 0 | 0 | 0 | 0 | 23.1583 | 3.57517 | 71.75 | 3.55916 |
| cmb_or_anti | 0.5 | 0.0325 | 0.00075 | 0.182 | 0 | 0.00075 | 0.0325 | 22.6067 | 3.66379 | 55.6193 | 1.90222 |
| cmb_or_anti | 1 | 0.56425 | 0.096 | 0.515 | 0 | 0.096 | 0.47975 | 17.0351 | 4.35302 | 33.02 | 1.26885 |
| cmb_or_anti | 1.5 | 0.91 | 0.42275 | 0.7495 | 0.23575 | 0.187 | 0.497 | 4.01196 | 5.64238 | 21.5301 | 1.12561 |
| cmb_or_anti | 2 | 0.98925 | 0.80625 | 0.95625 | 0.434 | 0.37225 | 0.186 | -14.9173 | 7.14553 | 16.255 | 1.07569 |
| locked | 0.25 | 0 | 0 | 0 | 0 | 0 | 0 | 4.1956 | 5.63384 | 55.1399 | 5.35606 |
| locked | 0.5 | 0.00175 | 0 | 0.0775 | 0 | 0 | 0.00175 | -2.53832 | 6.20902 | 47.4724 | 2.9865 |
| locked | 1 | 0.12275 | 0.0335 | 0.798 | 0 | 0.0335 | 0.09875 | -19.4633 | 7.45239 | 35.8729 | 1.88218 |
| locked | 1.5 | 0.39275 | 0.25975 | 0.97525 | 0.1935 | 0.06625 | 0.194 | -43.4217 | 8.91693 | 28.132 | 1.53421 |
| locked | 2 | 0.63325 | 0.6355 | 0.9965 | 0.6215 | 0.014 | 0.127 | -71.4272 | 10.373 | 23.6102 | 1.37663 |
| random | 0.25 | 0 | 0 | 0 | 0 | 0 | 0 | 30.574 | 2.30193 | 54.3618 | 2.58861 |
| random | 0.5 | 0.0281667 | 0.00833333 | 0.245 | 0 | 0.00833333 | 0.0258333 | 27.3921 | 2.9061 | 36.9515 | 1.61624 |
| random | 1 | 0.606333 | 0.429083 | 0.786417 | 0.00183333 | 0.42725 | 0.222333 | 15.4708 | 4.52911 | 21.93 | 1.22224 |
| random | 1.5 | 0.9135 | 0.818833 | 0.971583 | 0.305917 | 0.512917 | 0.11825 | -4.11005 | 6.33606 | 15.2453 | 1.12003 |
| random | 2 | 0.981667 | 0.952 | 0.9965 | 0.737833 | 0.214167 | 0.0363333 | -32.1916 | 8.26891 | 11.8667 | 1.08168 |

## Locked-Direction 1x Breakdown

| direction_label | bcut_deg | matched_pass_rate | scanlaw_strict_physical_rate | scanlaw_loose_physical_rate | scanlaw_colour_pass_rate | strict_physical_bic_rejected_rate | high_snr_rate | amp_too_high_rate | sep_too_high_rate | delta_bic_p50 | scanlaw_snr_p50 | scanlaw_sep_deg_p50 | scanlaw_amp_ratio_p50 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| locked_b15 | 10 | 0 | 0.002 | 0.766 | 0 | 0.002 | 1 | 0.998 | 0.882 | -47.2377 | 9.15164 | 37.6291 | 2.14131 |
| locked_b15 | 15 | 0.042 | 0.014 | 0.902 | 0 | 0.014 | 1 | 0.936 | 0.656 | -22.402 | 7.66605 | 33.4016 | 1.85412 |
| locked_b15 | 25 | 0.17 | 0.032 | 0.866 | 0 | 0.032 | 1 | 0.848 | 0.616 | -13.5179 | 7.02493 | 32.7205 | 1.78797 |
| locked_b15 | 35 | 0.296 | 0.07 | 0.792 | 0 | 0.07 | 0.998 | 0.702 | 0.584 | -0.455736 | 5.95802 | 32.876 | 1.67006 |
| locked_b25 | 10 | 0 | 0 | 0.592 | 0 | 0 | 1 | 1 | 0.978 | -57.709 | 9.70711 | 41.4754 | 2.24738 |
| locked_b25 | 15 | 0.016 | 0.002 | 0.85 | 0 | 0.002 | 1 | 0.948 | 0.822 | -27.7862 | 8.00961 | 36.966 | 1.90337 |
| locked_b25 | 25 | 0.178 | 0.052 | 0.866 | 0 | 0.052 | 1 | 0.788 | 0.594 | -10.5156 | 6.80784 | 32.2863 | 1.73094 |
| locked_b25 | 35 | 0.28 | 0.096 | 0.75 | 0 | 0.096 | 0.996 | 0.58 | 0.62 | 2.95869 | 5.66418 | 34.0527 | 1.56699 |

## Largest 1x BIC False-Negative Conditions

| direction_label | bcut_deg | matched_pass_rate | scanlaw_strict_physical_rate | scanlaw_loose_physical_rate | scanlaw_colour_pass_rate | strict_physical_bic_rejected_rate | high_snr_rate | amp_too_high_rate | sep_too_high_rate | delta_bic_p50 | scanlaw_snr_p50 | scanlaw_sep_deg_p50 | scanlaw_amp_ratio_p50 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| random_01 | 10 | 0.828 | 0.756 | 0.914 | 0 | 0.756 | 0.924 | 0.082 | 0.12 | 18.0366 | 4.29704 | 16.3705 | 1.13787 |
| random_05 | 15 | 0.808 | 0.72 | 0.864 | 0 | 0.72 | 0.884 | 0.016 | 0.208 | 18.0347 | 4.28083 | 19.4959 | 1.00983 |
| random_02 | 10 | 0.888 | 0.704 | 0.99 | 0.002 | 0.702 | 0.992 | 0.244 | 0.05 | 3.77149 | 5.72108 | 13.6948 | 1.32381 |
| random_01 | 15 | 0.768 | 0.686 | 0.798 | 0 | 0.686 | 0.804 | 0.042 | 0.122 | 21.9571 | 3.79512 | 17.5212 | 1.04883 |
| random_05 | 25 | 0.728 | 0.654 | 0.928 | 0 | 0.654 | 0.96 | 0.118 | 0.224 | 13.4919 | 4.72604 | 20.476 | 1.19677 |
| random_04 | 25 | 0.732 | 0.624 | 0.878 | 0 | 0.624 | 0.904 | 0.102 | 0.252 | 17.0254 | 4.3361 | 22.1268 | 1.09815 |
| random_05 | 10 | 0.776 | 0.606 | 0.77 | 0 | 0.606 | 0.81 | 0.004 | 0.28 | 21.3799 | 3.8886 | 22.5421 | 0.917058 |
| random_06 | 10 | 0.726 | 0.602 | 0.726 | 0 | 0.602 | 0.734 | 0.02 | 0.22 | 23.8591 | 3.55555 | 19.0793 | 0.93542 |
| random_06 | 15 | 0.742 | 0.592 | 0.732 | 0 | 0.592 | 0.744 | 0.052 | 0.142 | 22.6183 | 3.70698 | 16.9368 | 1.02537 |
| random_02 | 15 | 0.848 | 0.544 | 0.994 | 0 | 0.544 | 1 | 0.416 | 0.05 | -1.12054 | 6.12229 | 13.3875 | 1.45437 |
| random_06 | 25 | 0.616 | 0.432 | 0.624 | 0 | 0.432 | 0.638 | 0.134 | 0.19 | 24.6817 | 3.33844 | 17.9559 | 1.11499 |
| random_01 | 25 | 0.568 | 0.41 | 0.594 | 0 | 0.41 | 0.616 | 0.1 | 0.264 | 24.7646 | 3.32601 | 20.978 | 1.06311 |

## Interpretation

For random injected directions at 1x amplitude, scan-law+colour often still recovers a recognisable physical dipole, but the BIC term rejects almost all of those recoveries. This says the current BIC threshold is too stringent for a promotion gate unless calibrated by injection-recovery power.

For locked observed-direction injections at 1x amplitude, the problem is not mainly a BIC penalty. The full controlled fit is high-SNR, but the direction and amplitude are distorted: median separations are commonly just above 30 deg and median amplitude ratios are often above 1.5. That points to collinearity or template leakage between the locked angular mode and the scan-law+colour controls.

The follow-up gate should therefore use a two-part diagnostic: first calibrate detection power and false-negative rate with injections, then inspect template collinearity for the locked z ~ 1-1.5 mode before treating a non-pass as a physical rejection.
