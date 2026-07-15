# SU2 / Quaia Scan-Law + Colour Over-Absorption Probe

Date: 2026-07-16

## Bottom Line

The scan-law+colour gate is over-rejecting injected angular modes at the current 1x amplitude. The strict physical recovery rate is 0.2833, but the final scan-law+colour pass rate is only 0.0011; 0.2823 of all trials are physically recovered and then rejected only by the BIC gate. The locked observed-direction injections show a separate failure mode: they are usually high-SNR after controls, but are rotated or amplitude-inflated enough to fail the strict recovery cuts.

## 1x Gate Decomposition

- matched-quality recovery rate: `0.5012`
- strict scan-law+colour physical recovery rate: `0.2833`
- loose scan-law+colour physical recovery rate: `0.7345`
- final scan-law+colour pass rate: `0.0011`
- strict physical recovery rejected only by BIC: `0.2823`
- matched recovery lost before strict physical scan-law recovery: `0.2491`

## 1x By Direction Family

| category | amp_scale | matched_pass_rate | scanlaw_strict_physical_rate | scanlaw_loose_physical_rate | scanlaw_colour_pass_rate | strict_physical_bic_rejected_rate | matched_to_strict_physical_loss_rate | delta_bic_p50 | scanlaw_snr_p50 | scanlaw_sep_deg_p50 | scanlaw_amp_ratio_p50 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| cmb_or_anti | 1 | 0.56425 | 0.096 | 0.515 | 0 | 0.096 | 0.47975 | 17.0351 | 4.35302 | 33.02 | 1.26885 |
| locked | 1 | 0.12275 | 0.0335 | 0.798 | 0 | 0.0335 | 0.09875 | -19.4633 | 7.45239 | 35.8729 | 1.88218 |
| random | 1 | 0.606333 | 0.429083 | 0.786417 | 0.00183333 | 0.42725 | 0.222333 | 15.4708 | 4.52911 | 21.93 | 1.22224 |
