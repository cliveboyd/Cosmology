# FR LINX Anchor Sensitivity Reweight

This is a likelihood-only reweight of the completed abundance catalogue. It does not rerun LINX, add nuclear theory covariance, or replace the registered gate.

The combined-current case uses D/H x 1e5 = `2.533 +/- 0.024`, Yp = `0.2458 +/- 0.0013`, and neutron lifetime = `878.4 +/- 0.5 s`. Lithium remains `1.45 +/- 0.25` in units of Li/H x 1e10.

## Combined-Current Readout

| Scenario | Rows | D+He pass | D+He+Li pass | Best reweighted chi2 | Min Li given D+He |
|---|---:|---:|---:|---:|---:|
| su2_expansion_only | 7 | 0 | 0 | 281.7 | NA |
| su2_plus_modest_rate_controls | 906 | 99 | 0 | 230.5 | 5.122 |
| su2_plus_scanned_rate_controls | 1088 | 124 | 0 | 226.6 | 5.009 |
| all_scanned_controls | 30503 | 2872 | 0 | 214 | 4.741 |

## Interpretation

The updated anchors are a sensitivity analysis, not a new prior default for an FR cosmology. They retain the standard-background baryon treatment and do not marginalise the fixed deuterium-sensitive reaction pulls. The table should therefore be read as the direction and scale of likelihood sensitivity only.
