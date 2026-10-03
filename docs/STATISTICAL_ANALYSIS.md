# Statistical analysis for the final manuscript experiment

## Scope

The final experiment uses **20 paired computational stochastic replicates** (seeds
`20261200`-`20261219`). The same seed ID is reused across states and concentrations.
Statistical procedures therefore preserve pairing rather than treating simulations as
independent observations.

These replicates are **not biological samples**. P-values quantify the consistency of
responses across the specified simulation ensemble and should not be generalized to
animal, patient, or experimental-population inference.

## Baseline-state comparison

Baseline firing rates for healthy, input-driven, and NMDA-driven states are compared with
a Friedman repeated-measures rank test. Kendall's W is reported as the global effect-size
measure.

If the global test is significant, the three paired state contrasts are evaluated with
two-sided Wilcoxon signed-rank tests and Holm correction for family-wise error. Matched
rank-biserial correlation is reported as an effect-size measure.

## Concentration-response inference

For each drug and state, firing rate across the complete tested concentration grid is tested
with a Friedman repeated-measures rank test. Kendall's W is reported.

Planned comparisons of each non-zero concentration against the same-seed zero-concentration
condition use two-sided Wilcoxon signed-rank tests. Holm correction is applied separately
within each drug-state family. These tests complement, but do not replace, the 4PL
concentration-response analysis.

## Mechanism-dependent pathology contrast

At each non-zero concentration, pathological suppression in the NMDA-driven state is
compared with pathological suppression in the input-driven state using a paired Wilcoxon
signed-rank test matched by seed. Holm correction is applied within each drug across its
concentration grid. Matched rank-biserial correlation and paired suppression differences are
reported.

## 4PL uncertainty

The state-specific 4PL fits remain separate from the rank-based inferential tests. Functional
EC50 uncertainty is estimated with **500 paired-seed bootstrap resamples**, sampling whole
seed profiles with replacement. The bootstrap RNG seed is `20261001`, and percentile
95% confidence intervals use the 2.5th and 97.5th percentiles.

## Output tables

Running:

```bash
python scripts/run_stat_analysis.py
```

writes:

- `results/publication_final/Table9_baseline_inferential_statistics.csv`
- `results/publication_final/Table10_concentration_global_statistics.csv`
- `results/publication_final/Table11_pathology_state_contrast_summary.csv`
- `results/publication_final/TableS4_concentration_vs_zero_wilcoxon.csv`
- `results/publication_final/TableS5_pathology_state_suppression_contrasts.csv`

The script reads the paired replicate-level data in `results/reference/`.
