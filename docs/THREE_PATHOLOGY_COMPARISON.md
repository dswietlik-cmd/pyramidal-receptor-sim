# Three-mechanism pathology comparison (v0.9.5)

## States

The v0.9.5 comparison includes four states:

1. **Healthy**: reference configuration.
2. **Input-driven hyperexcitability**: 1.65x excitatory input drive; calibrated moderate phenotype.
3. **NMDA-driven excitotoxicity-like state**: 4.10x NMDA pathology multiplier; calibrated moderate phenotype.
4. **GABA-A complete-loss boundary**: 0x GABA-A pathology efficacy; failed to achieve a matched hyperexcitable phenotype and is therefore treated as a stress-test rather than a matched pathology.

All simulations use n=20 paired seeds and identical drug concentration grids.

## Baseline comparison

| State | Mean firing (Hz) | SD (Hz) | Mean paired delta vs healthy (Hz) |
|---|---:|---:|---:|
| Healthy | 10.965 | 0.906 | 0.000 |
| Input-driven hyperexcitability | 24.605 | 1.089 | 13.640 |
| NMDA-driven excitotoxicity-like | 23.605 | 5.616 | 12.640 |
| GABA-A complete loss | 11.280 | 0.738 | 0.315 |

The first two pathological mechanisms therefore provide approximately matched moderate hyperexcitability with different mechanisms. Isolated GABA-A loss does not.

## Pharmacological stress test

At concentrations that preserved at least ~80% of healthy firing, the largest observed within-state suppression was approximately:

| State | Perampanel | Memantine | Diazepam |
|---|---:|---:|---:|
| Input-driven | 7.8% at 8 nM | 2.6% at 50 nM | 1.5% at 5 nM |
| NMDA-driven | 12.7% at 8 nM | 7.7% at 50 nM | 0.7% at 5 nM |
| GABA-A complete-loss boundary | 13.0% at 8 nM | 5.7% at 50 nM | 0% at 0 nM |

These are model-specific functional results. They are not clinical therapeutic windows.

## Interpretation

The comparison strengthens two conclusions for the current pre-publication model:

- Drug response depends on the mechanism generating the abnormal state, not only on baseline firing rate.
- The current GABA-A branch is too weak at baseline for isolated loss of inhibition to reproduce a matched hyperexcitable phenotype. That limitation should be reported rather than hidden by additional compensatory parameter changes.

The GABA-A finding is therefore a calibration boundary that informs future model refinement.
