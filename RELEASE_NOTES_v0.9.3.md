# v0.9.3 - Healthy vs pathological pharmacology

This pre-publication release adds the first paired healthy-versus-moderate-hyperexcitability concentration-response experiment.

## Added

- paired n=20 healthy and moderate-hyperexcitability simulations;
- perampanel, memantine, and diazepam concentration-response grids in both states;
- Pathological Activity Normalization (PAN);
- exploratory Therapeutic Selectivity Index (TSI);
- extended receptor/synaptic outputs;
- a Numba-accelerated simulator checked against the reference Python model;
- `docs/HEALTHY_PATHOLOGICAL_DRUGS.md`;
- reference CSV outputs for the paired experiment.

## Interpretation

Within the present pathology definition (1.65x excitatory input drive), receptor modulation suppresses the healthy near-threshold state more readily than the moderate hyperexcitability state. This is a model result and not a clinical efficacy claim.
