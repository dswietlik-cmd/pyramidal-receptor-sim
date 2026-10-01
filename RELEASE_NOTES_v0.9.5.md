# v0.9.5 - GABA-A disinhibition boundary analysis

Pre-publication development release.

## Added

- A pathology-specific `gabaa_pathology_multiplier` independent of diazepam `gabaa_gain`.
- n=20 calibration of isolated GABA-A hypofunction/disinhibition from 1.0 to complete loss.
- Three-mechanism pharmacology comparison including healthy, input-driven, NMDA-driven, and GABA-A complete-loss states.
- `docs/GABAA_DISINHIBITION.md` and `docs/THREE_PATHOLOGY_COMPARISON.md`.
- Tests for separation and bounds of GABA-A pathology versus drug modulation.

## Key finding

Under the current CA1-like reference calibration, isolated removal of GABA-A inhibition does not generate the pre-specified 18/24/30 Hz hyperexcitability phenotypes. Mean firing changes from 10.965 Hz at the healthy reference to only 11.280 Hz after complete modeled GABA-A loss. The complete-loss state is therefore retained only as a boundary/stress-test condition and is not labeled as a matched moderate pathology.

## Validation

13 automated tests pass in the local validation environment.
