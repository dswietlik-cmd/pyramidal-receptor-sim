# GABA-A disinhibition / hypofunction boundary analysis (v0.9.5)

## Purpose

Version 0.9.5 tested whether isolated loss of GABA-A-mediated inhibitory efficacy can generate a hyperexcitable phenotype comparable with the previously calibrated input-driven and NMDA-driven moderate pathological states.

The disease parameter is `gabaa_pathology_multiplier`, which is intentionally separate from the pharmacological `gabaa_gain` used by diazepam. During a simulation, the two factors combine multiplicatively at GABA-A synaptic input insertion.

- `gabaa_pathology_multiplier = 1.0`: healthy reference inhibitory efficacy.
- `gabaa_pathology_multiplier = 0.0`: complete loss of modeled GABA-A inhibitory efficacy.
- Diazepam remains a separate concentration-dependent positive allosteric modulation in `gabaa_gain`.

No AMPA gain, NMDA gain, input frequency, spike threshold, or synaptic amplitude is changed by this pathology parameter.

## Calibration protocol

The same 20 paired jittered input realizations used elsewhere in the pre-publication workflow were evaluated across a coarse GABA-A pathology grid from 1.00 to 0.00 in steps of 0.05. Each run lasted 10 s with `dt = 0.5 ms` and `jitter = 0.5 ms`.

Pre-specified phenotype targets were approximately 18 Hz, 24 Hz, and 30 Hz, matching the mild, moderate, and severe hyperexcitability targets used for the other pathological mechanisms.

## Result

Isolated GABA-A hypofunction did **not** produce the required hyperexcitable phenotype in the present reference model.

| GABA-A pathology multiplier | Mean firing (Hz) | SD (Hz) | Mean GABA PSP component (mV) |
|---:|---:|---:|---:|
| 1.00 | 10.965 | 0.906 | -0.549 |
| 0.75 | 11.085 | 0.819 | -0.368 |
| 0.50 | 11.215 | 0.800 | -0.235 |
| 0.25 | 11.230 | 0.803 | -0.117 |
| 0.00 | 11.280 | 0.738 | 0.000 |

Even complete modeled loss of GABA-A inhibition increased mean firing by only about 0.315 Hz relative to the paired healthy baseline. Therefore no mild/moderate/severe GABA-A disinhibition presets are declared in this release.

## Interpretation

This is a model-identification result, not evidence that biological GABAergic disinhibition is unimportant. It shows that, under the current CA1-like reference calibration, inhibitory GABA-A input contributes too little to the operating-point firing rate for its removal alone to reproduce the matched 18/24/30 Hz phenotypes.

Forcing the GABA-A branch to the target firing rates would require changing additional model parameters and would therefore confound the intended one-mechanism pathology comparison.

Accordingly, v0.9.5 retains complete GABA-A loss only as a **boundary/stress-test condition**, not as a matched moderate pathology.

## Consequence for drug comparisons

The full n=20 drug screen includes the GABA-A complete-loss boundary condition alongside healthy, input-driven hyperexcitability, and NMDA-driven excitotoxicity-like states. Because its baseline firing differs from healthy by less than 2 Hz, Pathological Activity Normalization (PAN) and the exploratory therapeutic selectivity index are marked not estimable for this boundary state.

Diazepam cannot restore a GABA-A response when the pathology multiplier is exactly zero because the disease and drug factors combine multiplicatively. This is expected behavior for the mathematical formulation and should not be interpreted as a clinical prediction.

## Reproduce

```bash
pip install -e ".[fast]"
python scripts/calibrate_gaba_disinhibition_fast.py
python scripts/compare_three_pathologies_drugs_fast.py
```

Reference outputs are stored in `results/reference/`.
