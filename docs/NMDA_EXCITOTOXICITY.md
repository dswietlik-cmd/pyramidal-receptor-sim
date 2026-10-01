# NMDA-driven excitotoxicity-like pathology

## Purpose

This pathology is a mechanistic comparator to the input-driven hyperexcitability phenotype.
It increases NMDA-mediated synaptic drive while keeping presynaptic event frequencies,
AMPA function, GABA-A function, membrane thresholds, and pharmacological receptor gains
at their calibrated reference values.

The disease-level multiplier (`nmda_pathology_multiplier`) is independent of the drug-level
`nmda_gain`. This separation allows memantine to act on the NMDA target without defining
the pathological state by the pharmacological parameter itself.

## Implementation

For an excitatory event occurring while the NMDA gate is open, the inserted NMDA kernel is

`NMDA_input = nmda_pathology_multiplier * nmda_gain * NMDA_kernel`.

Healthy reference: `nmda_pathology_multiplier = 1.0`.

The model remains voltage gated through the calibrated NMDA activation threshold and the
explicit NMDA postsynaptic contribution remains coupled through `nmda_psp_scale`.

## Calibration (n=20 paired seeds)

The selected states were calibrated against the same target firing phenotypes used for the
input-driven pathology:

| State | NMDA pathology multiplier | Target | Mean firing ± SD |
|---|---:|---:|---:|
| Healthy | 1.00 | ~12 Hz | 11.015 ± 0.824 Hz |
| Mild | 3.00 | ~18 Hz | 18.530 ± 1.931 Hz |
| Moderate | 4.10 | ~24 Hz | 23.405 ± 5.556 Hz |
| Severe | 4.70 | ~30 Hz | 30.375 ± 6.147 Hz |

The moderate state is the main comparator to the previously calibrated input-driven moderate
hyperexcitability state (~24.6 Hz). Thus, the two pathologies have similar mean output firing
but differ in mechanism.

## Interpretation

This is an **excitotoxicity-like computational phenotype**, not a direct measurement of
extracellular glutamate, calcium concentration, neuronal injury, or clinical disease severity.
The raw plasticity state rises by orders of magnitude as NMDA drive increases; therefore,
its absolute magnitude should not be interpreted as a biological concentration. It is an
internal state variable whose logarithm contributes to synaptic gain.

The historical publications used a different phenomenological excitotoxicity formalism based
on increasing the `powerA/powerB` parameter. The new pathology does not numerically map
those historical values into the current refactored model; it preserves the same biological
concept (increased NMDA-associated excitotoxic drive) while keeping the new receptor-level
architecture explicit and independently drug-modulatable.
