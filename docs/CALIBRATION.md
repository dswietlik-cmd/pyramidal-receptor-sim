# Calibration and sensitivity analysis

## Calibration layers

The repository separates three calibration layers:

1. **Baseline neuronal state**: membrane, synaptic and input parameters defining the drug-free model.
2. **Receptor scaling**: relative AMPA, explicit NMDA and GABA-A contributions to the synaptic response.
3. **Drug pharmacodynamics**: concentration-response mappings that convert a drug concentration into a receptor gain.

These layers should not be conflated. A change in baseline excitability can alter the functional
drug-response curve even when the molecular concentration-response mapping is unchanged.

## AMPA

The reference AMPA EPSP amplitude is treated as an effective clustered Schaffer-collateral-like input.
The current reference value is stored in `configs/ca1_reference.txt` as `AMPA_EPSP_AMPLITUDE_MV`.

## NMDA

The NMDA register contributes explicitly to the somatic postsynaptic potential through
`NMDA_PSP_SCALE`. The reference value was chosen to anchor the native model's NMDA/AMPA kernel
relationship to an experimental CA1 NMDA/AMPA response ratio while preserving voltage-dependent
NMDA gating through `NMDA_ACTIVATION_THRESHOLD_MV`.

This calibration does not imply that the NMDA contribution is constant across membrane potentials.
The modeled magnesium-block mechanism continues to limit NMDA activation at hyperpolarized states.

## GABA-A

The inhibitory amplitude is calibrated as an effective compound inhibitory input relative to the
reference excitatory amplitude. It should not be interpreted as the amplitude of a single unitary
GABAergic synapse.

## Drug concentration-response parameters

The literature-anchored values used by the current release are defined only in `drugs.py`. Keeping
these parameters in a separate module makes the assumptions auditable and prevents hidden
pharmacological constants inside the neuronal equations.

## Functional EC50

The publication workflow fits a four-parameter logistic curve to firing suppression. The resulting
functional EC50 is an emergent property of the full configured model and depends on:

- receptor scaling;
- `SPIKE_THRESHOLD_MV`;
- `NMDA_ACTIVATION_THRESHOLD_MV`;
- E/I balance;
- input timing and synchrony;
- the pharmacodynamic mapping.

It is therefore not interchangeable with molecular IC50/EC50, receptor affinity, clinical dose, or
measured human exposure.

## Sensitivity analysis

At minimum, publication analyses should examine sensitivity to:

- AMPA/EPSP strength;
- GABA-A/IPSP strength;
- explicit NMDA coupling (`NMDA_PSP_SCALE`);
- NMDA activation threshold (`NMDA_ACTIVATION_THRESHOLD_MV`);
- spike threshold (`SPIKE_THRESHOLD_MV`).

The reference outputs under `results/reference/` should be regenerated after any change to the
model equations, configuration, input generator, or pharmacodynamic parameters.
