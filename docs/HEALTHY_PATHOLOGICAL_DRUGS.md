# Healthy vs moderate-hyperexcitability pharmacology

This analysis compares the same receptor pharmacology in two calibrated model states:

- **healthy:** excitatory drive multiplier = 1.00;
- **moderate hyperexcitability:** excitatory drive multiplier = 1.65.

No receptor gain, inhibitory input frequency, synaptic amplitude, or membrane threshold is changed to define the pathological state. The purpose is to test whether the pharmacological response depends on neuronal state rather than being built into the pathology definition.

## Design

For each state, drug, concentration, and replicate:

- 20 independent seed-defined input realizations are used;
- the same seed identifiers are retained across healthy and pathological states;
- the same jitter parameter (0.5 ms) and 10 s simulation duration are used;
- concentration is converted to receptor gain by the receptor pharmacodynamic model before simulation.

Drug grids (nM):

- perampanel: 0, 6, 8, 10, 12, 14, 16, 18, 22;
- memantine: 0, 50, 100, 125, 150, 175, 200, 250, 350, 500;
- diazepam: 0, 2.5, 5, 7.5, 10, 12.5, 15, 20, 30, 60.

## Endpoints

The main paired endpoint is **Pathological Activity Normalization (PAN)**:

```text
PAN = 1 - |F_path,drug - F_healthy,0| / |F_path,0 - F_healthy,0|
```

PAN = 1 indicates that pathological firing has been restored exactly to the drug-free healthy firing level. PAN = 0 indicates no movement toward the healthy reference. Negative values indicate movement farther away or marked overshoot beyond the healthy reference.

An exploratory **Therapeutic Selectivity Index (TSI)** is also reported:

```text
TSI = pathological suppression / (max(healthy suppression, 0) + 0.01)
```

TSI is an exploratory model index, not a clinical therapeutic index.

Secondary outputs include firing rate, mean ISI, AMPA/NMDA/GABA-A postsynaptic components, NMDA-open fraction, and mean plasticity state.

## Main result of the current model

The moderate hyperexcitability state remains comparatively resistant to receptor modulation, while the healthy near-threshold state is suppressed more strongly. In the tested ranges, no drug produced strong normalization of moderate hyperexcitability while simultaneously preserving >=80% of healthy firing.

Under the >=80% healthy-firing-preservation criterion, the largest observed pathological suppression was approximately:

- perampanel 8 nM: 7.6% pathological suppression;
- memantine 50 nM: 2.5% pathological suppression;
- diazepam 5 nM: 1.4% pathological suppression.

This result should be treated as a mechanistic property of the present model definition, not as a clinical efficacy statement. It indicates that hyperexcitability defined exclusively by increased excitatory input drive can be less drug-sensitive than the healthy state in this model.

## Reproducibility

Run:

```bash
pip install -e ".[fast]"
python scripts/run_healthy_pathological_drugs_fast.py
```

The script verifies the Numba implementation against the reference Python implementation for representative drug/state conditions before executing the full experiment.
