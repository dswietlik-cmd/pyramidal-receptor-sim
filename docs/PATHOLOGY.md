# Pathological hyperexcitability phenotype

## Rationale

The main pathological state is defined by **increased excitatory synaptic input drive**, while calibrated receptor gains, inhibitory input frequencies, membrane parameters, and synaptic amplitudes remain unchanged. This separates the disease phenotype from the receptor targets later manipulated pharmacologically.

The calibration used the same 10 s simulation duration, 0.5 ms time step, bounded ±0.5 ms jitter, and 20 paired seeds used in the publication workflow. Only the 13 excitatory input frequencies were multiplied; the 3 inhibitory inputs were unchanged.

## Calibrated states

| State | Excitatory drive multiplier | Target firing | Calibrated firing, mean ± SD (n=20) |
|---|---:|---:|---:|
| Healthy reference | 1.00 | ~12 Hz | 10.965 ± 0.906 Hz |
| Mild hyperexcitability | 1.30 | ~18 Hz | 18.055 ± 0.199 Hz |
| Moderate hyperexcitability | 1.65 | ~24 Hz | 24.605 ± 1.089 Hz |
| Severe hyperexcitability | 1.90 | ~30 Hz | 29.640 ± 0.289 Hz |

The pathological labels describe **phenotype severity in this computational model**, not clinical disease stages.

## Calibration outputs

`results/reference/hyperexcitability_calibration_replicates.csv` contains all 620 simulations (31 drive levels × 20 seeds).

`results/reference/hyperexcitability_calibration_summary.csv` contains mean/SD values for firing, ISI, AMPA/NMDA/GABA PSP components, NMDA gate-open fraction, and the model plasticity state.

`results/reference/hyperexcitability_selected_states.csv` contains the selected mild/moderate/severe calibration points.

## Important observation

The firing-rate response to input-drive scaling is not perfectly monotonic at every intermediate multiplier. This is expected from the nonlinear threshold dynamics and temporally jittered multi-input structure. Therefore state selection was based on the closest n=20 mean firing rate to each prespecified phenotype target rather than assuming a monotonic linear mapping.
