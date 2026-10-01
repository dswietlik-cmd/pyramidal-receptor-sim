# Final experiment specification

- States: healthy; input-driven moderate hyperexcitability (1.65x excitatory drive); NMDA-driven moderate excitotoxicity-like state (4.10x NMDA pathology multiplier).
- Replicates: n=20 paired seeds per state and concentration.
- Base seed: 20261200.
- Input jitter: 0.5 ms.
- Drugs: perampanel, memantine, diazepam.
- Concentration grids are defined in scripts/compare_two_pathologies_drugs_fast.py.
- Primary functional outcome: firing rate.
- Secondary model outcomes: ISI, AMPA/NMDA/GABA PSP components, NMDA gate-open fraction, and mean plasticity state.
- Pathological Activity Normalization (PAN) is computed against same-seed healthy and pathological baselines.
- A descriptive healthy-preserving window is defined on the tested concentration grid as concentrations preserving >=80% of mean healthy baseline firing.
- State-specific functional concentration-response curves are summarized with a 4-parameter logistic model and 500 paired-seed bootstrap resamples.

See `results/publication_final/Results_final_n20_two_pathologies.md` for the current Results draft.
