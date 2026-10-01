# v0.9.0 — Pre-publication research version

Public pre-publication version of **pyramidal-receptor-sim**.

This version is intended for transparency, collaborative inspection, and reproducible development while the associated manuscript is still being prepared. It is **not** the frozen manuscript-associated software release.

## Included

- CA1-like pyramidal neuron model with AMPA, explicit NMDA, and GABA-A components
- concentration-to-receptor mappings for memantine, perampanel, and diazepam
- deterministic and jittered synaptic input generators
- scripts for the current n=20 concentration-response workflow
- 4PL fitting workflow
- calibration, pharmacology, model, and reproducibility documentation
- automated tests and GitHub Actions configuration

## Reproducibility status

The current workflow is reproducible for the code and parameter set in this tag. However, model equations, parameterization, concentration grids, analysis scripts, and manuscript endpoints may still change before v1.0.0.

## Archival plan

Do **not** treat v0.9.0 as the final software citation for the manuscript. After the model and results are frozen, a manuscript-associated **v1.0.0** release will be created and archived with Zenodo (or another long-term repository) to obtain a DOI.

## License

BSD 3-Clause License. See `LICENSE`.
