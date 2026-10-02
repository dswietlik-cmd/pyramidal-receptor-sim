# v0.9.8 - Reproducibility and manuscript-pipeline consistency update

- Aligns package metadata, `__version__`, `CITATION.cff`, and README development status with v0.9.8.
- Replaces the legacy `reproduce_paper.py` path with the final paired healthy-vs-pathology manuscript pipeline.
- Ensures the reproduction entry point uses the final state-specific publication analysis with 500 paired-seed bootstrap resamples.
- Routes sensitivity regeneration to `run_sensitivity_oat_final_fast.py`, which reads the final healthy-state functional EC50 values directly from the final 4PL table.
- Adds a `publication` optional dependency group for `pandas`, `matplotlib`, and `numba`, and aligns `requirements.txt` with the full manuscript workflow.
- Updates the reproducibility documentation and README instructions to the final paired n=20 experiment.
- No neuronal equations, pharmacodynamic parameters, concentration grids, pathology definitions, primary results, or final manuscript numerical results are changed by this release.

This remains a pre-publication research release.
