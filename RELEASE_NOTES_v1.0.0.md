# v1.0.0 - Manuscript-associated frozen release

Version 1.0.0 freezes the computational model and manuscript-associated analysis workflow used for the final paper package.

## Included in the frozen release

- CA1-like pyramidal-neuron model with AMPA, NMDA, and GABA-A receptor components.
- Paired n=20 healthy, input-driven, and NMDA-driven concentration-response experiments.
- State-specific 4PL analysis with 500 paired-seed bootstrap resamples.
- Friedman repeated-measures tests, Holm-corrected paired Wilcoxon tests, Kendall's W, and matched rank-biserial effect sizes.
- Publication Tables 1-5 and 9-11 plus Supplementary Tables S3-S5.
- Reproducible English and Polish publication figures in SVG, PDF, and 800-dpi PNG formats.
- Language-specific figure directories under `results/publication_final/figures_professional/`.
- Final manuscript/repository alignment for Supplementary Tables S4 and S5.
- Corrected sensitivity-analysis and reproduction-pipeline paths.
- Updated reproducibility and manuscript-workflow documentation.

## Scientific freeze

No neuronal equations, pharmacodynamic parameters, pathology definitions, concentration grids, paired random seeds, primary concentration-response data, state-specific 4PL results, or final OAT sensitivity data are changed relative to the finalized v0.9.9 manuscript analysis.

The v1.0.0 tag should identify the exact software state cited by the manuscript. After release, archive the tagged release in Zenodo (or another long-term repository) and cite the resulting DOI in the manuscript.
