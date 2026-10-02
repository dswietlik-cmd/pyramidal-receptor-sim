# v0.9.9 - Inferential statistics and publication-figure reproducibility

- Adds `scripts/run_stat_analysis.py` for the final paired n=20 manuscript analysis.
- Adds Friedman repeated-measures tests, Holm-corrected paired Wilcoxon tests, Kendall's W,
  and matched rank-biserial effect sizes.
- Adds version-controlled statistical outputs: Tables 9-11 and Supplementary Tables S4-S5.
- Adds `scripts/make_professional_figures.py` for reproducible journal-ready figures.
- Professional figures are exported as vector SVG/PDF and 800-dpi PNG; optional 800-dpi
  TIFF export is available with `--tiff`.
- Integrates the statistical analysis and professional figure generation into
  `scripts/reproduce_paper.py`.
- Adds `docs/STATISTICAL_ANALYSIS.md` and updates the manuscript reproduction documentation.
- No neuronal equations, pharmacodynamic parameters, pathology definitions, concentration
  grids, random seeds, primary concentration-response data, 4PL fits, or sensitivity data are
  changed by this release.

The inferential analyses use computational stochastic replicates and are not intended as
biological population inference. This remains a pre-publication research release.
