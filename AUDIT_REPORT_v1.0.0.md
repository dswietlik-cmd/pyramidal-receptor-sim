# Final pre-release audit for v1.0.0

## Verified

- Python test suite: 19/19 tests pass.
- Default manuscript reproduction pipeline completes successfully.
- Final OAT sensitivity analysis completes successfully; regenerated Table S3 is byte-identical to the committed final Table S3 (SHA256 `5c01b9d2927fa9bf3fc2f7afc528e42b9135d7a78eea313c1b7a44b093093cd3`).
- Package version declarations are aligned at 1.0.0 in `pyproject.toml`, `src/pyramidal_receptor_sim/__init__.py`, and `CITATION.cff`.
- Canonical publication tables include Tables 1-5, 9-11, and Supplementary Tables S3-S5.
- English and Polish professional figure directories contain Figure 1-5 and Figure S1 in PNG, SVG, and PDF formats.
- The old duplicate `results/publication_final/en/` and `results/publication_final/pl/` directories are absent.
- The misplaced root sensitivity table and sensitivity script are absent; their canonical locations are under `results/publication_final/` and `scripts/`.
- The default reproduction pipeline does not recreate legacy Figures 1-5 in the top level of `results/publication_final/`.
- Active documentation paths are aligned with `results/publication_final/figures_professional/en/` and `/pl/`.
- No textual markers such as `ChatGPT`, `OpenAI`, `AI-generated`, or similar were found in the release-candidate repository files searched.
- The Polish manuscript was rendered to 35 pages and visually checked after correcting Supplementary Tables S4/S5.

## Scientific consistency

- Supplementary Table S4 is the within-drug/state paired comparison of each non-zero concentration versus 0 nM.
- Supplementary Table S5 is the paired NMDA-driven versus input-driven suppression contrast at matched concentrations.
- The manuscript cross-reference to pathology-state suppression contrasts points to Table S5.
- No neuronal equations, pharmacodynamic parameters, pathology definitions, concentration grids, paired seed design, primary concentration-response data, state-specific 4PL results, or final OAT sensitivity data are changed by the v1.0.0 release preparation.

## Before tagging v1.0.0

1. Apply the release-preparation patch.
2. Delete the obsolete root files listed in `DELETE_THESE_FILES.txt` if they are still present.
3. Commit and push to `main`.
4. Confirm GitHub Actions is green for Python 3.10, 3.11, and 3.12.
5. Record the final commit hash.
6. Create tag/release `v1.0.0` from that exact commit.
7. Archive the release in Zenodo (or another long-term repository) and obtain a DOI.
8. Update the manuscript Code and Data Availability section to cite `v1.0.0`, the final commit hash, and the DOI.
