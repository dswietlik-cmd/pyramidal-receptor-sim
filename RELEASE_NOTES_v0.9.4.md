# v0.9.4 — matched-mechanism pathology comparison

- Added an explicit disease-level `nmda_pathology_multiplier`, separate from pharmacological `nmda_gain`.
- Added calibrated NMDA-driven excitotoxicity-like states (healthy, mild, moderate, severe).
- Added n=20 calibration outputs for the NMDA-driven phenotype.
- Added a paired n=20 comparison of healthy, input-driven moderate hyperexcitability, and NMDA-driven moderate excitotoxicity-like states across perampanel, memantine, and diazepam concentration grids.
- Added PAN and exploratory selectivity outputs for both pathologies.
- Added documentation and a therapeutic-window comparison table.
- Expanded tests for pathology/drug parameter separation.
