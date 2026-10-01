# Pharmacodynamic parameterization

Drug concentration-response parameters are kept in `src/pyramidal_receptor_sim/drugs.py` and are
separate from the neuronal equations.

## Memantine / NMDA

The reference mapping uses `IC50 = 0.79 uM` and Hill coefficient `n = 0.92` for NMDA-receptor
inhibition. The model applies the resulting fractional block to the NMDA functional gain.

Reference identifier: PubMed PMID 19371579.

## Perampanel / AMPA

The reference mapping uses `IC50 = 0.56 uM` and Hill coefficient `n = 0.80` for AMPA-receptor
inhibition. The model applies the resulting fractional block to the AMPA functional gain.

Reference identifiers used during parameterization include PubMed PMID 25229608 and related
recombinant-receptor electrophysiology summarized in the project calibration records.

## Diazepam / GABA-A

The reference mapping uses `EC50 = 0.025 uM` for positive allosteric modulation and a bounded
maximal response ratio of 6.1. The current release uses Hill coefficient `n = 1.0` as a declared
modeling assumption; this parameter should be included in sensitivity analysis when the GABA-A
mapping is central to an experiment.

Reference identifier: PMC3852889.

## Interpretation

These values define receptor-level pharmacodynamic mappings under specific experimental conditions.
They are not interchangeable with clinical dose, free brain concentration, plasma exposure, or a
molecular affinity measured in a different assay. The publication workflow therefore distinguishes
literature IC50/EC50 inputs from emergent functional EC50 values obtained from neuronal firing.
