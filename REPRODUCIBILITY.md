# Reproducing the final manuscript analysis

## Environment

Use Python 3.10 or newer. From the repository root:

```bash
pip install -e ".[publication,dev]"
pytest -q
```

The `publication` extra installs the analysis dependencies required by the final workflow.

## Canonical workflow

For re-analysis of the version-controlled paired n=20 reference datasets:

```bash
python scripts/reproduce_paper.py
```

For full regeneration of the paired concentration-response simulations and the final OAT
sensitivity analysis:

```bash
python scripts/reproduce_paper.py --simulate --sensitivity
```

The canonical main analysis is implemented in `scripts/final_publication_analysis.py` and
uses 500 paired-seed bootstrap resamples with bootstrap RNG seed `20261001`.

## Paired stochastic design

The final experiment uses seeds `20261200`-`20261219`. The same seed IDs are reused
across concentrations so that changes in receptor modulation are not confounded by a new
stochastic input realization. Input trains are rate preserving with +/-0.5 ms jitter.

## Final states

- Healthy: excitatory drive 1.00x; NMDA pathology multiplier 1.00x.
- Input-driven moderate hyperexcitability: excitatory drive 1.65x; NMDA pathology multiplier 1.00x.
- NMDA-driven excitotoxicity-like state: excitatory drive 1.00x; NMDA pathology multiplier 4.10x.

## Changes requiring regeneration

Regenerate reference and publication outputs after modifying neuronal equations,
configuration values, drug concentration-response parameters, concentration grids, input
generation, simulation duration, random-seed logic, pathology definitions, or analysis code.

## Reporting

A reproducible analysis should report the software release/commit, configuration, drug
concentration grids, replicate count, seed strategy, simulation duration, jitter parameter,
4PL model, bootstrap settings, and replicate-level outputs. Functional EC50 values should
be described as model-level functional endpoints rather than molecular binding constants.
