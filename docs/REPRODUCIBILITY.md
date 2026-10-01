# Reproducing the reference analysis

## Environment

Use Python 3.10 or newer and install the package from the repository root with `pip install -e .`. For tests, install the `dev` extras and run `pytest`.

## Reference workflow

1. Keep `configs/ca1_reference.txt` unchanged.
2. Run `python scripts/run_publication_n20.py`.
3. Preserve the replicate-level CSV, including seeds and concentrations.
4. Run `python scripts/fit_4pl.py publication_n20_replicates.csv --output curve_fits.csv`.
5. Compare regenerated summaries and fits with `results/reference/`.

## Paired stochastic design

The publication workflow uses the same jittered input realization across concentrations within a replicate. This pairing is intentional: changing the drug concentration should not simultaneously change the stochastic input realization used for that paired comparison.

## Changes requiring regeneration

Regenerate all reference outputs after modifying neuronal equations, configuration values, drug concentration-response parameters, input generation, simulation duration, random-seed logic, or analysis code.

## Reporting

A reproducible analysis should report the software release/commit, configuration, drug concentration grid, replicate count, seed strategy, simulation duration, jitter parameter, fitting model, and replicate-level output. Functional EC50 values should be described as model-level functional endpoints rather than molecular binding constants.
