# Manuscript reproduction workflow

This repository separates **fast re-analysis of version-controlled reference data** from
**full computational regeneration**.

## 1. Fast reproducibility check

From the repository root:

```bash
pip install -e ".[dev]"
python scripts/reproduce_paper.py
```

This copies the committed `n=20` replicate dataset into `results/reproduced/`, refits the
4-parameter logistic curves, repeats the paired-seed bootstrap, and copies the reference
OAT sensitivity table. It is intended for rapid verification of the manuscript analysis.

## 2. Regenerate the n=20 concentration-response simulations

```bash
python scripts/reproduce_paper.py --simulate
```

The simulation uses:

- 20 independent replicates per concentration,
- base seed `20261200` (`seed = base + replicate`),
- rate-preserving jittered input trains,
- jitter `0.5 ms`,
- identical input trains across concentrations within a replicate (paired design),
- the receptor-calibrated CA1-like reference configuration.

This is substantially slower than re-analysis of committed data.

## 3. Regenerate sensitivity analysis

```bash
python scripts/reproduce_paper.py --simulate --sensitivity
```

One-at-a-time perturbations are evaluated at each drug's fitted functional EC50:

- AMPA/EPSP amplitude: ±10%,
- GABA-A/IPSP magnitude: ±10%,
- explicit NMDA PSP scale: ±20%,
- NMDA activation threshold: ±2 mV,
- spike threshold: ±2 mV.

For every scenario, the drug-treated neuron is compared with a matched drug-free control
using the same stochastic input seed. If the matched drug-free control is silent, relative
suppression is undefined and is reported as not estimable rather than imputed.

## 4. Bootstrap uncertainty

`bootstrap_4pl.py` resamples replicate identities with replacement, preserving the paired
concentration structure within each resampled replicate. The default is 200 bootstrap
resamples, matching the current manuscript analysis.

## Interpretation

The reported functional EC50 values are properties of the calibrated computational
system and stimulation protocol. They are not molecular binding constants, clinical doses,
or predicted human brain exposures.
