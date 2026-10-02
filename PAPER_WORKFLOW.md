# Manuscript reproduction workflow

This repository separates **re-analysis of version-controlled paired reference data** from
**full computational regeneration**. The canonical manuscript workflow is the paired
healthy-versus-two-pathology experiment used for the final publication analysis.

## 1. Environment

Use Python 3.10 or newer. From the repository root install the package with the publication
and development extras:

```bash
pip install -e ".[publication,dev]"
```

The `publication` extra installs the additional analysis dependencies used by the manuscript
workflow (`pandas`, `matplotlib`, and `numba`).

## 2. Fast manuscript re-analysis

From the repository root:

```bash
python scripts/reproduce_paper.py
```

This re-analyzes the committed paired n=20 reference datasets in `results/reference/` and
regenerates the canonical main-manuscript outputs in `results/publication_final/`:

- baseline-state summary,
- state-specific 4-parameter logistic fits,
- paired-seed bootstrap uncertainty,
- functional target concentrations,
- healthy-preserving therapeutic-window tables,
- Figures 1-5,
- the final Results draft.

The final analysis uses **500 paired-seed bootstrap resamples** and bootstrap RNG seed
`20261001`.

## 3. Regenerate the n=20 concentration-response simulations

```bash
python scripts/reproduce_paper.py --simulate
```

The paired simulation uses:

- 20 computational stochastic replicates per state and concentration,
- seeds `20261200` through `20261219`,
- a 0.5 ms integration step and 10 s simulation duration,
- rate-preserving input trains with +/-0.5 ms jitter,
- identical seed IDs across drug concentrations,
- healthy state: 1.00x excitatory drive and 1.00x NMDA pathology multiplier,
- input-driven pathology: 1.65x excitatory drive and 1.00x NMDA pathology multiplier,
- NMDA-driven pathology: 1.00x excitatory drive and 4.10x NMDA pathology multiplier.

The concentration grids are defined in `scripts/compare_two_pathologies_drugs_fast.py`.
Full regeneration rewrites the paired reference CSV files in `results/reference/` before
running the final analysis.

## 4. Regenerate the final sensitivity analysis

```bash
python scripts/reproduce_paper.py --sensitivity
```

For a complete regeneration from the model:

```bash
python scripts/reproduce_paper.py --simulate --sensitivity
```

The final OAT analysis uses the **final healthy-state functional EC50 values read directly
from `results/publication_final/Table2_state_specific_4PL_fits.csv`**. It evaluates:

- AMPA EPSP amplitude: +/-10%,
- GABA-A IPSP magnitude: +/-10%,
- explicit NMDA PSP scale: +/-20%,
- NMDA activation threshold: +/-2 mV,
- spike threshold: +/-2 mV.

For each perturbation, the drug-treated simulation is compared with a matched drug-free
control using the same stochastic input seed. If the matched control is silent, relative
suppression is not estimable rather than imputed. The outputs are:

- `results/publication_final/TableS3_sensitivity_oat_final_ec50.csv`,
- `results/publication_final/FigureS1_sensitivity_oat_final_ec50.png`.

## 5. Interpretation limits

Functional EC50 values are properties of the calibrated drug-receptor-synapse-neuron
system and simulation protocol. They are not molecular binding constants, clinical doses,
or predicted human brain exposures. The >=80% healthy-firing criterion is a descriptive
model threshold, not a clinically validated therapeutic index.
