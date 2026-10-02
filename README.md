# Pyramidal Receptor Simulator

[![Tests](https://github.com/dswietlik-cmd/pyramidal-receptor-sim/actions/workflows/tests.yml/badge.svg)](https://github.com/dswietlik-cmd/pyramidal-receptor-sim/actions/workflows/tests.yml)

**Pyramidal Receptor Simulator** is a reproducible Python model of a CA1-like pyramidal neuron for comparative *in silico* pharmacodynamics at **NMDA, AMPA, and GABA-A receptors**. The receptor is the pharmacological target; neuronal firing and synaptic responses are the functional readouts.

> **Research software.** This model is intended for computational research and method development. It is not a pharmacokinetic model, a clinical dose calculator, or a medical decision tool.

## Scientific question

The model is designed to compare how three receptor-level interventions alter neuronal excitability under a common computational framework:

| Drug | Target | Implemented action |
|---|---|---|
| Memantine | NMDA | concentration-dependent antagonism |
| Perampanel | AMPA | concentration-dependent antagonism |
| Diazepam | GABA-A | bounded positive allosteric modulation |

The computational chain is:

```text
drug concentration
       |
       v
receptor pharmacodynamics
       |
       v
AMPA / NMDA / GABA-A gain
       |
       v
synaptic response
       |
       v
somatic potential -> spikes -> firing-rate response
```

Only the current refactored implementation is distributed in this repository. Historical source code is not included.


## Healthy vs pathological pharmacology (v0.9.5)

The repository now includes a paired n=20 comparison of the same concentration-response protocols in the **healthy** state and the calibrated **moderate hyperexcitability** state (1.65x excitatory drive). The analysis reports firing, receptor/synaptic readouts, Pathological Activity Normalization (PAN), and an exploratory model selectivity index. See [`docs/HEALTHY_PATHOLOGICAL_DRUGS.md`](docs/HEALTHY_PATHOLOGICAL_DRUGS.md).

Run the complete paired experiment with:

```bash
pip install -e ".[fast]"
python scripts/run_healthy_pathological_drugs_fast.py
```

## Model at a glance

The reference CA1-like configuration contains 13 excitatory and 3 inhibitory inputs, a 0.5 ms integration step, separate AMPA/NMDA/GABA-A components, voltage-dependent NMDA gating, an explicit calibrated NMDA contribution to the somatic postsynaptic potential, and an NMDA-coupled plasticity state. Publication simulations use rate-preserving jittered input trains and recorded random seeds.

Outside the refractory period, the somatic postsynaptic potential is represented as

```text
V_post = V_rest + PSP_AMPA + PSP_NMDA + PSP_GABAA
```

Drug concentration is converted to receptor modulation before the neuronal simulation. Antagonist block uses a Hill-type mapping:

```text
block(C) = C^n / (IC50^n + C^n)
receptor_gain = 1 - block(C)
```

Diazepam uses a bounded Hill-type positive modulation of GABA-A. Full equations and assumptions are documented in [`docs/MODEL.md`](docs/MODEL.md), [`docs/PHARMACOLOGY.md`](docs/PHARMACOLOGY.md), and [`docs/CALIBRATION.md`](docs/CALIBRATION.md).

## Installation

Python 3.10 or newer is required.

```bash
git clone https://github.com/dswietlik-cmd/pyramidal-receptor-sim.git
cd pyramidal-receptor-sim
python -m venv .venv
source .venv/bin/activate       # Windows: .venv\Scripts\activate
pip install -e .
```

For development and automated tests:

```bash
pip install -e ".[dev]"
pytest
```

## Quick start

Run a single concentration:

```bash
python scripts/run_single.py --drug memantine --concentration 150 --unit nM
python scripts/run_single.py --drug perampanel --concentration 12 --unit nM
python scripts/run_single.py --drug diazepam --concentration 10 --unit nM
```

The output contains time-resolved spike state, somatic postsynaptic potential, receptor gains, and AMPA, NMDA, and GABA-A contributions.

## Reproduce the paper workflow

The canonical manuscript workflow uses **20 paired computational stochastic realizations per state and concentration**. Seed IDs are preserved across concentrations so that drug-induced changes are compared against matched stochastic input realizations.

Install the publication dependencies and run the final analysis with:

```bash
pip install -e ".[publication]"
python scripts/reproduce_paper.py
```

For full regeneration of the paired simulations and final OAT sensitivity analysis:

```bash
python scripts/reproduce_paper.py --simulate --sensitivity
```

The final main analysis uses **500 paired-seed bootstrap resamples** with bootstrap RNG seed `20261001`. Reference replicate-level inputs are stored under [`results/reference/`](results/reference/), while the canonical manuscript tables and figures are stored under [`results/publication_final/`](results/publication_final/).

Key paired input files are:

```text
healthy_two_pathologies_drug_replicates.csv
healthy_two_pathologies_paired_metrics.csv
```

Key manuscript outputs include:

```text
Table1_baseline_states.csv
Table2_state_specific_4PL_fits.csv
Table3_functional_target_concentrations.csv
Table4_therapeutic_window_grid.csv
Table5_best_effect_with_healthy_preservation_ge80.csv
Figure1_final_study_design.png ... Figure5_healthy_preserving_therapeutic_window.png
TableS3_sensitivity_oat_final_ec50.csv
FigureS1_sensitivity_oat_final_ec50.png
```

These files are reproducibility outputs, not immutable ground truth. Regenerate them after changing equations, parameters, drug mappings, pathology definitions, concentration grids, input statistics, or random-seed handling.

## Reference concentration-response results

For the final healthy-state analysis, the state-specific 4PL workflow produced the following **functional** EC50 estimates:

| Drug | Healthy-state functional EC50 | Bootstrap 95% CI |
|---|---:|---:|
| Perampanel | 12.70 nM | 11.37-13.99 nM |
| Memantine | 175.47 nM | 168.94-181.56 nM |
| Diazepam | 8.12 nM | 7.36-8.85 nM |

These values are emergent properties of this configured model. They are **not molecular binding constants, clinical target concentrations, or predictions of human brain exposure**. Functional EC50 is reported only when at least 50% mean suppression is observed within the tested concentration range.

## Sensitivity analysis

The publication workflow includes one-at-a-time sensitivity analyses for key excitability parameters. The current model is particularly sensitive to parameters that alter the operating point of the neuron, including NMDA activation threshold, spike threshold, and excitation/inhibition balance. A parameter set that silences the control neuron cannot yield an interpretable drug-induced firing-suppression estimate and should be reported as not estimable rather than forced into the analysis.

## Repository structure

```text
pyramidal-receptor-sim/
├── src/pyramidal_receptor_sim/
│   ├── model.py          # neuronal dynamics and receptor-specific PSP components
│   ├── drugs.py          # concentration -> receptor modulation
│   ├── inputs.py         # deterministic and jittered input generation
│   └── config.py         # reference model parameters
├── configs/              # reference configuration
├── scripts/              # simulation and 4PL fitting workflows
├── tests/                # automated tests
├── docs/                 # model, pharmacology, and calibration documentation
└── results/reference/    # reproducibility reference outputs
```

## Reproducibility checklist

For a manuscript analysis or tagged software release, record the repository commit or release tag, Python/package environment, configuration file, drug grid, number of replicates, base seed, jitter setting, simulation duration, and generated raw replicate table. Do not report only fitted EC50 values without retaining the underlying replicate-level output.

## Related publications

The computational framework builds on previous work on NMDA-receptor function, excitotoxicity, synaptic plasticity, hippocampal modeling, and *in silico* memantine therapy:

1. Swietlik D, Kusiak A, Krasny M, Bialowas J. *The Computer Simulation of Therapy with the NMDA Antagonist in Excitotoxic Neurodegeneration in an Alzheimer's Disease-like Pathology.* **J Clin Med.** 2022;11:1858. https://doi.org/10.3390/jcm11071858
2. Swietlik D, Kusiak A, Ossowska A. *Computational Modeling of Therapy with the NMDA Antagonist in Neurodegenerative Disease: Information Theory in the Mechanism of Action of Memantine.* **Int J Environ Res Public Health.** 2022;19:4727. https://doi.org/10.3390/ijerph19084727
3. Swietlik D, Bialowas J, Kusiak A, Krasny M. *Virtual Therapy with the NMDA Antagonist Memantine in Hippocampal Models of Moderate to Severe Alzheimer's Disease, In Silico Trials.* **Pharmaceuticals.** 2022;15:546. https://doi.org/10.3390/ph15050546

## Citation

A [`CITATION.cff`](CITATION.cff) file is included so GitHub can expose a **Cite this repository** action. The canonical repository is `https://github.com/dswietlik-cmd/pyramidal-receptor-sim`. For the pre-publication phase, cite the repository commit or the relevant v0.9.x tag when needed. The final manuscript-associated **v1.0.0** release will be archived with Zenodo (or another long-term repository), and its DOI will then be added here and to `CITATION.cff`.

## License

This project is released under the **BSD 3-Clause License**. See [`LICENSE`](LICENSE) for the full terms.

You may use, modify, and redistribute the software under those terms. Redistributions must retain the copyright notice and license conditions, and the name of the copyright holder may not be used to endorse or promote derived products without prior written permission.

## Development status

Current development release: **v0.9.9**. This is a pre-publication research version. The model, analysis workflow, parameterization, and reference results may still change before the manuscript-associated **v1.0.0** release.

## Manuscript reproduction pipeline

A single entry point is provided for the final pre-publication paired analysis:

```bash
pip install -e ".[publication]"
python scripts/reproduce_paper.py
```

For full regeneration of the paired concentration-response simulations and final OAT sensitivity analysis:

```bash
python scripts/reproduce_paper.py --simulate --sensitivity
```

See [`docs/PAPER_WORKFLOW.md`](docs/PAPER_WORKFLOW.md) for the exact paired-seed design,
500-resample bootstrap procedure, inferential statistical workflow, professional figure generation,
final sensitivity workflow, and interpretation limits. Statistical details are documented in
[`docs/STATISTICAL_ANALYSIS.md`](docs/STATISTICAL_ANALYSIS.md).

## Pathological hyperexcitability states

Version 0.9.2 adds a phenotype-level hyperexcitability calibration that changes only excitatory input drive, leaving AMPA, NMDA, GABA-A receptor gains and inhibitory input frequencies unchanged. The n=20 reference states are approximately 18 Hz (mild; 1.30× excitatory drive), 24 Hz (moderate; 1.65×), and 30 Hz (severe; 1.90×). See `docs/PATHOLOGY.md`.

## v0.9.4: matched pathological mechanisms

The repository now includes two moderate pathological phenotypes with similar mean firing
but different mechanisms: input-driven hyperexcitability and NMDA-driven excitotoxicity-like
activity. See `docs/NMDA_EXCITOTOXICITY.md` and `docs/TWO_PATHOLOGY_COMPARISON.md`.


## v0.9.5: GABA-A disinhibition boundary analysis

Version 0.9.5 adds a pathology-specific `gabaa_pathology_multiplier`, independent of diazepam pharmacological gain, and tests isolated GABA-A hypofunction across n=20 paired realizations. Under the current reference calibration, even complete modeled GABA-A loss increases mean firing only from about 10.97 to 11.28 Hz and therefore does **not** reproduce the pre-specified 18/24/30 Hz hyperexcitability phenotypes. The complete-loss condition is retained as a boundary/stress-test state rather than a matched moderate pathology. See `docs/GABAA_DISINHIBITION.md` and `docs/THREE_PATHOLOGY_COMPARISON.md`.


## v0.9.6: Final paired healthy-vs-pathology publication experiment

Version 0.9.6 freezes the current n=20 paired analysis comparing the healthy state with two matched moderate pathological states: input-driven hyperexcitability and NMDA-driven excitotoxicity-like pathology. It adds state-specific 4PL/paired-bootstrap analysis, a descriptive healthy-preserving therapeutic-window analysis, publication tables, Figures 1–5, and a Results draft in `results/publication_final/`. Functional EC50 values are reported only when 50% suppression is actually reached within the tested concentration range; extrapolated 4PL values are retained separately and are not treated as observed EC50 estimates.

## v0.9.7: Final sensitivity-aligned publication package

Version 0.9.7 adds the sensitivity analysis recalculated at the final healthy-state functional EC50 values used by the publication experiment, together with Supplementary Figure S1 and Table S3. It also removes an unnecessary pandas dependency from the publication-output pytest so the CI matrix remains lightweight. No primary healthy-versus-pathology concentration-response results were changed.

## v0.9.8: Reproducibility and manuscript-pipeline consistency

Version 0.9.8 aligns package metadata and documentation with the final pre-publication workflow. The single reproduction entry point now runs the final paired healthy/input-driven/NMDA-driven publication analysis with 500 paired-seed bootstrap resamples and routes sensitivity regeneration to the final EC50-aligned OAT script. The release also adds publication-analysis dependencies as an optional installation extra. No neuronal equations, pharmacodynamic parameters, pathology definitions, concentration grids, or manuscript numerical results are changed.
## v0.9.9: inferential statistics and publication-figure reproducibility

Version 0.9.9 adds the final paired inferential statistical analysis (Friedman tests, Holm-corrected paired Wilcoxon tests, Kendall's W, and matched rank-biserial effect sizes), version-controlled statistical Tables 9-11 and Supplementary Tables S4-S5, and a reproducible generator for professional publication figures. The reproduction entry point now regenerates the statistical outputs and journal-ready vector/PDF plus 800-dpi PNG figures. No neuronal equations, pharmacodynamic parameters, pathology definitions, concentration grids, primary concentration-response data, 4PL fits, or sensitivity data are changed.

