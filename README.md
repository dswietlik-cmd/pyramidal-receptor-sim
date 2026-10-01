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


## Healthy vs pathological pharmacology (v0.9.3)

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

The reference concentration-response workflow uses **20 paired stochastic realizations per concentration**. Within a replicate, the same jittered input realization is reused across concentrations so that receptor modulation is compared against an identical input train.

```bash
python scripts/run_publication_n20.py
python scripts/fit_4pl.py publication_n20_replicates.csv --output curve_fits.csv
```

Reference outputs are included under [`results/reference/`](results/reference/):

```text
publication_n20_replicates.csv
publication_n20_summary.csv
publication_n20_fits.csv
publication_sensitivity_oat_n20.csv
```

These files are reference validation outputs, not immutable ground truth. Regenerate them after changing equations, parameters, drug mappings, input statistics, or random-seed handling.

## Reference concentration-response results

For the current reference configuration, four-parameter logistic fits to the n=20 workflow produced the following **functional** EC50 estimates:

| Drug | Functional EC50 | Interpretation |
|---|---:|---|
| Perampanel | 11.75 nM | AMPA-mediated functional firing suppression |
| Memantine | 176.98 nM | NMDA-mediated functional firing suppression |
| Diazepam | 8.11 nM | GABA-A-mediated functional firing suppression |

These values are emergent properties of this configured model. They are **not molecular binding constants, clinical target concentrations, or predictions of human brain exposure**. They depend on receptor scaling, input timing, spike threshold, NMDA gating, and excitation/inhibition balance.

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

A [`CITATION.cff`](CITATION.cff) file is included so GitHub can expose a **Cite this repository** action. The canonical repository is `https://github.com/dswietlik-cmd/pyramidal-receptor-sim`. For the pre-publication phase, cite the repository commit or v0.9.1 tag when needed. The final manuscript-associated **v1.0.0** release will be archived with Zenodo (or another long-term repository), and its DOI will then be added here and to `CITATION.cff`.

## License

This project is released under the **BSD 3-Clause License**. See [`LICENSE`](LICENSE) for the full terms.

You may use, modify, and redistribute the software under those terms. Redistributions must retain the copyright notice and license conditions, and the name of the copyright holder may not be used to endorse or promote derived products without prior written permission.

## Development status

Current development release: **v0.9.1**. This is a pre-publication research version. The model, analysis workflow, parameterization, and reference results may still change before the manuscript-associated **v1.0.0** release.

## Manuscript reproduction pipeline

A single entry point is provided for the current pre-publication analysis:

```bash
python scripts/reproduce_paper.py
```

For full regeneration of the concentration-response simulations and OAT sensitivity analysis:

```bash
python scripts/reproduce_paper.py --simulate --sensitivity
```

See [`docs/PAPER_WORKFLOW.md`](docs/PAPER_WORKFLOW.md) for the exact paired-seed design,
bootstrap procedure, perturbation definitions, and interpretation limits.

## Pathological hyperexcitability states

Version 0.9.2 adds a phenotype-level hyperexcitability calibration that changes only excitatory input drive, leaving AMPA, NMDA, GABA-A receptor gains and inhibitory input frequencies unchanged. The n=20 reference states are approximately 18 Hz (mild; 1.30× excitatory drive), 24 Hz (moderate; 1.65×), and 30 Hz (severe; 1.90×). See `docs/PATHOLOGY.md`.
