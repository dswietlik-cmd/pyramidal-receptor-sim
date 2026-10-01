#!/usr/bin/env python3
"""Run one-at-a-time (OAT) sensitivity analysis at each drug's functional EC50.

The analysis uses the same 20 paired stochastic input seeds for drug-free and drug
conditions within every perturbation scenario. Outputs are suppression fractions relative
to a matched, perturbed, drug-free control.
"""
from __future__ import annotations

import argparse
import csv
import math
import statistics
from pathlib import Path

from pyramidal_receptor_sim import Neuron, configure_neuron, load_config
from pyramidal_receptor_sim.inputs import jittered_stream

N_REPLICATES = 20
BASE_SEED = 20261200
JITTER_MS = 0.5

PERTURBATIONS = {
    "baseline": lambda n: None,
    "EPSPd_-10%": lambda n: setattr(n, "ampa_epsp_amplitude", n.ampa_epsp_amplitude * 0.90),
    "EPSPd_+10%": lambda n: setattr(n, "ampa_epsp_amplitude", n.ampa_epsp_amplitude * 1.10),
    "IPSP_mag_-10%": lambda n: setattr(n, "gabaa_ipsp_amplitude", n.gabaa_ipsp_amplitude * 0.90),
    "IPSP_mag_+10%": lambda n: setattr(n, "gabaa_ipsp_amplitude", n.gabaa_ipsp_amplitude * 1.10),
    "NMDA_scale_-20%": lambda n: setattr(n, "nmda_psp_scale", n.nmda_psp_scale * 0.80),
    "NMDA_scale_+20%": lambda n: setattr(n, "nmda_psp_scale", n.nmda_psp_scale * 1.20),
    "CaMT_-2mV": lambda n: setattr(n, "nmda_activation_threshold", n.nmda_activation_threshold - 2.0),
    "CaMT_+2mV": lambda n: setattr(n, "nmda_activation_threshold", n.nmda_activation_threshold + 2.0),
    "Threshold_-2mV": lambda n: setattr(n, "spike_threshold", n.spike_threshold - 2.0),
    "Threshold_+2mV": lambda n: setattr(n, "spike_threshold", n.spike_threshold + 2.0),
}


def simulate(neuron: Neuron, ex, inh, dt_ms: float) -> float:
    spikes = 0
    for t in range(ex.shape[0]):
        spikes += int(neuron.step(ex[t].tolist(), inh[t].tolist()))
    return spikes / (ex.shape[0] * dt_ms / 1000.0)


def load_ec50s(path: Path) -> dict[str, float]:
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    return {row["drug"]: float(row["EC50_nM"]) for row in rows}


def make_neuron(cfg, scenario: str) -> Neuron:
    n = configure_neuron(Neuron(), cfg)
    PERTURBATIONS[scenario](n)
    # Rebuild kernels/registers after any parameter perturbation while preserving the change.
    n.reset_state()
    return n


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--config", default="configs/ca1_reference.txt")
    p.add_argument("--fits", default="results/reference/publication_n20_fits.csv")
    p.add_argument("--output", default="publication_sensitivity_oat_n20.csv")
    args = p.parse_args()

    cfg = load_config(args.config)
    ec50s = load_ec50s(Path(args.fits))
    dt = float(cfg["DT_MS"]); steps = int(cfg["STEPS"])
    ne = int(cfg["N_EXCITATORY_INPUTS"]); ni = int(cfg["N_INHIBITORY_INPUTS"])
    freqs = [float(cfg[f"INPUT_{i}"]) for i in range(1, ne + ni + 1)]

    summary = []
    for drug, concentration in ec50s.items():
        for scenario in PERTURBATIONS:
            values = []
            for rep in range(N_REPLICATES):
                seed = BASE_SEED + rep
                ex, inh = jittered_stream(freqs, steps, dt, ne, ni, seed, JITTER_MS)

                ctl = make_neuron(cfg, scenario)
                control_hz = simulate(ctl, ex, inh, dt)

                treated = make_neuron(cfg, scenario)
                treated.apply_drug_concentration(drug, concentration, "nM")
                drug_hz = simulate(treated, ex, inh, dt)

                if control_hz > 0:
                    values.append(1.0 - drug_hz / control_hz)

            if values:
                mean = statistics.mean(values)
                sd = statistics.stdev(values) if len(values) > 1 else 0.0
                n_valid = len(values)
            else:
                mean = float("nan"); sd = float("nan"); n_valid = 0
            summary.append({
                "drug": drug,
                "scenario": scenario,
                "concentration_nM": concentration,
                "n": n_valid,
                "mean_suppression_fraction": mean,
                "sd_suppression_fraction": sd,
            })

    with Path(args.output).open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=summary[0].keys())
        w.writeheader(); w.writerows(summary)
    print(f"wrote {args.output} ({len(summary)} rows)")


if __name__ == "__main__":
    main()
