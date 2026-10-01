#!/usr/bin/env python3
"""Run one drug concentration and save a time-series CSV."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from pyramidal_receptor_sim import Neuron, configure_neuron, load_config
from pyramidal_receptor_sim.inputs import deterministic_inputs


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--drug", required=True, choices=["memantine", "perampanel", "diazepam"])
    p.add_argument("--concentration", required=True, type=float)
    p.add_argument("--unit", default="nM", choices=["M", "mM", "uM", "nM", "pM"])
    p.add_argument("--config", default="configs/ca1_reference.txt")
    p.add_argument("--output", default="drug_run.csv")
    args = p.parse_args()

    cfg = load_config(args.config)
    neuron = configure_neuron(Neuron(), cfg)
    state = neuron.apply_drug_concentration(args.drug, args.concentration, args.unit)
    dt_ms = float(cfg["DT_MS"])
    steps = int(cfg["STEPS"])
    ne, ni = int(cfg["N_EXCITATORY_INPUTS"]), int(cfg["N_INHIBITORY_INPUTS"])
    freqs = [float(cfg[f"INPUT_{i}"]) for i in range(1, ne + ni + 1)]

    spikes = 0
    out = Path(args.output)
    with out.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["step", "time_ms", "spike", "membrane_potential_mV", "ampa_gain", "nmda_gain", "gabaa_gain", "ampa_psp_mV", "nmda_psp_mV", "gabaa_psp_mV"])
        for step in range(steps):
            ex, inh = deterministic_inputs(step, freqs, ne, ni, dt_ms)
            spike = neuron.step(ex, inh)
            spikes += int(spike)
            w.writerow([step, step * dt_ms, spike, neuron.membrane_potential, neuron.ampa_gain, neuron.nmda_gain, neuron.gabaa_gain, neuron.last_ampa_psp, neuron.last_nmda_psp, neuron.last_gaba_psp])

    duration_s = steps * dt_ms / 1000.0
    print(state)
    print(f"spikes={spikes}")
    print(f"firing_rate_hz={spikes / duration_s:.6f}")
    print(f"output={out}")


if __name__ == "__main__":
    main()
