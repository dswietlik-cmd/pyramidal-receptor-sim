#!/usr/bin/env python3
"""Reproduce the n=20 concentration-response experiment using the reference Python model.

This script prioritizes transparency and reproducibility over speed. It uses identical
jittered input seeds across concentrations within each replicate (paired design).
"""

from __future__ import annotations

import csv
import statistics
from pathlib import Path

from pyramidal_receptor_sim import Neuron, configure_neuron, load_config
from pyramidal_receptor_sim.inputs import jittered_stream

GRIDS_NM = {
    "perampanel": [0, 10, 12, 14, 16, 18],
    "memantine": [0, 100, 125, 150, 175, 200, 250],
    "diazepam": [0, 5, 7.5, 10, 12.5, 15, 20, 30],
}
N_REPLICATES = 20
BASE_SEED = 20261200
JITTER_MS = 0.5


def simulate(neuron, ex, inh, dt_ms):
    spikes = 0
    spike_times = []
    for t in range(ex.shape[0]):
        if neuron.step(ex[t].tolist(), inh[t].tolist()):
            spikes += 1
            spike_times.append(t * dt_ms)
    duration_s = ex.shape[0] * dt_ms / 1000.0
    isi = [b - a for a, b in zip(spike_times, spike_times[1:])]
    return spikes / duration_s, (statistics.mean(isi) if isi else float("nan"))


def main():
    cfg = load_config("configs/ca1_reference.txt")
    dt = float(cfg["DT_MS"]); steps = int(cfg["STEPS"])
    ne = int(cfg["N_EXCITATORY_INPUTS"]); ni = int(cfg["N_INHIBITORY_INPUTS"])
    freqs = [float(cfg[f"INPUT_{i}"]) for i in range(1, ne + ni + 1)]
    rows = []

    for rep in range(N_REPLICATES):
        seed = BASE_SEED + rep
        ex, inh = jittered_stream(freqs, steps, dt, ne, ni, seed, JITTER_MS)
        for drug, grid in GRIDS_NM.items():
            ctl = configure_neuron(Neuron(), cfg)
            ctl.apply_drug_concentration(drug, 0, "nM")
            control_hz, _ = simulate(ctl, ex, inh, dt)
            for concentration in grid:
                n = configure_neuron(Neuron(), cfg)
                n.apply_drug_concentration(drug, concentration, "nM")
                hz, mean_isi = simulate(n, ex, inh, dt)
                rows.append({
                    "drug": drug, "replicate": rep + 1, "seed": seed,
                    "concentration_nM": concentration, "firing_hz": hz,
                    "mean_isi_ms": mean_isi, "control_hz": control_hz,
                    "suppression_fraction": 1.0 - hz / control_hz if control_hz else float("nan"),
                })

    out = Path("publication_n20_replicates.csv")
    with out.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
    print(f"wrote {out} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
