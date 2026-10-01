#!/usr/bin/env python3
"""Paired n=20 drug-response comparison across three pathology mechanisms.

Mechanisms:
1) input-driven moderate hyperexcitability (1.65x excitatory input drive),
2) NMDA-driven moderate excitotoxicity-like state (4.10x NMDA pathology drive),
3) maximal isolated GABA-A hypofunction (0x GABA-A pathology efficacy).

The third mechanism is intentionally retained as a boundary/stress-test state:
calibration showed that isolated GABA-A loss does NOT reach the pre-specified
18/24/30 Hz hyperexcitability targets in this reference model.
"""
from __future__ import annotations

import csv
import math
import random
import statistics
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from pyramidal_receptor_sim import load_config
from pyramidal_receptor_sim.drugs import receptor_gain_for_drug
from compare_two_pathologies_drugs_fast import simulate_fast, pars, GRIDS_NM

CFG = load_config(ROOT / "configs/ca1_reference.txt")
NE = int(CFG["N_EXCITATORY_INPUTS"])
NI = int(CFG["N_INHIBITORY_INPUTS"])
N = NE + NI
DT = float(CFG["DT_MS"])
STEPS = int(CFG["STEPS"])
BASE_FREQS = np.array([float(CFG[f"INPUT_{i}"]) for i in range(1, N + 1)], dtype=np.float64)
JITTER_MS = 0.5
N_REPLICATES = 20
BASE_SEED = 20261200

STATES = {
    "healthy": {"drive": 1.00, "nmda_path": 1.00, "gaba_path": 1.00, "matched_hyperexcitable": False},
    "input_hyperexcitable": {"drive": 1.65, "nmda_path": 1.00, "gaba_path": 1.00, "matched_hyperexcitable": True},
    "nmda_excitotoxic": {"drive": 1.00, "nmda_path": 4.10, "gaba_path": 1.00, "matched_hyperexcitable": True},
    "gabaa_complete_loss": {"drive": 1.00, "nmda_path": 1.00, "gaba_path": 0.00, "matched_hyperexcitable": False},
}


def make_stream(seed: int, drive: float):
    freqs = BASE_FREQS.copy()
    freqs[:NE] *= float(drive)
    rng = random.Random(seed)
    js = int(round(JITTER_MS / DT))
    ex = np.zeros((STEPS, NE), np.uint8)
    inh = np.zeros((STEPS, NI), np.uint8)
    for idx, freq in enumerate(freqs):
        if freq <= 0:
            continue
        period = max(1, int(round(1000.0 / (freq * DT))))
        for base in range(0, STEPS, period):
            t = max(0, min(STEPS - 1, base + (rng.randint(-js, js) if js else 0)))
            if idx < NE:
                ex[t, idx] = 1
            else:
                inh[t, idx - NE] = 1
    return ex, inh


def gains_for(drug: str, concentration_nM: float):
    rec = receptor_gain_for_drug(drug, concentration_nM, "nM")
    ampa = nmda = gaba = 1.0
    if rec["target"] == "AMPA":
        ampa = rec["receptor_gain"]
    elif rec["target"] == "NMDA":
        nmda = rec["receptor_gain"]
    elif rec["target"] == "GABAA":
        gaba = rec["receptor_gain"]
    return ampa, nmda, gaba, rec


def write_csv(path: Path, rows: list[dict]):
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)


def main():
    p = pars()
    streams = {
        (state, rep): make_stream(BASE_SEED + rep, spec["drive"])
        for state, spec in STATES.items()
        for rep in range(N_REPLICATES)
    }

    baseline_by_seed: dict[int, dict[str, float]] = {}
    for rep in range(N_REPLICATES):
        seed = BASE_SEED + rep
        baseline_by_seed[seed] = {}
        for state, spec in STATES.items():
            ex, inh = streams[(state, rep)]
            out = simulate_fast(ex, inh, *p, 1.0, 1.0, spec["gaba_path"], spec["nmda_path"])
            baseline_by_seed[seed][state] = out[0]

    rows: list[dict] = []
    for state, spec in STATES.items():
        for rep in range(N_REPLICATES):
            seed = BASE_SEED + rep
            ex, inh = streams[(state, rep)]
            baseline_hz = baseline_by_seed[seed][state]
            for drug, grid in GRIDS_NM.items():
                for concentration in grid:
                    ampa, nmda, gaba_drug, rec = gains_for(drug, float(concentration))
                    # Disease-level GABA efficacy and diazepam gain are distinct;
                    # their functional effects combine multiplicatively.
                    effective_gaba = spec["gaba_path"] * gaba_drug
                    out = simulate_fast(ex, inh, *p, ampa, nmda, effective_gaba, spec["nmda_path"])
                    rows.append({
                        "state": state,
                        "matched_hyperexcitable": spec["matched_hyperexcitable"],
                        "excitatory_drive_multiplier": spec["drive"],
                        "nmda_pathology_multiplier": spec["nmda_path"],
                        "gabaa_pathology_multiplier": spec["gaba_path"],
                        "drug": drug,
                        "replicate": rep + 1,
                        "seed": seed,
                        "concentration_nM": concentration,
                        "receptor_gain": rec["receptor_gain"],
                        "effective_gabaa_gain": effective_gaba,
                        "target_effect": rec["target_effect"],
                        "firing_hz": out[0],
                        "mean_isi_ms": out[1],
                        "mean_ampa_psp_mv": out[2],
                        "mean_nmda_psp_mv": out[3],
                        "mean_gaba_psp_mv": out[4],
                        "nmda_open_fraction": out[5],
                        "mean_plasticity_state": out[6],
                        "state_baseline_hz": baseline_hz,
                        "within_state_suppression_fraction": 1.0 - out[0] / baseline_hz if baseline_hz else float("nan"),
                    })

    ref = ROOT / "results" / "reference"
    ref.mkdir(parents=True, exist_ok=True)
    write_csv(ref / "healthy_three_pathologies_drug_replicates.csv", rows)

    # Summary table.
    summary_rows: list[dict] = []
    for state in STATES:
        for drug, grid in GRIDS_NM.items():
            for concentration in grid:
                grp = [r for r in rows if r["state"] == state and r["drug"] == drug and r["concentration_nM"] == concentration]
                s = {"state": state, "drug": drug, "concentration_nM": concentration, "n": len(grp)}
                for key in ["firing_hz", "mean_isi_ms", "mean_ampa_psp_mv", "mean_nmda_psp_mv", "mean_gaba_psp_mv", "nmda_open_fraction", "mean_plasticity_state", "within_state_suppression_fraction"]:
                    vals = [float(r[key]) for r in grp if not math.isnan(float(r[key]))]
                    s[key + "_mean"] = statistics.mean(vals) if vals else float("nan")
                    s[key + "_sd"] = statistics.stdev(vals) if len(vals) > 1 else 0.0
                summary_rows.append(s)
    write_csv(ref / "healthy_three_pathologies_drug_summary.csv", summary_rows)

    # Baseline mechanism comparison.
    base_rows = []
    for state, spec in STATES.items():
        vals = [baseline_by_seed[BASE_SEED + rep][state] for rep in range(N_REPLICATES)]
        healthy_vals = [baseline_by_seed[BASE_SEED + rep]["healthy"] for rep in range(N_REPLICATES)]
        paired_delta = [a - b for a, b in zip(vals, healthy_vals)]
        base_rows.append({
            "state": state,
            "matched_hyperexcitable": spec["matched_hyperexcitable"],
            "n": N_REPLICATES,
            "mean_firing_hz": statistics.mean(vals),
            "sd_firing_hz": statistics.stdev(vals),
            "mean_paired_delta_vs_healthy_hz": statistics.mean(paired_delta),
            "sd_paired_delta_vs_healthy_hz": statistics.stdev(paired_delta),
        })
    write_csv(ref / "three_pathology_baseline_comparison.csv", base_rows)

    # PAN/TSI only when the pathology is actually separated from healthy enough
    # to make normalization meaningful. GABA complete loss is retained as a
    # boundary condition and therefore receives NA PAN/TSI when delta < 2 Hz.
    metrics: list[dict] = []
    for state in ["input_hyperexcitable", "nmda_excitotoxic", "gabaa_complete_loss"]:
        for drug, grid in GRIDS_NM.items():
            for concentration in grid:
                h_by_seed = {r["seed"]: r for r in rows if r["state"] == "healthy" and r["drug"] == drug and r["concentration_nM"] == concentration}
                p_by_seed = {r["seed"]: r for r in rows if r["state"] == state and r["drug"] == drug and r["concentration_nM"] == concentration}
                for seed in sorted(h_by_seed):
                    h0 = baseline_by_seed[seed]["healthy"]
                    p0 = baseline_by_seed[seed][state]
                    hd = h_by_seed[seed]["firing_hz"]
                    pd = p_by_seed[seed]["firing_hz"]
                    denom = abs(p0 - h0)
                    estimable = denom >= 2.0
                    pan = 1.0 - abs(pd - h0) / denom if estimable else float("nan")
                    path_supp = (p0 - pd) / p0 if p0 > 0 else float("nan")
                    healthy_supp = (h0 - hd) / h0 if h0 > 0 else float("nan")
                    tsi = path_supp / (max(healthy_supp, 0.0) + 0.01) if estimable and not math.isnan(path_supp) else float("nan")
                    metrics.append({
                        "pathology_state": state,
                        "normalization_estimable": estimable,
                        "drug": drug,
                        "concentration_nM": concentration,
                        "seed": seed,
                        "healthy_baseline_hz": h0,
                        "healthy_drug_hz": hd,
                        "pathological_baseline_hz": p0,
                        "pathological_drug_hz": pd,
                        "baseline_delta_hz": p0 - h0,
                        "pathological_activity_normalization": pan,
                        "pathological_suppression_fraction": path_supp,
                        "healthy_suppression_fraction": healthy_supp,
                        "therapeutic_selectivity_index_exploratory": tsi,
                    })
    write_csv(ref / "healthy_three_pathologies_paired_metrics.csv", metrics)

    print("state, mean_hz, sd_hz, mean_delta_vs_healthy")
    for r in base_rows:
        print(f"{r['state']},{r['mean_firing_hz']:.3f},{r['sd_firing_hz']:.3f},{r['mean_paired_delta_vs_healthy_hz']:.3f}")


if __name__ == "__main__":
    main()
