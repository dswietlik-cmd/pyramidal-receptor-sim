#!/usr/bin/env python3
"""Calibrate GABA-A hypofunction/disinhibition using a pathology-specific receptor multiplier.

The pathology multiplier is separate from the pharmacological ``gabaa_gain``.
During simulation their effects combine multiplicatively. This allows reduced
GABA-A efficacy to define disease state while diazepam remains a separate drug
perturbation.
"""
from __future__ import annotations

import argparse
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
from compare_two_pathologies_drugs_fast import simulate_fast, pars

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
TARGETS = {"mild": 18.0, "moderate": 24.0, "severe": 30.0}


def make_stream(seed: int):
    rng = random.Random(seed)
    js = int(round(JITTER_MS / DT))
    ex = np.zeros((STEPS, NE), np.uint8)
    inh = np.zeros((STEPS, NI), np.uint8)
    for idx, freq in enumerate(BASE_FREQS):
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


def write_csv(path: Path, rows: list[dict]):
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)


def summarize(rows: list[dict]) -> list[dict]:
    out = []
    factors = sorted({float(r["gabaa_pathology_multiplier"]) for r in rows}, reverse=True)
    metrics = [
        "firing_hz", "mean_isi_ms", "mean_ampa_psp_mv", "mean_nmda_psp_mv",
        "mean_gaba_psp_mv", "nmda_open_fraction", "mean_plasticity_state",
    ]
    for fac in factors:
        grp = [r for r in rows if float(r["gabaa_pathology_multiplier"]) == fac]
        s = {"gabaa_pathology_multiplier": fac, "n": len(grp)}
        for key in metrics:
            vals = [float(r[key]) for r in grp if not math.isnan(float(r[key]))]
            s[key + "_mean"] = statistics.mean(vals) if vals else float("nan")
            s[key + "_sd"] = statistics.stdev(vals) if len(vals) > 1 else 0.0
        out.append(s)
    return out


def closest_states(summary: list[dict]) -> list[dict]:
    selected = []
    # Healthy is explicitly the unmodified inhibitory drive.
    healthy = min(summary, key=lambda r: abs(float(r["gabaa_pathology_multiplier"]) - 1.0))
    selected.append({"state": "healthy", "target_firing_hz": 12.0, **healthy})
    for name, target in TARGETS.items():
        best = min(summary, key=lambda r: abs(float(r["firing_hz_mean"]) - target))
        selected.append({"state": name, "target_firing_hz": target, **best})
    return selected


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fine", action="store_true", help="Use a dense 0.01 GABA-A pathology grid.")
    args = ap.parse_args()
    if args.fine:
        factors = [round(x, 2) for x in np.arange(1.0, -0.001, -0.01)]
    else:
        factors = [round(x, 2) for x in np.arange(1.0, -0.001, -0.05)]

    p = pars()
    rows: list[dict] = []
    for rep in range(N_REPLICATES):
        seed = BASE_SEED + rep
        ex, inh = make_stream(seed)
        for fac in factors:
            out = simulate_fast(ex, inh, *p, 1.0, 1.0, fac, 1.0)
            rows.append({
                "gabaa_pathology_multiplier": fac,
                "replicate": rep + 1,
                "seed": seed,
                "firing_hz": out[0],
                "mean_isi_ms": out[1],
                "mean_ampa_psp_mv": out[2],
                "mean_nmda_psp_mv": out[3],
                "mean_gaba_psp_mv": out[4],
                "nmda_open_fraction": out[5],
                "mean_plasticity_state": out[6],
            })

    ref = ROOT / "results" / "reference"
    ref.mkdir(parents=True, exist_ok=True)
    suffix = "fine" if args.fine else "coarse"
    rep_path = ref / f"gabaa_hypofunction_calibration_{suffix}_replicates.csv"
    sum_path = ref / f"gabaa_hypofunction_calibration_{suffix}_summary.csv"
    sel_path = ref / f"gabaa_hypofunction_selected_states_{suffix}.csv"
    write_csv(rep_path, rows)
    summary = summarize(rows)
    write_csv(sum_path, summary)
    selected = closest_states(summary)
    write_csv(sel_path, selected)

    print("state,target_hz,gabaa_pathology_multiplier,mean_hz,sd_hz")
    for r in selected:
        print(f"{r['state']},{r['target_firing_hz']:.1f},{float(r['gabaa_pathology_multiplier']):.2f},{float(r['firing_hz_mean']):.3f},{float(r['firing_hz_sd']):.3f}")


if __name__ == "__main__":
    main()
