#!/usr/bin/env python3
"""Paired n=20 comparison of healthy, input-driven hyperexcitability, and NMDA-driven excitotoxicity-like states.

Runs concentration-response simulations for perampanel, memantine, and diazepam
in the calibrated healthy state (1.0x excitatory drive) and the calibrated moderate
hyperexcitability state (1.65x excitatory drive). Identical seed IDs are used across
states and concentrations. The fast Numba core is algebraically matched to the
reference model for this experiment.
"""
from __future__ import annotations

import csv
import math
import random
import statistics
import sys
from pathlib import Path

import numpy as np
from numba import njit

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pyramidal_receptor_sim import Neuron, configure_neuron, load_config
from pyramidal_receptor_sim.drugs import receptor_gain_for_drug
from pyramidal_receptor_sim.inputs import jittered_stream
from pyramidal_receptor_sim.pathology import HEALTHY, MODERATE, frequencies_for_state

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

GRIDS_NM = {
    "perampanel": [0, 6, 8, 10, 12, 14, 16, 18, 22],
    "memantine": [0, 50, 100, 125, 150, 175, 200, 250, 350, 500],
    "diazepam": [0, 2.5, 5, 7.5, 10, 12.5, 15, 20, 30, 60],
}
STATES = {
    "healthy": {"drive": 1.0, "nmda_path": 1.0},
    "input_hyperexcitable": {"drive": 1.65, "nmda_path": 1.0},
    "nmda_excitotoxic": {"drive": 1.0, "nmda_path": 4.10},
}


@njit(cache=True)
def influence(source, dest, n_total, lsw):
    if dest >= source:
        return 1.0 - (1.0 - lsw) / (n_total - 1.0) * (dest - source)
    if (source - dest) <= n_total / 4.0:
        return 1.0 - 4.0 / n_total * (source - dest)
    return 0.0


@njit(cache=True)
def simulate_fast(
    ex, inh, EPSPd, IPSPd, CaMT, NMDA_scale, ReP, RePB, Threshold, KEP,
    CRefrak, FQ, Clog, SleepC, Rsp, LSW, ampa_gain, nmda_gain, gaba_gain, nmda_pathology_multiplier
):
    ExTab = np.zeros(33); MexTab = np.zeros(33); InhTab = np.zeros(33)
    ExTab[0] = ExTab[1] = ExTab[2] = ExTab[32] = 0.0
    for i in range(1, 5): ExTab[2 + i] = EPSPd / 4.0 * i
    for i in range(1, 26): ExTab[6 + i] = EPSPd / 26.0 * (26 - i)
    MexTab[0] = MexTab[1] = MexTab[2] = MexTab[32] = 0.0
    for i in range(1, 5): MexTab[2 + i] = EPSPd / 20.0 * i
    for i in range(1, 26): MexTab[6 + i] = EPSPd / 99.0 * (26 - i)
    InhTab[0] = InhTab[1] = InhTab[2] = 0.0
    for i in range(1, 5): InhTab[2 + i] = IPSPd / 4.0 * i
    for i in range(1, 20): InhTab[6 + i] = IPSPd / 20.0 * (20 - i)

    STEx = np.full((NE, 33), ReP); STMex = np.full((NE, 33), ReP); STInh = np.full((NI, 33), ReP)
    MemVol = np.zeros(NE); SpinePot = np.full(N, ReP); CaPath = np.zeros(NE, np.int64)
    InhDel = np.zeros(NI, np.int64); SleepCout = np.zeros(NI, np.int64)
    Refrak = 0; PSP = ReP; SleepState = 0; spikes = 0; last_sp = -1.0; isi_sum = 0.0; isi_n = 0
    ampa_sum = 0.0; nmda_sum = 0.0; gaba_sum = 0.0; gate_sum = 0.0; plast_sum = 0.0

    for step in range(ex.shape[0]):
        for ii in range(NE):
            if ex[step, ii] != 0:
                adaptation = ((ReP - Rsp) - (ReP - SpinePot[ii])) / (ReP - Rsp)
                memory = 1.0 + (1.0 / 6.0) * math.log(MemVol[ii] + 1.0) / Clog
                factor = adaptation * memory
                for k in range(33): STEx[ii, k] += factor * ampa_gain * ExTab[k]
        for ii in range(NI):
            if inh[step, ii] != 0:
                InhDel[ii] = 6
                for k in range(33): STInh[ii, k] += gaba_gain * InhTab[k]
        for ii in range(NE):
            if ex[step, ii] != 0 and CaPath[ii] == 1:
                for k in range(33): STMex[ii, k] += nmda_pathology_multiplier * nmda_gain * MexTab[k]
            power = (STMex[ii, 0] - ReP) * 6.0
            MemVol[ii] += math.exp(power) - 1.0
            if MemVol[ii] > FQ: MemVol[ii] -= FQ
            else: MemVol[ii] = 0.0
            number = ii + 1; pot = ReP
            for m in range(NE): pot += (STEx[m, 0] - ReP) * influence(m + 1, number, N, LSW)
            for m in range(NI): pot += (STInh[m, 0] - ReP) * influence(m + 1 + NE, number, N, LSW)
            SpinePot[ii] = pot; CaPath[ii] = 1 if pot >= CaMT else 0
        for jj in range(NI):
            number = NE + jj + 1; pot = ReP
            for m in range(NE): pot += (STEx[m, 0] - ReP) * influence(m + 1, number, N, LSW)
            for m in range(NI): pot += (STInh[m, 0] - ReP) * influence(m + 1 + NE, number, N, LSW)
            SpinePot[NE + jj] = pot

        SOut = 0; ampa_comp = 0.0; nmda_comp = 0.0; gaba_comp = 0.0
        if Refrak == 0:
            PSP = ReP
            for ii in range(NE):
                number = ii + 1; weight = (1.0 - LSW) / (NE - 1.0) * (number - NE) + 1.0
                a = weight * (STEx[ii, 0] - ReP)
                nm = NMDA_scale * weight * (STMex[ii, 0] - ReP)
                ampa_comp += a; nmda_comp += nm; PSP += a + nm
            for ii in range(NI):
                g = STInh[ii, 0] - ReP; gaba_comp += g; PSP += g
            if PSP < KEP: PSP = KEP
            if PSP >= Threshold:
                SOut = 1; Refrak = CRefrak
                for ii in range(NE):
                    number = ii + 1; adrv = ReP - ((ReP - RePB) / NE) * (NE - number)
                    for k in range(33): STEx[ii, k] = adrv
                for ii in range(NI):
                    for k in range(33): STInh[ii, k] = ReP
            else:
                for ii in range(NE):
                    for k in range(32): STEx[ii, k] = STEx[ii, k + 1]
                    for k in range(32): STMex[ii, k] = STMex[ii, k + 1]
                if SleepState == 0:
                    for ii in range(NI):
                        for k in range(32): STInh[ii, k] = STInh[ii, k + 1]
                else:
                    for ii in range(NI):
                        if InhDel[ii] > 0:
                            InhDel[ii] -= 1
                            for k in range(32): STInh[ii, k] = STInh[ii, k + 1]
                        else:
                            SleepCout[ii] += 1
                            if SleepCout[ii] == SleepC:
                                SleepCout[ii] = 0
                                for k in range(32): STInh[ii, k] = STInh[ii, k + 1]
        else:
            for ii in range(NE):
                for k in range(32): STEx[ii, k] = STEx[ii, k + 1]
                for k in range(32): STMex[ii, k] = STMex[ii, k + 1]
            if SleepState == 0:
                for ii in range(NI):
                    for k in range(32): STInh[ii, k] = STInh[ii, k + 1]
            else:
                for ii in range(NI):
                    if InhDel[ii] > 0:
                        InhDel[ii] -= 1
                        for k in range(32): STInh[ii, k] = STInh[ii, k + 1]
                    else:
                        SleepCout[ii] += 1
                        if SleepCout[ii] == SleepC:
                            SleepCout[ii] = 0
                            for k in range(32): STInh[ii, k] = STInh[ii, k + 1]
            if Refrak == (CRefrak - 1): PSP = ReP
            Refrak -= 1
        SleepState = 1 if PSP < ReP else 0
        if SOut == 1:
            spikes += 1; tm = step * DT
            if last_sp >= 0: isi_sum += tm - last_sp; isi_n += 1
            last_sp = tm
        ampa_sum += ampa_comp; nmda_sum += nmda_comp; gaba_sum += gaba_comp
        gate_sum += np.sum(CaPath) / NE; plast_sum += np.mean(MemVol)

    duration_s = ex.shape[0] * DT / 1000.0
    return (
        spikes / duration_s,
        isi_sum / isi_n if isi_n > 0 else np.nan,
        ampa_sum / ex.shape[0], nmda_sum / ex.shape[0], gaba_sum / ex.shape[0],
        gate_sum / ex.shape[0], plast_sum / ex.shape[0]
    )


def pars():
    return (
        float(CFG["AMPA_EPSP_AMPLITUDE_MV"]), float(CFG["GABAA_IPSP_AMPLITUDE_MV"]),
        float(CFG["NMDA_ACTIVATION_THRESHOLD_MV"]), float(CFG["NMDA_PSP_SCALE"]),
        float(CFG["RESTING_POTENTIAL_MV"]), float(CFG["RESET_POTENTIAL_UPPER_MV"]),
        float(CFG["SPIKE_THRESHOLD_MV"]), float(CFG["MINIMUM_MEMBRANE_POTENTIAL_MV"]),
        int(CFG["REFRACTORY_STEPS"]), float(CFG["PLASTICITY_DECAY_STEP"]),
        float(CFG["PLASTICITY_LOG_SCALE"]), int(CFG["INHIBITORY_DECAY_INTERVAL"]),
        float(CFG["ADAPTATION_REFERENCE_POTENTIAL_MV"]), float(CFG["MINIMUM_SYNAPTIC_WEIGHT"]),
    )


def make_stream(seed: int, state_name: str):
    spec = STATES[state_name]
    freqs = BASE_FREQS.copy()
    freqs[:NE] *= float(spec["drive"])
    rng = random.Random(seed); js = int(round(JITTER_MS / DT))
    ex = np.zeros((STEPS, NE), np.uint8); inh = np.zeros((STEPS, NI), np.uint8)
    for idx, freq in enumerate(freqs):
        if freq <= 0: continue
        period = max(1, int(round(1000.0 / (freq * DT))))
        for base in range(0, STEPS, period):
            t = max(0, min(STEPS - 1, base + (rng.randint(-js, js) if js else 0)))
            if idx < NE: ex[t, idx] = 1
            else: inh[t, idx - NE] = 1
    return ex, inh


def gains_for(drug: str, concentration_nM: float):
    rec = receptor_gain_for_drug(drug, concentration_nM, "nM")
    ampa = nmda = gaba = 1.0
    if rec["target"] == "AMPA": ampa = rec["receptor_gain"]
    elif rec["target"] == "NMDA": nmda = rec["receptor_gain"]
    elif rec["target"] == "GABAA": gaba = rec["receptor_gain"]
    return ampa, nmda, gaba, rec


def reference_python_hz(seed: int, state_name: str, drug: str, concentration_nM: float):
    freqs = frequencies_for_state(BASE_FREQS, NE, STATES[state_name])
    ex, inh = jittered_stream(freqs, STEPS, DT, NE, NI, seed, JITTER_MS)
    neuron = configure_neuron(Neuron(), CFG)
    neuron.apply_drug_concentration(drug, concentration_nM, "nM")
    spikes = 0
    for t in range(STEPS): spikes += neuron.step(ex[t].tolist(), inh[t].tolist())
    return spikes / (STEPS * DT / 1000.0)


def write_csv(path: Path, rows: list[dict]):
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)


def main():
    p = pars()
    # Fast-core experiment. Pathology-specific NMDA drive is explicit in this script.
    # Cache one input train per state x seed; use it for all drugs/concentrations.
    streams = {(state, rep): make_stream(BASE_SEED + rep, state) for state in STATES for rep in range(N_REPLICATES)}

    rows: list[dict] = []
    for state_name, state in STATES.items():
        for rep in range(N_REPLICATES):
            seed = BASE_SEED + rep; ex, inh = streams[(state_name, rep)]
            baseline = simulate_fast(ex, inh, *p, 1.0, 1.0, 1.0, STATES[state_name]["nmda_path"])
            baseline_hz = baseline[0]
            for drug, grid in GRIDS_NM.items():
                for concentration in grid:
                    ampa, nmda, gaba, rec = gains_for(drug, float(concentration))
                    out = simulate_fast(ex, inh, *p, ampa, nmda, gaba, STATES[state_name]["nmda_path"])
                    rows.append({
                        "state": state_name,
                        "drive_multiplier": STATES[state_name]["drive"],
                        "nmda_pathology_multiplier": STATES[state_name]["nmda_path"],
                        "drug": drug,
                        "replicate": rep + 1,
                        "seed": seed,
                        "concentration_nM": concentration,
                        "receptor_gain": rec["receptor_gain"],
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

    ref = ROOT / "results" / "reference"; ref.mkdir(parents=True, exist_ok=True)
    rep_path = ref / "healthy_two_pathologies_drug_replicates.csv"; write_csv(rep_path, rows)

    # Paired metrics for each pathology versus the same-seed healthy baseline.
    baseline_by_seed = {}
    for rep in range(N_REPLICATES):
        seed = BASE_SEED + rep
        baseline_by_seed[seed] = {state: simulate_fast(*streams[(state, rep)], *p, 1.0, 1.0, 1.0, STATES[state]["nmda_path"])[0] for state in STATES}

    metrics_rows = []
    for pathology_state in ["input_hyperexcitable", "nmda_excitotoxic"]:
        for drug, grid in GRIDS_NM.items():
            for concentration in grid:
                hrows = [r for r in rows if r["state"] == "healthy" and r["drug"] == drug and r["concentration_nM"] == concentration]
                prows = [r for r in rows if r["state"] == pathology_state and r["drug"] == drug and r["concentration_nM"] == concentration]
                h_by_seed = {r["seed"]: r for r in hrows}; p_by_seed = {r["seed"]: r for r in prows}
                for seed in sorted(h_by_seed):
                    h0 = baseline_by_seed[seed]["healthy"]; p0 = baseline_by_seed[seed][pathology_state]
                    hd = h_by_seed[seed]["firing_hz"]; pd = p_by_seed[seed]["firing_hz"]
                    denom = abs(p0 - h0)
                    pan = 1.0 - abs(pd - h0) / denom if denom > 0 else float("nan")
                    path_supp = (p0 - pd) / p0 if p0 > 0 else float("nan")
                    healthy_supp = (h0 - hd) / h0 if h0 > 0 else float("nan")
                    tsi = path_supp / (max(healthy_supp, 0.0) + 0.01) if not math.isnan(path_supp) else float("nan")
                    metrics_rows.append({
                        "pathology_state": pathology_state, "drug": drug, "concentration_nM": concentration, "seed": seed,
                        "healthy_baseline_hz": h0, "healthy_drug_hz": hd,
                        "pathological_baseline_hz": p0, "pathological_drug_hz": pd,
                        "pathological_activity_normalization": pan,
                        "pathological_suppression_fraction": path_supp,
                        "healthy_suppression_fraction": healthy_supp,
                        "therapeutic_selectivity_index_exploratory": tsi,
                    })
    metric_path = ref / "healthy_two_pathologies_paired_metrics.csv"; write_csv(metric_path, metrics_rows)

    # Summaries across n=20.
    summary_rows = []
    for state_name in STATES:
        for drug, grid in GRIDS_NM.items():
            for concentration in grid:
                grp = [r for r in rows if r["state"] == state_name and r["drug"] == drug and r["concentration_nM"] == concentration]
                s = {"state": state_name, "drug": drug, "concentration_nM": concentration, "n": len(grp)}
                for key in ["firing_hz", "mean_isi_ms", "mean_ampa_psp_mv", "mean_nmda_psp_mv", "mean_gaba_psp_mv", "nmda_open_fraction", "mean_plasticity_state", "within_state_suppression_fraction"]:
                    vals = [float(r[key]) for r in grp if not math.isnan(float(r[key]))]
                    s[key + "_mean"] = statistics.mean(vals) if vals else float("nan")
                    s[key + "_sd"] = statistics.stdev(vals) if len(vals) > 1 else 0.0
                summary_rows.append(s)
    summary_path = ref / "healthy_two_pathologies_drug_summary.csv"; write_csv(summary_path, summary_rows)

    paired_summary = []
    for pathology_state in ["input_hyperexcitable", "nmda_excitotoxic"]:
        for drug, grid in GRIDS_NM.items():
            for concentration in grid:
                grp = [r for r in metrics_rows if r["pathology_state"] == pathology_state and r["drug"] == drug and r["concentration_nM"] == concentration]
                srow = {"pathology_state": pathology_state, "drug": drug, "concentration_nM": concentration, "n": len(grp)}
                for key in ["healthy_baseline_hz", "healthy_drug_hz", "pathological_baseline_hz", "pathological_drug_hz", "pathological_activity_normalization", "pathological_suppression_fraction", "healthy_suppression_fraction", "therapeutic_selectivity_index_exploratory"]:
                    vals = [float(r[key]) for r in grp if not math.isnan(float(r[key]))]
                    srow[key + "_mean"] = statistics.mean(vals) if vals else float("nan")
                    srow[key + "_sd"] = statistics.stdev(vals) if len(vals) > 1 else 0.0
                paired_summary.append(srow)
    paired_summary_path = ref / "healthy_two_pathologies_paired_summary.csv"; write_csv(paired_summary_path, paired_summary)

    print(f"wrote {rep_path.name}: {len(rows)} rows")
    print(f"wrote {metric_path.name}: {len(metrics_rows)} rows")
    print(f"wrote {summary_path.name}: {len(summary_rows)} rows")
    print(f"wrote {paired_summary_path.name}: {len(paired_summary)} rows")

    print("\nBest normalization (maximum mean PAN):")
    for pathology_state in ["input_hyperexcitable", "nmda_excitotoxic"]:
        print(pathology_state)
        for drug in GRIDS_NM:
            cand = [r for r in paired_summary if r["pathology_state"] == pathology_state and r["drug"] == drug and r["concentration_nM"] > 0]
            best = max(cand, key=lambda r: r["pathological_activity_normalization_mean"])
            print(" ", drug, best["concentration_nM"], "nM", "PAN=", round(best["pathological_activity_normalization_mean"],3), "healthy=", round(best["healthy_drug_hz_mean"],3), "Hz", "path=", round(best["pathological_drug_hz_mean"],3), "Hz")


if __name__ == "__main__":
    main()
