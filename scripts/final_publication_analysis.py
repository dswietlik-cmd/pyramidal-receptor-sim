#!/usr/bin/env python3
"""Final publication analysis for healthy vs two matched pathological states.

Uses the n=20 paired dataset produced by compare_two_pathologies_drugs_fast.py.
Outputs state-specific 4PL fits with paired-seed bootstrap, publication tables,
and a concise Results draft. Publication figures are generated separately by
make_professional_figures.py.
"""
from __future__ import annotations

import csv
import math
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "results" / "reference"
OUT = ROOT / "results" / "publication_final"
OUT.mkdir(parents=True, exist_ok=True)

REPS = REF / "healthy_two_pathologies_drug_replicates.csv"
PAIRED = REF / "healthy_two_pathologies_paired_metrics.csv"
BOOTSTRAP_N = 500
BOOTSTRAP_SEED = 20261001
STATES = ["healthy", "input_hyperexcitable", "nmda_excitotoxic"]
PATHOLOGIES = ["input_hyperexcitable", "nmda_excitotoxic"]
DRUGS = ["perampanel", "memantine", "diazepam"]
LABELS = {
    "healthy": "Healthy",
    "input_hyperexcitable": "Input-driven",
    "nmda_excitotoxic": "NMDA-driven",
    "perampanel": "Perampanel",
    "memantine": "Memantine",
    "diazepam": "Diazepam",
}


def four_pl(x, bottom, top, ec50, hill):
    x = np.asarray(x, dtype=float)
    return bottom + (top - bottom) * (x ** hill) / (ec50 ** hill + x ** hill)


def fit_curve(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    nz = x[x > 0]
    p0 = [0.0, min(1.0, float(np.nanmax(y))), float(np.median(nz)), 2.0]
    pars, _ = curve_fit(
        four_pl, x, y, p0=p0,
        bounds=([0, .05, .001, .1], [.25, 1.10, 10000, 25]),
        maxfev=100000,
    )
    pred = four_pl(x, *pars)
    ss_res = float(np.sum((y - pred) ** 2)); ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return pars, r2


def q(a, p):
    arr = np.asarray([v for v in a if np.isfinite(v)], float)
    return float(np.percentile(arr, p)) if len(arr) else float("nan")


def target_conc(pars, target):
    bottom, top, ec50, hill = pars
    if not (bottom < target < top): return float("nan")
    return float(ec50 * ((target-bottom)/(top-target)) ** (1.0/hill))


def write_csv(path, rows):
    if not rows: return
    with open(path, "w", newline="", encoding="utf-8") as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)


def main():
    df = pd.read_csv(REPS)
    pm = pd.read_csv(PAIRED)

    # Baseline table (one row per state x seed is duplicated by drug/concentration in source;
    # derive from concentration=0 and one drug only).
    base = df[(df.drug=="perampanel") & (df.concentration_nM==0)].copy()
    baseline_rows=[]
    for state in STATES:
        g=base[base.state==state]
        baseline_rows.append({
            "state": state,
            "label": LABELS[state],
            "n": len(g),
            "firing_hz_mean": g.firing_hz.mean(),
            "firing_hz_sd": g.firing_hz.std(ddof=1),
            "mean_isi_ms_mean": g.mean_isi_ms.mean(),
            "nmda_open_fraction_mean": g.nmda_open_fraction.mean(),
            "mean_plasticity_state_mean": g.mean_plasticity_state.mean(),
        })
    write_csv(OUT/"Table1_baseline_states.csv", baseline_rows)

    # State-specific concentration response fits. Use within-state suppression and paired-seed bootstrap.
    rng=np.random.default_rng(BOOTSTRAP_SEED)
    fit_rows=[]; target_rows=[]
    curves={}
    for drug in DRUGS:
        for state in STATES:
            g=df[(df.drug==drug)&(df.state==state)].copy()
            reps=sorted(g.replicate.unique())
            concs=np.array(sorted(g.concentration_nM.unique()), float)
            lookup={(int(r.replicate), float(r.concentration_nM)): float(r.within_state_suppression_fraction) for _,r in g.iterrows()}
            y=np.array([np.mean([lookup[(int(rep),c)] for rep in reps]) for c in concs])
            try:
                pars,r2=fit_curve(concs,y)
            except Exception:
                pars=np.array([np.nan]*4);r2=np.nan
            boot=[]; btargets=defaultdict(list)
            for _ in range(BOOTSTRAP_N):
                sampled=rng.choice(reps,size=len(reps),replace=True)
                yb=np.array([np.mean([lookup[(int(rep),c)] for rep in sampled]) for c in concs])
                try: bp,_=fit_curve(concs,yb)
                except Exception: continue
                boot.append(bp)
                for t in (0.25,0.50,0.75): btargets[t].append(target_conc(bp,t))
            barr=np.asarray(boot,float) if boot else np.empty((0,4))
            max_observed=float(np.nanmax(y))
            ec50_in_range = pars[2] if max_observed >= 0.50 else np.nan
            fit_rows.append({
                "drug":drug,"state":state,"n_replicates":len(reps),"bootstrap_n":len(boot),
                "max_observed_suppression":max_observed,
                "bottom":pars[0],"Emax_4PL":pars[1],"functional_EC50_nM":ec50_in_range,
                "extrapolated_4PL_EC50_nM":pars[2],"Hill":pars[3],"R2":r2,
                "EC50_bootstrap_median_nM":float(np.median(barr[:,2])) if (len(barr) and max_observed >= 0.50) else np.nan,
                "EC50_CI2.5_nM":q(barr[:,2],2.5) if (len(barr) and max_observed >= 0.50) else np.nan,
                "EC50_CI97.5_nM":q(barr[:,2],97.5) if (len(barr) and max_observed >= 0.50) else np.nan,
                "EC50_estimable_within_tested_range": bool(max_observed >= 0.50),
            })
            for t in (0.25,0.50,0.75):
                vals=btargets[t] if max_observed >= t else []
                target_rows.append({"drug":drug,"state":state,"suppression_target":t,
                                    "max_observed_suppression":max_observed,
                                    "concentration_median_nM":float(np.nanmedian(vals)) if vals else np.nan,
                                    "CI2.5_nM":q(vals,2.5),"CI97.5_nM":q(vals,97.5),"n_boot":sum(np.isfinite(vals)),
                                    "estimable_within_tested_range": bool(max_observed >= t)})
            curves[(drug,state)] = (concs, y, pars)
    write_csv(OUT/"Table2_state_specific_4PL_fits.csv", fit_rows)
    write_csv(OUT/"Table3_functional_target_concentrations.csv", target_rows)

    # Healthy-preserving therapeutic window on tested grids.
    tw=[]
    for pathology in PATHOLOGIES:
        for drug in DRUGS:
            g=pm[(pm.pathology_state==pathology)&(pm.drug==drug)].copy()
            for conc in sorted(g.concentration_nM.unique()):
                gg=g[g.concentration_nM==conc]
                hpres=(gg.healthy_drug_hz/gg.healthy_baseline_hz).mean()
                psupp=gg.pathological_suppression_fraction.mean()
                pan=gg.pathological_activity_normalization.mean()
                diff=(gg.pathological_suppression_fraction-gg.healthy_suppression_fraction).mean()
                tw.append({"pathology_state":pathology,"drug":drug,"concentration_nM":conc,
                           "healthy_preservation_fraction":hpres,"pathological_suppression_fraction":psupp,
                           "PAN_mean":pan,"selectivity_difference":diff,
                           "eligible_healthy_ge_0.80":bool(hpres>=0.80)})
    write_csv(OUT/"Table4_therapeutic_window_grid.csv", tw)
    best=[]
    for pathology in PATHOLOGIES:
        for drug in DRUGS:
            cand=[r for r in tw if r["pathology_state"]==pathology and r["drug"]==drug and r["eligible_healthy_ge_0.80"] and r["concentration_nM"]>0]
            if cand:
                b=max(cand,key=lambda r:r["pathological_suppression_fraction"])
                best.append(b)
    write_csv(OUT/"Table5_best_effect_with_healthy_preservation_ge80.csv", best)

    # Publication figures are generated separately by scripts/make_professional_figures.py
    # (English) and scripts/make_professional_figures_pl.py (Polish). Keeping figure
    # generation out of this analysis script prevents obsolete flat figure files from
    # being recreated in results/publication_final/.

    # Results draft
    bmap={r['state']:r for r in baseline_rows}
    lines=[]
    lines.append('# Results draft — final n=20 experiment\n')
    lines.append('## Baseline matching')
    lines.append(f"The healthy state fired at {bmap['healthy']['firing_hz_mean']:.2f} ± {bmap['healthy']['firing_hz_sd']:.2f} Hz. The input-driven and NMDA-driven pathological states were closely matched for baseline hyperexcitability at {bmap['input_hyperexcitable']['firing_hz_mean']:.2f} ± {bmap['input_hyperexcitable']['firing_hz_sd']:.2f} Hz and {bmap['nmda_excitotoxic']['firing_hz_mean']:.2f} ± {bmap['nmda_excitotoxic']['firing_hz_sd']:.2f} Hz, respectively (n=20 paired seeds per state).\n")
    lines.append('## State-dependent pharmacodynamic response')
    for drug in DRUGS:
        fr=[r for r in fit_rows if r['drug']==drug]
        lines.append(f"**{LABELS[drug]}.** " + '; '.join([f"{LABELS[r['state']]} functional EC50 {r['functional_EC50_nM']:.2f} nM (bootstrap 95% CI {r['EC50_CI2.5_nM']:.2f}–{r['EC50_CI97.5_nM']:.2f})" if np.isfinite(r['functional_EC50_nM']) else f"{LABELS[r['state']]} EC50 not estimable within the tested concentration range (maximum observed suppression {100*r['max_observed_suppression']:.1f}%)" for r in fr]) + '.\n')
    lines.append('## Healthy-preserving window')
    for pathology in PATHOLOGIES:
        lines.append(f"For {LABELS[pathology].lower()} pathology, the largest mean pathological suppression observed on the tested grid while preserving at least 80% of mean healthy firing was:")
        for drug in DRUGS:
            rows=[r for r in best if r['pathology_state']==pathology and r['drug']==drug]
            if rows:
                r=rows[0]; lines.append(f"- {LABELS[drug]}: {100*r['pathological_suppression_fraction']:.1f}% suppression at {r['concentration_nM']:g} nM, with {100*r['healthy_preservation_fraction']:.1f}% healthy firing preserved.")
            else: lines.append(f"- {LABELS[drug]}: no non-zero tested concentration met the >=80% healthy-preservation criterion.")
        lines.append('')
    lines.append('These functional EC50 values are properties of the calibrated drug–receptor–neuron system and should not be interpreted as molecular binding constants or clinical exposure targets. The healthy-preserving analysis is restricted to the tested concentration grids and is therefore descriptive rather than a clinical therapeutic-index estimate.')
    (OUT/'Results_final_n20_two_pathologies.md').write_text('\n'.join(lines),encoding='utf-8')

    # Methods note
    methods=f'''# Final experiment specification\n\n- States: healthy; input-driven moderate hyperexcitability (1.65x excitatory drive); NMDA-driven moderate excitotoxicity-like state (4.10x NMDA pathology multiplier).\n- Replicates: n=20 paired seeds per state and concentration.\n- Base seed: 20261200.\n- Input jitter: 0.5 ms.\n- Drugs: perampanel, memantine, diazepam.\n- Concentration grids are defined in scripts/compare_two_pathologies_drugs_fast.py.\n- Primary functional outcome: firing rate.\n- Secondary model outcomes: ISI, AMPA/NMDA/GABA PSP components, NMDA gate-open fraction, and mean plasticity state.\n- Pathological Activity Normalization (PAN) is computed against same-seed healthy and pathological baselines.\n- A descriptive healthy-preserving window is defined on the tested concentration grid as concentrations preserving >=80% of mean healthy baseline firing.\n- State-specific functional concentration-response curves are summarized with a 4-parameter logistic model and 500 paired-seed bootstrap resamples.\n'''
    (OUT/'METHODS_FINAL_EXPERIMENT.md').write_text(methods,encoding='utf-8')

    print('wrote final publication outputs to',OUT)

if __name__=='__main__':
    main()
