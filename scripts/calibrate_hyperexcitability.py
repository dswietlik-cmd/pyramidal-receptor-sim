#!/usr/bin/env python3
from __future__ import annotations
import csv, math, statistics
from pathlib import Path
from pyramidal_receptor_sim import Neuron, configure_neuron, load_config
from pyramidal_receptor_sim.inputs import jittered_stream

N_REPLICATES=20
BASE_SEED=20261200
JITTER_MS=0.5
DRIVE_GRID=[round(1.0+i*0.05,2) for i in range(0,31)]  # 1.00..2.50
TARGETS=[18.0,24.0,30.0]

def simulate(neuron, ex, inh, dt_ms):
    spikes=0; spike_times=[]; ampa=[]; nmda=[]; gaba=[]; gate_open=0; gate_total=0; plast=[]
    for t in range(ex.shape[0]):
        if neuron.step(ex[t].tolist(), inh[t].tolist()):
            spikes += 1; spike_times.append(t*dt_ms)
        ampa.append(neuron.last_ampa_psp); nmda.append(neuron.last_nmda_psp); gaba.append(neuron.last_gaba_psp)
        gate_open += sum(neuron.nmda_gate_open[1:neuron.n_excitatory_inputs+1]); gate_total += neuron.n_excitatory_inputs
        plast.append(statistics.mean(neuron.plasticity_state[1:neuron.n_excitatory_inputs+1]))
    duration_s=ex.shape[0]*dt_ms/1000.0
    isi=[b-a for a,b in zip(spike_times, spike_times[1:])]
    return {
      'firing_hz': spikes/duration_s,
      'mean_isi_ms': statistics.mean(isi) if isi else float('nan'),
      'mean_ampa_psp_mv': statistics.mean(ampa),
      'mean_nmda_psp_mv': statistics.mean(nmda),
      'mean_gaba_psp_mv': statistics.mean(gaba),
      'nmda_open_fraction': gate_open/gate_total if gate_total else float('nan'),
      'mean_plasticity_state': statistics.mean(plast),
    }

def main():
    root=Path(__file__).resolve().parents[1]
    cfg=load_config(root/'configs/ca1_reference.txt')
    dt=float(cfg['DT_MS']); steps=int(cfg['STEPS']); ne=int(cfg['N_EXCITATORY_INPUTS']); ni=int(cfg['N_INHIBITORY_INPUTS'])
    base=[float(cfg[f'INPUT_{i}']) for i in range(1,ne+ni+1)]
    rows=[]
    for drive in DRIVE_GRID:
      freqs=[f*drive if idx<ne else f for idx,f in enumerate(base)]
      for rep in range(N_REPLICATES):
        seed=BASE_SEED+rep
        ex,inh=jittered_stream(freqs,steps,dt,ne,ni,seed,JITTER_MS)
        n=configure_neuron(Neuron(),cfg)
        out=simulate(n,ex,inh,dt)
        rows.append({'drive_multiplier':drive,'replicate':rep+1,'seed':seed,**out})
      print('drive',drive,'done')
    rep_path=root/'results/reference/hyperexcitability_calibration_replicates.csv'
    rep_path.parent.mkdir(parents=True,exist_ok=True)
    with rep_path.open('w',newline='',encoding='utf-8') as f:
      w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
    summary=[]
    for drive in DRIVE_GRID:
      grp=[r for r in rows if r['drive_multiplier']==drive]
      item={'drive_multiplier':drive}
      for k in ['firing_hz','mean_isi_ms','mean_ampa_psp_mv','mean_nmda_psp_mv','mean_gaba_psp_mv','nmda_open_fraction','mean_plasticity_state']:
        vals=[r[k] for r in grp if not (isinstance(r[k],float) and math.isnan(r[k]))]
        item[k+'_mean']=statistics.mean(vals) if vals else float('nan')
        item[k+'_sd']=statistics.stdev(vals) if len(vals)>1 else 0.0
      summary.append(item)
    sum_path=root/'results/reference/hyperexcitability_calibration_summary.csv'
    with sum_path.open('w',newline='',encoding='utf-8') as f:
      w=csv.DictWriter(f,fieldnames=summary[0].keys()); w.writeheader(); w.writerows(summary)
    selected=[]
    for target in TARGETS:
      best=min(summary,key=lambda r:abs(r['firing_hz_mean']-target))
      selected.append({'target_firing_hz':target,**best})
    sel_path=root/'results/reference/hyperexcitability_selected_states.csv'
    with sel_path.open('w',newline='',encoding='utf-8') as f:
      w=csv.DictWriter(f,fieldnames=selected[0].keys()); w.writeheader(); w.writerows(selected)
    print('selected:')
    for r in selected: print(r['target_firing_hz'],r['drive_multiplier'],r['firing_hz_mean'],r['firing_hz_sd'])

if __name__=='__main__': main()
