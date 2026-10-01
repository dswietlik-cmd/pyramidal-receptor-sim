#!/usr/bin/env python3
from __future__ import annotations
import csv, math, random, statistics, sys
from pathlib import Path
import numpy as np
from numba import njit

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from pyramidal_receptor_sim import Neuron, configure_neuron, load_config

CFG=load_config(ROOT/'configs/ca1_reference.txt')
NE=int(CFG['N_EXCITATORY_INPUTS']); NI=int(CFG['N_INHIBITORY_INPUTS']); N=NE+NI
DT=float(CFG['DT_MS']); STEPS=int(CFG['STEPS'])
BASE_FREQS=np.array([float(CFG[f'INPUT_{i}']) for i in range(1,N+1)],dtype=np.float64)
JITTER_MS=0.5; NREP=20; BASE_SEED=20261200
DRIVE_GRID=np.round(np.arange(1.0,2.501,0.05),2)
TARGETS=(18.0,24.0,30.0)

@njit(cache=True)
def influence(source,dest,N,LSW):
    if dest>=source: return 1.0-(1.0-LSW)/(N-1.0)*(dest-source)
    if (source-dest)<=N/4.0: return 1.0-4.0/N*(source-dest)
    return 0.0

@njit(cache=True)
def simulate(ex,inh, EPSPd,IPSPd,CaMT,NMDA_scale,ReP,RePB,Threshold,KEP,CRefrak,FQ,Clog,SleepC,Rsp,LSW):
    ExTab=np.zeros(33); MexTab=np.zeros(33); InhTab=np.zeros(33)
    ExTab[0]=ExTab[1]=ExTab[2]=ExTab[32]=0.0
    for i in range(1,5): ExTab[2+i]=EPSPd/4.0*i
    for i in range(1,26): ExTab[6+i]=EPSPd/26.0*(26-i)
    MexTab[0]=MexTab[1]=MexTab[2]=MexTab[32]=0.0
    for i in range(1,5): MexTab[2+i]=EPSPd/20.0*i
    for i in range(1,26): MexTab[6+i]=EPSPd/99.0*(26-i)
    InhTab[0]=InhTab[1]=InhTab[2]=0.0
    for i in range(1,5): InhTab[2+i]=IPSPd/4.0*i
    for i in range(1,20): InhTab[6+i]=IPSPd/20.0*(20-i)
    STEx=np.full((NE,33),ReP); STMex=np.full((NE,33),ReP); STInh=np.full((NI,33),ReP)
    MemVol=np.zeros(NE); SpinePot=np.full(N,ReP); CaPath=np.zeros(NE,np.int64)
    InhDel=np.zeros(NI,np.int64); SleepCout=np.zeros(NI,np.int64)
    Refrak=0; PSP=ReP; SleepState=0; spikes=0; last_sp=-1.0; isi_sum=0.0; isi_n=0
    ampa_sum=0.0; nmda_sum=0.0; gaba_sum=0.0; gate_sum=0.0; plast_sum=0.0
    for step in range(ex.shape[0]):
        for ii in range(NE):
            if ex[step,ii]!=0:
                adaptation=((ReP-Rsp)-(ReP-SpinePot[ii]))/(ReP-Rsp)
                memory=1.0+(1.0/6.0)*math.log(MemVol[ii]+1.0)/Clog
                factor=adaptation*memory
                for k in range(33): STEx[ii,k]+=factor*ExTab[k]
        for ii in range(NI):
            if inh[step,ii]!=0:
                InhDel[ii]=6
                for k in range(33): STInh[ii,k]+=InhTab[k]
        for ii in range(NE):
            if ex[step,ii]!=0 and CaPath[ii]==1:
                for k in range(33): STMex[ii,k]+=MexTab[k]
            power=(STMex[ii,0]-ReP)*6.0
            MemVol[ii]+=math.exp(power)-1.0
            if MemVol[ii]>FQ: MemVol[ii]-=FQ
            else: MemVol[ii]=0.0
            number=ii+1; pot=ReP
            for m in range(NE): pot+=(STEx[m,0]-ReP)*influence(m+1,number,N,LSW)
            for m in range(NI): pot+=(STInh[m,0]-ReP)*influence(m+1+NE,number,N,LSW)
            SpinePot[ii]=pot; CaPath[ii]=1 if pot>=CaMT else 0
        for jj in range(NI):
            number=NE+jj+1; pot=ReP
            for m in range(NE): pot+=(STEx[m,0]-ReP)*influence(m+1,number,N,LSW)
            for m in range(NI): pot+=(STInh[m,0]-ReP)*influence(m+1+NE,number,N,LSW)
            SpinePot[NE+jj]=pot
        SOut=0; ampa_comp=0.0; nmda_comp=0.0; gaba_comp=0.0
        if Refrak==0:
            PSP=ReP
            for ii in range(NE):
                number=ii+1; weight=(1.0-LSW)/(NE-1.0)*(number-NE)+1.0
                a=weight*(STEx[ii,0]-ReP); n=NMDA_scale*weight*(STMex[ii,0]-ReP)
                ampa_comp += a; nmda_comp += n; PSP += a+n
            for ii in range(NI):
                g=STInh[ii,0]-ReP; gaba_comp += g; PSP += g
            if PSP<KEP: PSP=KEP
            if PSP>=Threshold:
                SOut=1; Refrak=CRefrak
                for ii in range(NE):
                    number=ii+1; adrv=ReP-((ReP-RePB)/NE)*(NE-number)
                    for k in range(33): STEx[ii,k]=adrv
                for ii in range(NI):
                    for k in range(33): STInh[ii,k]=ReP
            else:
                for ii in range(NE):
                    for k in range(32): STEx[ii,k]=STEx[ii,k+1]
                    for k in range(32): STMex[ii,k]=STMex[ii,k+1]
                if SleepState==0:
                    for ii in range(NI):
                        for k in range(32): STInh[ii,k]=STInh[ii,k+1]
                else:
                    for ii in range(NI):
                        if InhDel[ii]>0:
                            InhDel[ii]-=1
                            for k in range(32): STInh[ii,k]=STInh[ii,k+1]
                        else:
                            SleepCout[ii]+=1
                            if SleepCout[ii]==SleepC:
                                SleepCout[ii]=0
                                for k in range(32): STInh[ii,k]=STInh[ii,k+1]
        else:
            for ii in range(NE):
                for k in range(32): STEx[ii,k]=STEx[ii,k+1]
                for k in range(32): STMex[ii,k]=STMex[ii,k+1]
            if SleepState==0:
                for ii in range(NI):
                    for k in range(32): STInh[ii,k]=STInh[ii,k+1]
            else:
                for ii in range(NI):
                    if InhDel[ii]>0:
                        InhDel[ii]-=1
                        for k in range(32): STInh[ii,k]=STInh[ii,k+1]
                    else:
                        SleepCout[ii]+=1
                        if SleepCout[ii]==SleepC:
                            SleepCout[ii]=0
                            for k in range(32): STInh[ii,k]=STInh[ii,k+1]
            if Refrak==(CRefrak-1): PSP=ReP
            Refrak-=1
        SleepState=1 if PSP<ReP else 0
        if SOut==1:
            spikes+=1; tm=step*DT
            if last_sp>=0: isi_sum += tm-last_sp; isi_n += 1
            last_sp=tm
        ampa_sum += ampa_comp; nmda_sum += nmda_comp; gaba_sum += gaba_comp
        gate_sum += np.sum(CaPath)/NE; plast_sum += np.mean(MemVol)
    dur=ex.shape[0]*DT/1000.0
    return spikes/dur, (isi_sum/isi_n if isi_n>0 else np.nan), ampa_sum/ex.shape[0], nmda_sum/ex.shape[0], gaba_sum/ex.shape[0], gate_sum/ex.shape[0], plast_sum/ex.shape[0]

def make_stream(seed,drive):
    freqs=BASE_FREQS.copy(); freqs[:NE]*=drive
    r=random.Random(seed); js=int(round(JITTER_MS/DT)); ex=np.zeros((STEPS,NE),np.uint8); inh=np.zeros((STEPS,NI),np.uint8)
    for idx,f in enumerate(freqs):
        if f<=0: continue
        per=max(1,int(round(1000.0/(f*DT))))
        for base in range(0,STEPS,per):
            t=max(0,min(STEPS-1,base+(r.randint(-js,js) if js else 0)))
            if idx<NE: ex[t,idx]=1
            else: inh[t,idx-NE]=1
    return ex,inh

def pars():
    return (float(CFG['AMPA_EPSP_AMPLITUDE_MV']),float(CFG['GABAA_IPSP_AMPLITUDE_MV']),float(CFG['NMDA_ACTIVATION_THRESHOLD_MV']),float(CFG['NMDA_PSP_SCALE']),float(CFG['RESTING_POTENTIAL_MV']),float(CFG['RESET_POTENTIAL_UPPER_MV']),float(CFG['SPIKE_THRESHOLD_MV']),float(CFG['MINIMUM_MEMBRANE_POTENTIAL_MV']),int(CFG['REFRACTORY_STEPS']),float(CFG['PLASTICITY_DECAY_STEP']),float(CFG['PLASTICITY_LOG_SCALE']),int(CFG['INHIBITORY_DECAY_INTERVAL']),float(CFG['ADAPTATION_REFERENCE_POTENTIAL_MV']),float(CFG['MINIMUM_SYNAPTIC_WEIGHT']))

def py_hz(seed,drive):
    from pyramidal_receptor_sim.inputs import jittered_stream
    freqs=BASE_FREQS.copy(); freqs[:NE]*=drive
    ex,inh=jittered_stream(freqs,STEPS,DT,NE,NI,seed,JITTER_MS)
    n=configure_neuron(Neuron(),CFG); s=0
    for t in range(STEPS): s+=n.step(ex[t].tolist(),inh[t].tolist())
    return s/(STEPS*DT/1000.0)

def main():
    p=pars(); ex,inh=make_stream(BASE_SEED,1.0); _=simulate(ex,inh,*p)
    fast=simulate(ex,inh,*p)[0]; py=py_hz(BASE_SEED,1.0)
    print('equivalence baseline:',fast,py)
    if abs(fast-py)>1e-9: raise SystemExit('FAST/PYTHON MISMATCH')
    rows=[]
    for drive in DRIVE_GRID:
        vals=[]
        for rep in range(NREP):
            seed=BASE_SEED+rep; ex,inh=make_stream(seed,float(drive)); out=simulate(ex,inh,*p)
            row={'drive_multiplier':float(drive),'replicate':rep+1,'seed':seed,'firing_hz':out[0],'mean_isi_ms':out[1],'mean_ampa_psp_mv':out[2],'mean_nmda_psp_mv':out[3],'mean_gaba_psp_mv':out[4],'nmda_open_fraction':out[5],'mean_plasticity_state':out[6]}
            rows.append(row); vals.append(out[0])
        print(f'{drive:.2f}: {statistics.mean(vals):.3f} +/- {statistics.stdev(vals):.3f} Hz')
    ref=ROOT/'results/reference'; ref.mkdir(parents=True,exist_ok=True)
    with (ref/'hyperexcitability_calibration_replicates.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
    summary=[]
    metrics=['firing_hz','mean_isi_ms','mean_ampa_psp_mv','mean_nmda_psp_mv','mean_gaba_psp_mv','nmda_open_fraction','mean_plasticity_state']
    for drive in DRIVE_GRID:
        grp=[r for r in rows if r['drive_multiplier']==float(drive)]; s={'drive_multiplier':float(drive),'n':len(grp)}
        for k in metrics:
            vals=[r[k] for r in grp if not (isinstance(r[k],float) and math.isnan(r[k]))]
            s[k+'_mean']=statistics.mean(vals); s[k+'_sd']=statistics.stdev(vals) if len(vals)>1 else 0.0
        summary.append(s)
    with (ref/'hyperexcitability_calibration_summary.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=summary[0].keys()); w.writeheader(); w.writerows(summary)
    selected=[]
    for target in TARGETS:
        best=min(summary,key=lambda r:abs(r['firing_hz_mean']-target))
        selected.append({'state':{18.0:'mild',24.0:'moderate',30.0:'severe'}[target],'target_firing_hz':target,**best})
    with (ref/'hyperexcitability_selected_states.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=selected[0].keys()); w.writeheader(); w.writerows(selected)
    print('SELECTED')
    for r in selected: print(r['state'],r['drive_multiplier'],r['firing_hz_mean'],r['firing_hz_sd'])

if __name__=='__main__': main()
