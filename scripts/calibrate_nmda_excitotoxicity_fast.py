#!/usr/bin/env python3
from __future__ import annotations
import csv, math, random, statistics, sys
from pathlib import Path
import numpy as np
from numba import njit
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/'src'))
from pyramidal_receptor_sim import load_config
CFG=load_config(ROOT/'configs/ca1_reference.txt')
NE=int(CFG['N_EXCITATORY_INPUTS']); NI=int(CFG['N_INHIBITORY_INPUTS']); N=NE+NI
DT=float(CFG['DT_MS']); STEPS=int(CFG['STEPS']); JITTER_MS=0.5; BASE_SEED=20261200
BASE_FREQS=np.array([float(CFG[f'INPUT_{i}']) for i in range(1,N+1)],dtype=np.float64)
@njit(cache=True)
def influence(source,dest,n_total,lsw):
    if dest>=source:return 1-(1-lsw)/(n_total-1)*(dest-source)
    if (source-dest)<=n_total/4:return 1-4/n_total*(source-dest)
    return 0.0
@njit(cache=True)
def sim(ex,inh,EPSPd,IPSPd,CaMT,NMDA_scale,ReP,RePB,Threshold,KEP,CRefrak,FQ,Clog,SleepC,Rsp,LSW,nmda_path):
    ExTab=np.zeros(33); MexTab=np.zeros(33); InhTab=np.zeros(33)
    for i in range(1,5): ExTab[2+i]=EPSPd/4*i; MexTab[2+i]=EPSPd/20*i; InhTab[2+i]=IPSPd/4*i
    for i in range(1,26): ExTab[6+i]=EPSPd/26*(26-i); MexTab[6+i]=EPSPd/99*(26-i)
    for i in range(1,20): InhTab[6+i]=IPSPd/20*(20-i)
    STEx=np.full((NE,33),ReP); STMex=np.full((NE,33),ReP); STInh=np.full((NI,33),ReP)
    MemVol=np.zeros(NE); SpinePot=np.full(N,ReP); CaPath=np.zeros(NE,np.int64); InhDel=np.zeros(NI,np.int64); SleepCout=np.zeros(NI,np.int64)
    Refrak=0; PSP=ReP; SleepState=0; spikes=0; last=-1.; isi_sum=0.; isi_n=0; ampa_sum=nmda_sum=gaba_sum=gate_sum=plast_sum=0.
    for step in range(ex.shape[0]):
        for ii in range(NE):
            if ex[step,ii]:
                adaptation=((ReP-Rsp)-(ReP-SpinePot[ii]))/(ReP-Rsp); memory=1+(1/6)*math.log(MemVol[ii]+1)/Clog; factor=adaptation*memory
                for k in range(33): STEx[ii,k]+=factor*ExTab[k]
        for ii in range(NI):
            if inh[step,ii]:
                InhDel[ii]=6
                for k in range(33):STInh[ii,k]+=InhTab[k]
        for ii in range(NE):
            if ex[step,ii] and CaPath[ii]==1:
                for k in range(33): STMex[ii,k]+=nmda_path*MexTab[k]
            power=(STMex[ii,0]-ReP)*6.; MemVol[ii]+=math.exp(power)-1
            if MemVol[ii]>FQ:MemVol[ii]-=FQ
            else:MemVol[ii]=0
            number=ii+1; pot=ReP
            for m in range(NE):pot+=(STEx[m,0]-ReP)*influence(m+1,number,N,LSW)
            for m in range(NI):pot+=(STInh[m,0]-ReP)*influence(m+1+NE,number,N,LSW)
            SpinePot[ii]=pot; CaPath[ii]=1 if pot>=CaMT else 0
        for jj in range(NI):
            number=NE+jj+1; pot=ReP
            for m in range(NE):pot+=(STEx[m,0]-ReP)*influence(m+1,number,N,LSW)
            for m in range(NI):pot+=(STInh[m,0]-ReP)*influence(m+1+NE,number,N,LSW)
            SpinePot[NE+jj]=pot
        SOut=0; ac=nc=gc=0.
        if Refrak==0:
            PSP=ReP
            for ii in range(NE):
                number=ii+1; w=(1-LSW)/(NE-1)*(number-NE)+1; a=w*(STEx[ii,0]-ReP); nm=NMDA_scale*w*(STMex[ii,0]-ReP); ac+=a; nc+=nm; PSP+=a+nm
            for ii in range(NI):g=STInh[ii,0]-ReP;gc+=g;PSP+=g
            if PSP<KEP:PSP=KEP
            if PSP>=Threshold:
                SOut=1;Refrak=CRefrak
                for ii in range(NE):
                    number=ii+1;adrv=ReP-((ReP-RePB)/NE)*(NE-number)
                    for k in range(33):STEx[ii,k]=adrv
                for ii in range(NI):
                    for k in range(33):STInh[ii,k]=ReP
            else:
                for ii in range(NE):
                    for k in range(32):STEx[ii,k]=STEx[ii,k+1];STMex[ii,k]=STMex[ii,k+1]
                for ii in range(NI):
                    for k in range(32):STInh[ii,k]=STInh[ii,k+1]
        else:
            for ii in range(NE):
                for k in range(32):STEx[ii,k]=STEx[ii,k+1];STMex[ii,k]=STMex[ii,k+1]
            for ii in range(NI):
                for k in range(32):STInh[ii,k]=STInh[ii,k+1]
            if Refrak==(CRefrak-1):PSP=ReP
            Refrak-=1
        SleepState=1 if PSP<ReP else 0
        if SOut:
            spikes+=1;tm=step*DT
            if last>=0:isi_sum+=tm-last;isi_n+=1
            last=tm
        ampa_sum+=ac;nmda_sum+=nc;gaba_sum+=gc;gate_sum+=np.sum(CaPath)/NE;plast_sum+=np.mean(MemVol)
    dur=ex.shape[0]*DT/1000
    return spikes/dur, isi_sum/isi_n if isi_n else np.nan,ampa_sum/ex.shape[0],nmda_sum/ex.shape[0],gaba_sum/ex.shape[0],gate_sum/ex.shape[0],plast_sum/ex.shape[0]
def pars():return (float(CFG['AMPA_EPSP_AMPLITUDE_MV']),float(CFG['GABAA_IPSP_AMPLITUDE_MV']),float(CFG['NMDA_ACTIVATION_THRESHOLD_MV']),float(CFG['NMDA_PSP_SCALE']),float(CFG['RESTING_POTENTIAL_MV']),float(CFG['RESET_POTENTIAL_UPPER_MV']),float(CFG['SPIKE_THRESHOLD_MV']),float(CFG['MINIMUM_MEMBRANE_POTENTIAL_MV']),int(CFG['REFRACTORY_STEPS']),float(CFG['PLASTICITY_DECAY_STEP']),float(CFG['PLASTICITY_LOG_SCALE']),int(CFG['INHIBITORY_DECAY_INTERVAL']),float(CFG['ADAPTATION_REFERENCE_POTENTIAL_MV']),float(CFG['MINIMUM_SYNAPTIC_WEIGHT']))
def stream(seed):
    rng=random.Random(seed);js=int(round(JITTER_MS/DT));ex=np.zeros((STEPS,NE),np.uint8);inh=np.zeros((STEPS,NI),np.uint8)
    for idx,freq in enumerate(BASE_FREQS):
        period=max(1,int(round(1000/(freq*DT))))
        for base in range(0,STEPS,period):
            t=max(0,min(STEPS-1,base+(rng.randint(-js,js) if js else 0)))
            if idx<NE:ex[t,idx]=1
            else:inh[t,idx-NE]=1
    return ex,inh
def write(path,rows):
    with path.open('w',newline='',encoding='utf-8') as f:w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
def main():
    factors=[1.0,3.0,4.1,4.7]
    p=pars(); rows=[]
    for rep in range(20):
        ex,inh=stream(BASE_SEED+rep)
        for fac in factors:
            o=sim(ex,inh,*p,fac);rows.append({'nmda_pathology_multiplier':fac,'replicate':rep+1,'seed':BASE_SEED+rep,'firing_hz':o[0],'mean_isi_ms':o[1],'mean_ampa_psp_mv':o[2],'mean_nmda_psp_mv':o[3],'mean_gaba_psp_mv':o[4],'nmda_open_fraction':o[5],'mean_plasticity_state':o[6]})
    ref=ROOT/'results'/'reference';ref.mkdir(parents=True,exist_ok=True);write(ref/'nmda_excitotoxicity_calibration_replicates.csv',rows)
    summary=[]
    for fac in factors:
        grp=[r for r in rows if r['nmda_pathology_multiplier']==fac]
        sr={'nmda_pathology_multiplier':fac,'n':len(grp)}
        for key in ['firing_hz','mean_isi_ms','mean_ampa_psp_mv','mean_nmda_psp_mv','mean_gaba_psp_mv','nmda_open_fraction','mean_plasticity_state']:
            vals=[float(r[key]) for r in grp if not math.isnan(float(r[key]))]
            sr[key+'_mean']=statistics.mean(vals) if vals else float('nan')
            sr[key+'_sd']=statistics.stdev(vals) if len(vals)>1 else 0.0
        summary.append(sr)
    write(ref/'nmda_excitotoxicity_calibration_summary.csv',summary)
    print('factor meanHz')
    for fac in factors:
        vals=[r['firing_hz'] for r in rows if r['nmda_pathology_multiplier']==fac];print(fac,statistics.mean(vals))
if __name__=='__main__':main()
