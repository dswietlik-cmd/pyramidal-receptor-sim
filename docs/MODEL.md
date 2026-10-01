# Model equations and implementation

## 1. Biological scope

The model represents a CA1-like pyramidal neuron receiving 13 excitatory and 3 inhibitory inputs.
Excitatory transmission is decomposed into AMPA- and NMDA-mediated components. Inhibitory
transmission is represented by a GABA-A-mediated component. The reference simulation step is
0.5 ms.

The model is intentionally receptor-centered: pharmacological interventions first modify receptor
function, and neuronal firing is treated as a downstream functional readout.

## 2. Synaptic response kernels

For each excitatory input `i`, the model stores separate time-dependent response registers
`AMPA_i(t)` and `NMDA_i(t)`. Each inhibitory input has a `GABAA_i(t)` register. Incoming spikes
insert a predefined rise-and-decay kernel into the corresponding register.

Receptor gains are dimensionless:

```text
g_AMPA = 1.0    drug-free AMPA reference
g_NMDA = 1.0    drug-free NMDA reference
g_GABAA = 1.0   drug-free GABA-A reference
```

Antagonism reduces the relevant excitatory gain below 1.0. Positive allosteric modulation of GABA-A
can increase `g_GABAA` above 1.0.

## 3. Somatic postsynaptic potential

Outside the refractory period, the model computes

```text
V_post(t) = V_rest + PSP_AMPA(t) + PSP_NMDA(t) + PSP_GABAA(t)
```

The AMPA term is

```text
PSP_AMPA(t) = sum_i w_i [AMPA_i(t) - V_rest]
```

where the AMPA gain is applied when the AMPA response kernel is inserted.

The explicit NMDA term is

```text
PSP_NMDA(t) = s_NMDA * sum_i w_i [NMDA_i(t) - V_rest]
```

where `s_NMDA = nmda_psp_scale` is the externally calibrated coupling between the NMDA register and
the somatic potential. The NMDA gain is applied when the NMDA kernel is inserted.

The inhibitory term is

```text
PSP_GABAA(t) = sum_j [GABAA_j(t) - V_rest]
```

with the GABA-A gain applied during inhibitory-kernel insertion.

## 4. Voltage-dependent NMDA gating

The NMDA component is not inserted solely because an excitatory presynaptic spike occurs. The model
first evaluates a local spine potential. NMDA activation is permitted when

```text
V_spine >= V_NMDA_threshold
```

where `V_NMDA_threshold` is configured by `NMDA_ACTIVATION_THRESHOLD_MV` and controls removal of the modeled magnesium block. This
provides the voltage dependence needed for NMDA activation.

The previous-step gating state determines whether an incoming excitatory event inserts an NMDA
response. The gating state is then recalculated for the next step.

## 5. NMDA-coupled plasticity state

NMDA activity also updates an internal memory/plasticity variable `M_i`. Its current implementation
uses an exponential accumulation followed by a forgetting term. The synaptic gain used for later
excitatory input is

```text
plasticity_gain_i = 1 + ln(M_i + 1) / (6 * plasticity_log_scale)
```

Thus NMDA activity has two distinct modeled consequences:

1. an acute explicit contribution to the postsynaptic potential;
2. a slower contribution to synaptic plasticity through the memory state.

These mechanisms should be reported separately in analyses.

## 6. Spike generation and refractory state

If the somatic potential reaches the firing threshold,

```text
V_post >= V_spike_threshold
```

a spike is emitted and the model enters a refractory period configured by `REFRACTORY_STEPS`.
Synaptic registers continue to evolve during refractoriness so that synaptic time courses are not
artificially frozen.

## 7. Drug-receptor mapping

Drug concentration-response parameters are defined in `drugs.py`, separate from the neuronal core.
For an antagonist,

```text
B(C) = C^n / (IC50^n + C^n)
g_receptor(C) = 1 - B(C)
```

For the GABA-A positive allosteric modulator,

```text
H(C) = C^n / (EC50^n + C^n)
g_GABAA(C) = 1 + Emax_increase * H(C)
```

The concentration and IC50/EC50 are converted to the same units before the Hill relationship is
evaluated.

## 8. Functional readouts

The reference workflows record:

- spike count and firing rate;
- interspike interval;
- somatic postsynaptic potential;
- AMPA, NMDA and GABA-A PSP components;
- receptor gains;
- NMDA voltage-gating state and plasticity variables when requested.

A functional EC50 estimated from firing suppression is an emergent systems-level quantity and should
not be interpreted as a receptor-binding EC50 or IC50.
