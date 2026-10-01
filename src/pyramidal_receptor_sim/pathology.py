"""Phenotype-level pathological input-drive presets.

The presets change only excitatory input frequency. Receptor gains, inhibitory
input frequencies, membrane parameters, and synaptic amplitudes remain at their
calibrated reference values. This avoids defining pathology by directly altering
the pharmacological targets later tested in the study.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

@dataclass(frozen=True)
class HyperexcitabilityState:
    name: str
    excitatory_drive_multiplier: float
    target_firing_hz: float
    calibrated_mean_firing_hz: float
    calibrated_sd_firing_hz: float

HEALTHY = HyperexcitabilityState("healthy", 1.00, 12.0, 10.965, 0.906)
MILD = HyperexcitabilityState("mild", 1.30, 18.0, 18.055, 0.199)
MODERATE = HyperexcitabilityState("moderate", 1.65, 24.0, 24.605, 1.089)
SEVERE = HyperexcitabilityState("severe", 1.90, 30.0, 29.640, 0.289)

STATES = {s.name: s for s in (HEALTHY, MILD, MODERATE, SEVERE)}

def apply_excitatory_drive(frequencies: Iterable[float], n_excitatory_inputs: int, multiplier: float) -> list[float]:
    if multiplier <= 0:
        raise ValueError("excitatory drive multiplier must be > 0")
    freqs = list(map(float, frequencies))
    return [f * multiplier if i < n_excitatory_inputs else f for i, f in enumerate(freqs)]

def frequencies_for_state(frequencies: Iterable[float], n_excitatory_inputs: int, state: str | HyperexcitabilityState) -> list[float]:
    preset = STATES[state] if isinstance(state, str) else state
    return apply_excitatory_drive(frequencies, n_excitatory_inputs, preset.excitatory_drive_multiplier)


@dataclass(frozen=True)
class NMDAExcitotoxicityState:
    """Calibrated excitotoxicity-like phenotype driven by NMDA synaptic input."""
    name: str
    nmda_pathology_multiplier: float
    target_firing_hz: float
    calibrated_mean_firing_hz: float
    calibrated_sd_firing_hz: float

NMDA_HEALTHY = NMDAExcitotoxicityState("healthy", 1.00, 12.0, 11.015, 0.824)
NMDA_MILD = NMDAExcitotoxicityState("mild", 3.00, 18.0, 18.530, 1.931)
NMDA_MODERATE = NMDAExcitotoxicityState("moderate", 4.10, 24.0, 23.405, 5.556)
NMDA_SEVERE = NMDAExcitotoxicityState("severe", 4.70, 30.0, 30.375, 6.147)
NMDA_EXCITOTOXICITY_STATES = {s.name: s for s in (NMDA_HEALTHY, NMDA_MILD, NMDA_MODERATE, NMDA_SEVERE)}

def apply_nmda_excitotoxicity(neuron, state: str | NMDAExcitotoxicityState):
    preset = NMDA_EXCITOTOXICITY_STATES[state] if isinstance(state, str) else state
    neuron.set_nmda_pathology_multiplier(preset.nmda_pathology_multiplier)
    return neuron

@dataclass(frozen=True)
class GABAADisinhibitionState:
    """GABA-A hypofunction state defined by reduced postsynaptic inhibitory efficacy."""
    name: str
    gabaa_pathology_multiplier: float
    target_firing_hz: float | None
    calibrated_mean_firing_hz: float
    calibrated_sd_firing_hz: float

# Calibration in v0.9.5 showed that isolated loss of GABA-A efficacy cannot
# generate the pre-specified 18/24/30 Hz hyperexcitability phenotypes in the
# current reference model. Complete loss is retained as a boundary condition.
GABAA_HEALTHY = GABAADisinhibitionState("healthy", 1.00, 12.0, 10.965, 0.906)
GABAA_COMPLETE_LOSS = GABAADisinhibitionState("complete_loss", 0.00, None, 11.280, 0.738)
GABAA_DISINHIBITION_STATES = {s.name: s for s in (GABAA_HEALTHY, GABAA_COMPLETE_LOSS)}

def apply_gabaa_disinhibition(neuron, state: str | GABAADisinhibitionState):
    preset = GABAA_DISINHIBITION_STATES[state] if isinstance(state, str) else state
    neuron.set_gabaa_pathology_multiplier(preset.gabaa_pathology_multiplier)
    return neuron
