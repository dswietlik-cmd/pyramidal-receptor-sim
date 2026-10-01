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
