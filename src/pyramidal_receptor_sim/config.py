"""Reference-configuration loading and application."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping

from .model import Neuron


def parse_number(text: str):
    value = text.strip().replace(",", ".")
    if any(c in value.lower() for c in (".", "e")):
        return float(value)
    return int(value)


def load_config(path: str | Path) -> dict:
    config = {}
    with Path(path).open("r", encoding="utf-8-sig") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            key, value = line.split("=", 1)
            config[key.strip()] = parse_number(value.split("#", 1)[0].strip())
    return config


def configure_neuron(neuron: Neuron, cfg: Mapping[str, float | int]) -> Neuron:
    """Apply a reference configuration and reset all dynamic state."""
    neuron.plasticity_decay_step = float(cfg["PLASTICITY_DECAY_STEP"])
    neuron.n_excitatory_inputs = int(cfg["N_EXCITATORY_INPUTS"])
    neuron.n_inhibitory_inputs = int(cfg["N_INHIBITORY_INPUTS"])
    neuron.n_inputs = neuron.n_excitatory_inputs + neuron.n_inhibitory_inputs
    neuron.spike_threshold = float(cfg["SPIKE_THRESHOLD_MV"])
    neuron.minimum_membrane_potential = float(cfg["MINIMUM_MEMBRANE_POTENTIAL_MV"])
    neuron.resting_potential = float(cfg["RESTING_POTENTIAL_MV"])
    neuron.reset_potential_upper = float(cfg["RESET_POTENTIAL_UPPER_MV"])
    neuron.ampa_epsp_amplitude = float(cfg["AMPA_EPSP_AMPLITUDE_MV"])
    neuron.gabaa_ipsp_amplitude = float(cfg["GABAA_IPSP_AMPLITUDE_MV"])
    neuron.refractory_steps = int(cfg["REFRACTORY_STEPS"])
    neuron.minimum_synaptic_weight = float(cfg["MINIMUM_SYNAPTIC_WEIGHT"])
    neuron.nmda_activation_threshold = float(cfg["NMDA_ACTIVATION_THRESHOLD_MV"])
    neuron.plasticity_log_scale = float(cfg["PLASTICITY_LOG_SCALE"])
    neuron.inhibitory_decay_interval = int(cfg["INHIBITORY_DECAY_INTERVAL"])
    neuron.adaptation_reference_potential = float(
        cfg.get("ADAPTATION_REFERENCE_POTENTIAL_MV", neuron.adaptation_reference_potential)
    )
    neuron.nmda_psp_scale = float(cfg["NMDA_PSP_SCALE"])
    neuron.reset_state()
    return neuron
