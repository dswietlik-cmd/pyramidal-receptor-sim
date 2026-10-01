"""CA1-like receptor pharmacodynamics simulator."""

from .config import configure_neuron, load_config
from .drugs import DRUG_MODELS, DrugModel, apply_drug
from .model import Neuron

__all__ = [
    "Neuron",
    "DrugModel",
    "DRUG_MODELS",
    "apply_drug",
    "load_config",
    "configure_neuron",
]

__version__ = "0.9.3"

from .pathology import HyperexcitabilityState, STATES, HEALTHY, MILD, MODERATE, SEVERE
