"""Pharmacodynamic mappings used by the receptor simulation model.

The values in this module describe concentration-response abstractions at the
receptor level. They are not pharmacokinetic models and must not be interpreted
as clinical doses or predicted human brain concentrations.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping


@dataclass(frozen=True)
class DrugModel:
    """Parameters for a receptor-level concentration-response model."""

    name: str
    target: str
    mechanism: str
    ec50_uM: float
    hill_coefficient: float
    maximal_response_ratio: float | None = None


DRUG_MODELS: Mapping[str, DrugModel] = {
    "memantine": DrugModel(
        name="Memantine",
        target="NMDA",
        mechanism="antagonist",
        ec50_uM=0.79,
        hill_coefficient=0.92,
    ),
    "perampanel": DrugModel(
        name="Perampanel",
        target="AMPA",
        mechanism="antagonist",
        ec50_uM=0.56,
        hill_coefficient=0.80,
    ),
    "diazepam": DrugModel(
        name="Diazepam",
        target="GABAA",
        mechanism="positive_allosteric_modulator",
        ec50_uM=0.025,
        hill_coefficient=1.0,
        maximal_response_ratio=6.1,
    ),
}


def concentration_to_uM(concentration: float, unit: str = "uM") -> float:
    """Convert a non-negative concentration to micromolar."""
    if concentration < 0:
        raise ValueError("concentration must be >= 0")

    normalized = unit.strip().replace("µ", "u").lower()
    factors = {
        "m": 1e6,
        "mm": 1e3,
        "um": 1.0,
        "nm": 1e-3,
        "pm": 1e-6,
    }
    if normalized not in factors:
        raise ValueError("unit must be one of M, mM, uM, nM, pM")
    return float(concentration) * factors[normalized]


def hill_fraction(concentration: float, ec50: float, hill_coefficient: float = 1.0) -> float:
    """Return fractional receptor effect in the interval [0, 1]."""
    if concentration < 0:
        raise ValueError("concentration must be >= 0")
    if ec50 <= 0:
        raise ValueError("ec50 must be > 0")
    if hill_coefficient <= 0:
        raise ValueError("hill_coefficient must be > 0")
    if concentration == 0:
        return 0.0

    c = concentration ** hill_coefficient
    e = ec50 ** hill_coefficient
    return c / (c + e)


def receptor_gain_for_drug(drug: str, concentration: float, unit: str = "uM") -> dict:
    """Map drug concentration to receptor gain and return an audit record."""
    key = drug.strip().lower()
    if key not in DRUG_MODELS:
        supported = ", ".join(sorted(DRUG_MODELS))
        raise ValueError(f"unknown drug {drug!r}; supported: {supported}")

    model = DRUG_MODELS[key]
    c_uM = concentration_to_uM(concentration, unit)
    fraction = hill_fraction(c_uM, model.ec50_uM, model.hill_coefficient)

    if model.mechanism == "antagonist":
        receptor_gain = 1.0 - fraction
        target_effect = fraction
    elif model.mechanism == "positive_allosteric_modulator":
        assert model.maximal_response_ratio is not None
        maximal_increase = model.maximal_response_ratio - 1.0
        receptor_gain = 1.0 + maximal_increase * fraction
        target_effect = receptor_gain - 1.0
    else:
        raise RuntimeError(f"unsupported mechanism: {model.mechanism}")

    return {
        "drug": model.name,
        "target": model.target,
        "mechanism": model.mechanism,
        "input_concentration": float(concentration),
        "input_unit": unit,
        "concentration_uM": c_uM,
        "hill_fraction": fraction,
        "target_effect": target_effect,
        "receptor_gain": receptor_gain,
        "ec50_or_ic50_uM": model.ec50_uM,
        "hill_coefficient": model.hill_coefficient,
    }


def apply_drug(neuron, drug: str, concentration: float, unit: str = "nM") -> dict:
    """Apply a supported drug concentration to an existing neuron instance."""
    return neuron.apply_drug_concentration(drug, concentration, unit)
