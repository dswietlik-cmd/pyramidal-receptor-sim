"""Deterministic and jittered synaptic input generators."""

from __future__ import annotations

import random
from typing import Iterable

import numpy as np


def period_steps(freq_hz: float, dt_ms: float) -> int:
    if freq_hz <= 0:
        return 0
    return max(int(round(1000.0 / (freq_hz * dt_ms))), 1)


def deterministic_inputs(step: int, frequencies: Iterable[float], ne: int, ni: int, dt_ms: float):
    ex = [0] * ne
    inh = [0] * ni
    for idx, freq in enumerate(frequencies):
        pulse = 1 if freq > 0 and step % period_steps(float(freq), dt_ms) == 0 else 0
        if idx < ne:
            ex[idx] = pulse
        else:
            inh[idx - ne] = pulse
    return ex, inh


def jittered_stream(
    frequencies: Iterable[float],
    steps: int,
    dt_ms: float,
    ne: int,
    ni: int,
    seed: int,
    jitter_ms: float = 0.5,
):
    """Generate rate-preserving periodic input trains with bounded temporal jitter."""
    freqs = list(map(float, frequencies))
    rng = random.Random(seed)
    jitter_steps = int(round(jitter_ms / dt_ms))
    ex = np.zeros((steps, ne), dtype=np.uint8)
    inh = np.zeros((steps, ni), dtype=np.uint8)
    for idx, freq in enumerate(freqs):
        if freq <= 0:
            continue
        period = period_steps(freq, dt_ms)
        for base in range(0, steps, period):
            delta = rng.randint(-jitter_steps, jitter_steps) if jitter_steps else 0
            t = max(0, min(steps - 1, base + delta))
            if idx < ne:
                ex[t, idx] = 1
            else:
                inh[t, idx - ne] = 1
    return ex, inh
