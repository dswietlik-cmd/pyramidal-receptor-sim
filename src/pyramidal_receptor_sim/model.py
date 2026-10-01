import math
from dataclasses import dataclass
from typing import List

from .drugs import receptor_gain_for_drug


@dataclass
class InterspikePhase:
    current_interval_steps: int = 0
    previous_interval_steps: int = 0


class Neuron:
    """CA1-like pyramidal-neuron model with AMPA, NMDA and GABA-A components.

    The model uses discrete 0.5 ms simulation steps in the reference
    configuration. Excitatory inputs contain separate AMPA and NMDA response
    registers; inhibitory inputs represent GABA-A-mediated transmission. NMDA
    activation is voltage gated and contributes both to the somatic
    postsynaptic potential and to the model's synaptic-plasticity state.

    Receptor gains are dimensionless multipliers. A value of 1.0 represents
    the calibrated drug-free reference state.
    """

    def __init__(self, all_values: bool = True):
        self.spike_steps: List[int] = []
        self.interval_phase = InterspikePhase()
        self.init(all_values)

    def init_values(self):
        self.plasticity_decay_step = 1
        self.n_excitatory_inputs = 13
        self.n_inhibitory_inputs = 3
        self.n_inputs = self.n_excitatory_inputs + self.n_inhibitory_inputs
        self.spike_threshold = -50
        self.minimum_membrane_potential = -90
        self.resting_potential = -80
        self.reset_potential_upper = -79
        # Externally anchored effective clustered synaptic amplitudes for a CA1-like model.
        # AMPA/EPSP: 5 mV lies within experimentally used Schaffer-collateral evoked
        # CA1 pyramidal EPSPs (~4-7 mV) and close to a reported mean ~4.24 mV.
        self.ampa_epsp_amplitude = 5.0
        # GABA-A/IPSP amplitude is calibrated to the experimental compound E/I ratio
        # in CA1 pyramidal neurons: EPSP 2.1 mV / IPSP 0.93 mV = 2.258.
        # With a 5 mV AMPA reference amplitude, the calibrated GABA-A amplitude is -2.2143 mV.
        self.gabaa_ipsp_amplitude = -2.2142857142857144
        self.refractory_steps = 3
        self.minimum_synaptic_weight = 0.1
        self.nmda_activation_threshold = -68
        self.plasticity_log_scale = 2.3026
        self.inhibitory_decay_interval = 15
        self.adaptation_reference_potential = -10

        # Pharmacological receptor modulation.
        # 1.0 = baseline receptor function.
        # AMPA and NMDA antagonism use values in [0, 1].
        # GABA-A positive modulation can use values > 1.
        self.ampa_gain = 1.0
        self.nmda_gain = 1.0
        self.gabaa_gain = 1.0

        # Explicit NMDA contribution to the postsynaptic membrane potential.
        # Externally calibrated for a CA1 Schaffer-collateral-like synapse.
        # Otmakhova et al. reported an NMDA/AMPA EPSP-area ratio of 1.46 +/- 0.3
        # for Schaffer collateral input in low-Mg2+ ACSF under current clamp.
        # The native model kernel-area ratio is ~0.2522, so scale = 1.46/0.2522
        # = ~5.79. Voltage-dependent CaPath gating remains active in normal runs.
        self.nmda_psp_scale = 5.789

        # Per-step decomposition of somatic PSP for audit/analysis.
        self.last_ampa_psp = 0.0
        self.last_nmda_psp = 0.0
        self.last_gaba_psp = 0.0

    def init_tables(self):
        self.ampa_kernel = [0.0] * 33
        self.nmda_kernel = [0.0] * 33
        self.gabaa_kernel = [0.0] * 33

        self.plasticity_state = [0.0] * (self.n_excitatory_inputs + self.n_inhibitory_inputs + 1)
        self.spine_potential = [0.0] * (self.n_excitatory_inputs + self.n_inhibitory_inputs + 1)
        self.excitatory_input_state = [0] * (self.n_excitatory_inputs + 1)
        self.inhibitory_input_state = [0] * (self.n_inhibitory_inputs + 1)
        self.ampa_registers = [[0.0] * 33 for _ in range(self.n_excitatory_inputs + 1)]
        self.nmda_registers = [[0.0] * 33 for _ in range(self.n_excitatory_inputs + 1)]
        self.gabaa_registers = [[0.0] * 33 for _ in range(self.n_inhibitory_inputs + 1)]
        self.inhibitory_delay_counter = [0] * (self.n_inhibitory_inputs + 1)
        self.inhibitory_decay_counter = [0] * (self.n_inhibitory_inputs + 1)
        self.nmda_gate_open = [0] * (self.n_excitatory_inputs + 1)

    def fill_vectors(self):
        for i in range(1, self.n_excitatory_inputs + 1):
            self.excitatory_input_state[i] = 0
        for i in range(1, self.n_inhibitory_inputs + 1):
            self.inhibitory_input_state[i] = 0

        for i in range(1, self.n_inhibitory_inputs + 1):
            self.inhibitory_decay_counter[i] = 0
        self.inhibitory_slow_decay_state = 0

        self.ampa_kernel[0] = self.ampa_kernel[1] = self.ampa_kernel[2] = self.ampa_kernel[32] = 0.0
        for i in range(1, 5):
            self.ampa_kernel[2 + i] = self.ampa_epsp_amplitude / 4 * i
        for i in range(1, 26):
            self.ampa_kernel[6 + i] = self.ampa_epsp_amplitude / 26 * (26 - i)

        self.nmda_kernel[0] = self.nmda_kernel[1] = self.nmda_kernel[2] = self.nmda_kernel[32] = 0.0
        for i in range(1, 5):
            self.nmda_kernel[2 + i] = self.ampa_epsp_amplitude / 20 * i
        for i in range(1, 26):
            self.nmda_kernel[6 + i] = self.ampa_epsp_amplitude / 99 * (26 - i)

        self.gabaa_kernel[0] = self.gabaa_kernel[1] = self.gabaa_kernel[2] = 0.0
        for i in range(1, 5):
            self.gabaa_kernel[2 + i] = self.gabaa_ipsp_amplitude / 4 * i
        for i in range(1, 20):
            self.gabaa_kernel[6 + i] = self.gabaa_ipsp_amplitude / 20 * (20 - i)
        for i in range(1, 8):
            self.gabaa_kernel[25 + i] = 0.0

        for idx in range(33):
            for k in range(1, self.n_excitatory_inputs + 1):
                self.ampa_registers[k][idx] = self.resting_potential
                self.nmda_registers[k][idx] = self.resting_potential
            for k in range(1, self.n_inhibitory_inputs + 1):
                self.gabaa_registers[k][idx] = self.resting_potential

        # Dynamic arrays are initialized explicitly so repeated experiments are reproducible.

    def init(self, all_values: bool = True):
        if all_values:
            self.init_values()
        self.init_tables()
        self.fill_vectors()
        self.refractory_counter = 0
        self.interspike_counter = 0
        self.spike_output = 0
        self.membrane_potential = float(self.resting_potential)
        self.spike_steps = []
        self.interval_phase = InterspikePhase()

    def reset_state(self):
        """Reset all dynamic state while preserving current model parameters."""
        # Rebuild dynamic registers using the current input counts and parameters.
        self.init_tables()
        self.fill_vectors()
        for i in range(1, self.n_excitatory_inputs + 1):
            self.plasticity_state[i] = 0.0
            self.nmda_gate_open[i] = 0
            self.spine_potential[i] = float(self.resting_potential)
        for i in range(self.n_excitatory_inputs + 1, self.n_inputs + 1):
            self.spine_potential[i] = float(self.resting_potential)
        self.refractory_counter = 0
        self.interspike_counter = 0
        self.spike_output = 0
        self.membrane_potential = float(self.resting_potential)
        self.inhibitory_slow_decay_state = 0
        self.spike_steps = []
        self.interval_phase = InterspikePhase()
        self.last_ampa_psp = 0.0
        self.last_nmda_psp = 0.0
        self.last_gaba_psp = 0.0

    def influence(self, source: int, dest: int) -> float:
        if dest >= source:
            return 1 - (1 - self.minimum_synaptic_weight) / (self.n_inputs - 1) * (dest - source)
        if (source - dest) <= (self.n_inputs / 4):
            return 1 - 4 / self.n_inputs * (source - dest)
        return 0.0

    def update_spine_potential(self, number: int) -> int:
        if not (1 <= number <= self.n_inputs):
            return 1
        pot = float(self.resting_potential)
        for m in range(1, self.n_excitatory_inputs + 1):
            pot += (self.ampa_registers[m][0] - self.resting_potential) * self.influence(m, number)
        for m in range(1, self.n_inhibitory_inputs + 1):
            pot += (self.gabaa_registers[m][0] - self.resting_potential) * self.influence(m + self.n_excitatory_inputs, number)
        self.spine_potential[number] = pot
        return 0

    def nmda_threshold_crossed(self, number: int) -> int:
        self.update_spine_potential(number)
        return 1 if self.spine_potential[number] >= self.nmda_activation_threshold else 0

    def update_nmda_memory(self, number: int):
        power = (self.nmda_registers[number][0] - self.resting_potential) * 6
        self.plasticity_state[number] += math.exp(power) - 1
        if self.plasticity_state[number] > self.plasticity_decay_step:
            self.plasticity_state[number] -= self.plasticity_decay_step
        else:
            self.plasticity_state[number] = 0.0

    def synaptic_weight(self, number: int) -> float:
        return (1 - self.minimum_synaptic_weight) / (self.n_excitatory_inputs - 1) * (number - self.n_excitatory_inputs) + 1

    def plasticity_gain(self, number: int) -> float:
        return 1 + (1 / 6) * math.log(self.plasticity_state[number] + 1) / self.plasticity_log_scale

    def post_spike_reset_potential(self, number: int) -> float:
        return self.resting_potential - ((self.resting_potential - self.reset_potential_upper) / self.n_excitatory_inputs) * (self.n_excitatory_inputs - number)

    def adaptation_factor(self, i: int) -> float:
        return ((self.resting_potential - self.adaptation_reference_potential) - (self.resting_potential - self.spine_potential[i])) / (self.resting_potential - self.adaptation_reference_potential)

    def set_receptor_modulation(
        self,
        ampa_gain: float = 1.0,
        nmda_gain: float = 1.0,
        gaba_gain: float = 1.0,
        nmda_psp_scale: float | None = None,
    ):
        """Set receptor-level functional gains without resetting dynamic state."""
        if ampa_gain < 0 or nmda_gain < 0 or gaba_gain < 0:
            raise ValueError("receptor gains must be >= 0")
        if nmda_psp_scale is not None and nmda_psp_scale < 0:
            raise ValueError("nmda_psp_scale must be >= 0")

        self.ampa_gain = float(ampa_gain)
        self.nmda_gain = float(nmda_gain)
        self.gabaa_gain = float(gaba_gain)
        if nmda_psp_scale is not None:
            self.nmda_psp_scale = float(nmda_psp_scale)

    def reset_receptor_gains(self):
        """Return AMPA, NMDA and GABA-A gains to the drug-free reference state."""
        self.ampa_gain = 1.0
        self.nmda_gain = 1.0
        self.gabaa_gain = 1.0

    def apply_drug_concentration(self, drug: str, concentration: float, unit: str = "uM") -> dict:
        """Apply a supported drug using its receptor-level concentration-response model.

        This method is a phenomenological pharmacodynamic mapping. It does not
        model absorption, distribution, metabolism, elimination, clinical dose,
        or human brain exposure.
        """
        result = receptor_gain_for_drug(drug, concentration, unit)
        target = result["target"]
        gain = result["receptor_gain"]

        if target == "NMDA":
            self.nmda_gain = gain
        elif target == "AMPA":
            self.ampa_gain = gain
        elif target == "GABAA":
            self.gabaa_gain = gain
        else:
            raise RuntimeError(f"unsupported receptor target: {target}")
        return result

    def apply_memantine(self, concentration: float, unit: str = "uM") -> dict:
        return self.apply_drug_concentration("memantine", concentration, unit)

    def apply_perampanel(self, concentration: float, unit: str = "uM") -> dict:
        return self.apply_drug_concentration("perampanel", concentration, unit)

    def apply_diazepam(self, concentration: float, unit: str = "nM") -> dict:
        return self.apply_drug_concentration("diazepam", concentration, unit)

    def set_inputs(self, ex: List[int], inh: List[int]):
        if len(ex) != self.n_excitatory_inputs:
            raise ValueError(f"Expected {self.n_excitatory_inputs} excitatory inputs, got {len(ex)}")
        if len(inh) != self.n_inhibitory_inputs:
            raise ValueError(f"Expected {self.n_inhibitory_inputs} inhibitory inputs, got {len(inh)}")
        for i, value in enumerate(ex, start=1):
            self.excitatory_input_state[i] = 1 if value else 0
        for i, value in enumerate(inh, start=1):
            self.inhibitory_input_state[i] = 1 if value else 0

    @staticmethod
    def _shift_left(row: List[float]):
        for k in range(32):
            row[k] = row[k + 1]

    def run(self) -> int:
        self.spike_output = 0

        # Excitatory input update.
        for i in range(1, self.n_excitatory_inputs + 1):
            if self.excitatory_input_state[i] == 1:
                factor = self.adaptation_factor(i) * self.plasticity_gain(i)
                for k in range(33):
                    self.ampa_registers[i][k] += factor * self.ampa_gain * self.ampa_kernel[k]

        # Inhibitory input update.
        for i in range(1, self.n_inhibitory_inputs + 1):
            if self.inhibitory_input_state[i] == 1:
                self.inhibitory_delay_counter[i] = 6
                for k in range(33):
                    self.gabaa_registers[i][k] += self.gabaa_gain * self.gabaa_kernel[k]

        # NMDA/plasticity update: the previous voltage-gating state controls
        # insertion of the NMDA response before the gating state is recalculated.
        for i in range(1, self.n_excitatory_inputs + 1):
            if self.excitatory_input_state[i] == 1 and self.nmda_gate_open[i] == 1:
                for k in range(33):
                    self.nmda_registers[i][k] += self.nmda_gain * self.nmda_kernel[k]
            self.update_nmda_memory(i)
            self.nmda_gate_open[i] = self.nmda_threshold_crossed(i)

        # Calculate inhibitory/somatic compartment local potentials.
        for i in range(self.n_excitatory_inputs + 1, self.n_inputs + 1):
            self.update_spine_potential(i)

        if self.refractory_counter == 0:
            self.membrane_potential = float(self.resting_potential)
            self.last_ampa_psp = 0.0
            self.last_nmda_psp = 0.0
            self.last_gaba_psp = 0.0
            for i in range(1, self.n_excitatory_inputs + 1):
                ampa_term = self.synaptic_weight(i) * (self.ampa_registers[i][0] - self.resting_potential)
                # The NMDA register now contributes explicitly to the somatic PSP,
                # using the same spatial weight as the corresponding excitatory input.
                # The NMDA gain acts when the NMDA kernel is inserted; nmda_psp_scale
                # controls the calibrated coupling of that register to membrane potential.
                nmda_term = self.nmda_psp_scale * self.synaptic_weight(i) * (self.nmda_registers[i][0] - self.resting_potential)
                self.last_ampa_psp += ampa_term
                self.last_nmda_psp += nmda_term
                self.membrane_potential += ampa_term + nmda_term
            for i in range(1, self.n_inhibitory_inputs + 1):
                gaba_term = self.gabaa_registers[i][0] - self.resting_potential
                self.last_gaba_psp += gaba_term
                self.membrane_potential += gaba_term

            if self.membrane_potential < self.minimum_membrane_potential:
                self.membrane_potential = float(self.minimum_membrane_potential)

            if self.membrane_potential >= self.spike_threshold:
                self.spike_output = 1
                self.refractory_counter = self.refractory_steps
                for idx in range(33):
                    for k in range(1, self.n_excitatory_inputs + 1):
                        self.ampa_registers[k][idx] = self.post_spike_reset_potential(k)
                    for k in range(1, self.n_inhibitory_inputs + 1):
                        self.gabaa_registers[k][idx] = self.resting_potential
            else:
                for i in range(1, self.n_excitatory_inputs + 1):
                    self._shift_left(self.ampa_registers[i])
                    self._shift_left(self.nmda_registers[i])

                if self.inhibitory_slow_decay_state == 0:
                    for i in range(1, self.n_inhibitory_inputs + 1):
                        self._shift_left(self.gabaa_registers[i])
                else:
                    for i in range(1, self.n_inhibitory_inputs + 1):
                        if self.inhibitory_delay_counter[i] > 0:
                            self.inhibitory_delay_counter[i] -= 1
                            self._shift_left(self.gabaa_registers[i])
                        else:
                            self.inhibitory_decay_counter[i] += 1
                            if self.inhibitory_decay_counter[i] == self.inhibitory_decay_interval:
                                self.inhibitory_decay_counter[i] = 0
                                self._shift_left(self.gabaa_registers[i])
        else:
            # During the refractory period, synaptic responses must continue
            # to evolve in time. Otherwise excitatory registers would be artificially frozen,
            # allowing excitation to accumulate and causing longer refractory
            # periods to paradoxically increase firing rate.
            for i in range(1, self.n_excitatory_inputs + 1):
                self._shift_left(self.ampa_registers[i])
                self._shift_left(self.nmda_registers[i])

            # Inhibitory synaptic responses should also continue to decay/shift
            # according to the same sleep-state rules used outside refractoriness.
            if self.inhibitory_slow_decay_state == 0:
                for i in range(1, self.n_inhibitory_inputs + 1):
                    self._shift_left(self.gabaa_registers[i])
            else:
                for i in range(1, self.n_inhibitory_inputs + 1):
                    if self.inhibitory_delay_counter[i] > 0:
                        self.inhibitory_delay_counter[i] -= 1
                        self._shift_left(self.gabaa_registers[i])
                    else:
                        self.inhibitory_decay_counter[i] += 1
                        if self.inhibitory_decay_counter[i] == self.inhibitory_decay_interval:
                            self.inhibitory_decay_counter[i] = 0
                            self._shift_left(self.gabaa_registers[i])

            if self.refractory_counter == (self.refractory_steps - 1):
                self.membrane_potential = float(self.resting_potential)
            self.refractory_counter -= 1

        self.inhibitory_slow_decay_state = 1 if self.membrane_potential < self.resting_potential else 0

        if self.spike_output == 1:
            self.spike_steps.append(self.interspike_counter)
            k = len(self.spike_steps)
            if k > 2:
                current = self.spike_steps[k - 1]
                aprev_step = self.spike_steps[k - 2]
                pprev_step = self.spike_steps[k - 3]
                self.interval_phase.current_interval_steps = current - aprev_step
                self.interval_phase.previous_interval_steps = aprev_step - pprev_step
        else:
            self.interspike_counter += 1

        return self.spike_output

    def step(self, ex: List[int], inh: List[int]) -> int:
        self.set_inputs(ex, inh)
        return self.run()

    def debug_state(self) -> dict:
        return {
            "spike": self.spike_output,
            "psp": self.membrane_potential,
            "refrak": self.refractory_counter,
            "sleep_state": self.inhibitory_slow_decay_state,
            "mem_vol": self.plasticity_state[1:self.n_excitatory_inputs + 1],
            "memory_gain": [self.plasticity_gain(i) for i in range(1, self.n_excitatory_inputs + 1)],
            "ca_path": self.nmda_gate_open[1:self.n_excitatory_inputs + 1],
            "spine_pot": self.spine_potential[1:self.n_inputs + 1],
            "ampa_gain": self.ampa_gain,
            "nmda_gain": self.nmda_gain,
            "gaba_gain": self.gabaa_gain,
            "nmda_psp_scale": self.nmda_psp_scale,
            "ampa_psp_component": self.last_ampa_psp,
            "nmda_psp_component": self.last_nmda_psp,
            "gaba_psp_component": self.last_gaba_psp,
        }

