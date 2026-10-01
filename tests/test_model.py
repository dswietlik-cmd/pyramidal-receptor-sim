from pyramidal_receptor_sim import Neuron


def test_zero_concentration_keeps_baseline_gain():
    n = Neuron()
    n.apply_memantine(0, "nM")
    assert n.nmda_gain == 1.0

    n = Neuron()
    n.apply_perampanel(0, "nM")
    assert n.ampa_gain == 1.0


def test_diazepam_increases_gabaa_gain():
    n = Neuron()
    n.apply_diazepam(10, "nM")
    assert n.gabaa_gain > 1.0


def test_antagonists_reduce_target_gain():
    n = Neuron()
    n.apply_memantine(150, "nM")
    assert 0 < n.nmda_gain < 1

    n = Neuron()
    n.apply_perampanel(12, "nM")
    assert 0 < n.ampa_gain < 1


def test_public_parameter_names_are_descriptive():
    n = Neuron()
    assert hasattr(n, "spike_threshold")
    assert hasattr(n, "nmda_activation_threshold")
    assert hasattr(n, "ampa_epsp_amplitude")
    assert hasattr(n, "gabaa_ipsp_amplitude")
