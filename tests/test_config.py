from pyramidal_receptor_sim import Neuron, configure_neuron, load_config


def test_reference_config_loads():
    cfg = load_config("configs/ca1_reference.txt")
    n = configure_neuron(Neuron(), cfg)
    assert n.n_excitatory_inputs == 13
    assert n.n_inhibitory_inputs == 3
    assert n.nmda_psp_scale > 0
    assert n.spike_threshold == -52
