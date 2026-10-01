from pyramidal_receptor_sim.pathology import STATES, apply_excitatory_drive

def test_pathology_presets_are_ordered():
    assert STATES['healthy'].excitatory_drive_multiplier < STATES['mild'].excitatory_drive_multiplier < STATES['moderate'].excitatory_drive_multiplier < STATES['severe'].excitatory_drive_multiplier

def test_only_excitatory_frequencies_are_scaled():
    freqs=[10,20,30,5]
    out=apply_excitatory_drive(freqs,3,1.5)
    assert out == [15.0,30.0,45.0,5.0]


def test_nmda_pathology_multiplier_is_separate_from_drug_gain():
    from pyramidal_receptor_sim import Neuron
    n = Neuron()
    n.set_nmda_pathology_multiplier(4.1)
    assert n.nmda_pathology_multiplier == 4.1
    assert n.nmda_gain == 1.0

def test_invalid_nmda_pathology_multiplier_rejected():
    import pytest
    from pyramidal_receptor_sim import Neuron
    n = Neuron()
    with pytest.raises(ValueError):
        n.set_nmda_pathology_multiplier(0)


def test_gabaa_pathology_multiplier_is_separate_from_drug_gain():
    from pyramidal_receptor_sim import Neuron
    n = Neuron()
    n.set_gabaa_pathology_multiplier(0.5)
    assert n.gabaa_pathology_multiplier == 0.5
    assert n.gabaa_gain == 1.0
    n.apply_diazepam(10.0, "nM")
    assert n.gabaa_pathology_multiplier == 0.5
    assert n.gabaa_gain > 1.0


def test_gabaa_pathology_multiplier_bounds():
    import pytest
    from pyramidal_receptor_sim import Neuron
    n = Neuron()
    with pytest.raises(ValueError):
        n.set_gabaa_pathology_multiplier(-0.01)
    with pytest.raises(ValueError):
        n.set_gabaa_pathology_multiplier(1.01)
