from pyramidal_receptor_sim.pathology import STATES, apply_excitatory_drive

def test_pathology_presets_are_ordered():
    assert STATES['healthy'].excitatory_drive_multiplier < STATES['mild'].excitatory_drive_multiplier < STATES['moderate'].excitatory_drive_multiplier < STATES['severe'].excitatory_drive_multiplier

def test_only_excitatory_frequencies_are_scaled():
    freqs=[10,20,30,5]
    out=apply_excitatory_drive(freqs,3,1.5)
    assert out == [15.0,30.0,45.0,5.0]
