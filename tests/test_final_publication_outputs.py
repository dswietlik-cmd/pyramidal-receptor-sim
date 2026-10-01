from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

def test_final_publication_outputs_exist():
    out = ROOT / "results" / "publication_final"
    assert (out / "Table1_baseline_states.csv").exists()
    assert (out / "Table2_state_specific_4PL_fits.csv").exists()
    assert (out / "Figure5_healthy_preserving_therapeutic_window.png").exists()

def test_ec50_is_not_reported_when_half_suppression_not_observed():
    df = pd.read_csv(ROOT / "results" / "publication_final" / "Table2_state_specific_4PL_fits.csv")
    bad = df[df["max_observed_suppression"] < 0.50]
    assert bad["functional_EC50_nM"].isna().all()
