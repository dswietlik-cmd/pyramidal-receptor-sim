import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_final_publication_outputs_exist():
    out = ROOT / "results" / "publication_final"
    en = out / "figures_professional" / "en"
    pl = out / "figures_professional" / "pl"

    assert (out / "Table1_baseline_states.csv").exists()
    assert (out / "Table2_state_specific_4PL_fits.csv").exists()
    assert (en / "Figure5_therapeutic_window_professional.png").exists()
    assert (en / "Figure5_therapeutic_window_professional.svg").exists()
    assert (en / "Figure5_therapeutic_window_professional.pdf").exists()
    assert (pl / "Figure5_therapeutic_window_professional.png").exists()


def test_ec50_is_not_reported_when_half_suppression_not_observed():
    path = ROOT / "results" / "publication_final" / "Table2_state_specific_4PL_fits.csv"
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    for row in rows:
        max_suppression = float(row["max_observed_suppression"])
        if max_suppression < 0.50:
            ec50 = row["functional_EC50_nM"].strip().lower()
            assert ec50 in {"", "nan"}
