import csv
from pathlib import Path


def test_reference_results_present():
    root = Path(__file__).resolve().parents[1]
    required = [
        "publication_n20_replicates.csv",
        "publication_n20_summary.csv",
        "publication_n20_fits.csv",
        "publication_n20_bootstrap_fits.csv",
        "publication_n20_bootstrap_targets.csv",
        "publication_sensitivity_oat_n20.csv",
    ]
    for name in required:
        assert (root / "results" / "reference" / name).is_file()


def test_reference_contains_three_drugs():
    root = Path(__file__).resolve().parents[1]
    path = root / "results" / "reference" / "publication_n20_fits.csv"
    with path.open(newline="", encoding="utf-8") as f:
        drugs = {row["drug"] for row in csv.DictReader(f)}
    assert drugs == {"perampanel", "memantine", "diazepam"}
