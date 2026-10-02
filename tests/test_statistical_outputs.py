import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "publication_final"


def read_rows(name):
    with (OUT / name).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def test_inferential_statistical_outputs_exist():
    expected = [
        "Table9_baseline_inferential_statistics.csv",
        "Table10_concentration_global_statistics.csv",
        "Table11_pathology_state_contrast_summary.csv",
        "TableS4_concentration_vs_zero_wilcoxon.csv",
        "TableS5_pathology_state_suppression_contrasts.csv",
    ]
    for name in expected:
        assert (OUT / name).exists(), name


def test_baseline_statistics_have_paired_n20_design():
    rows = read_rows("Table9_baseline_inferential_statistics.csv")
    assert len(rows) == 4
    overall = rows[0]
    assert overall["analysis"] == "overall"
    assert int(overall["n"]) == 20
    assert overall["effect_size"] == "Kendall W"
    assert 0.0 <= float(overall["effect"]) <= 1.0


def test_concentration_global_statistics_cover_all_drug_state_pairs():
    rows = read_rows("Table10_concentration_global_statistics.csv")
    assert len(rows) == 9
    pairs = {(row["drug"], row["state"]) for row in rows}
    assert len(pairs) == 9
    assert {int(row["n"]) for row in rows} == {20}


def test_pathology_contrast_summary_contains_three_drugs():
    rows = read_rows("Table11_pathology_state_contrast_summary.csv")
    assert {row["drug"] for row in rows} == {"perampanel", "memantine", "diazepam"}
