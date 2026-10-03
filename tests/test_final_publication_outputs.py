import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_final_publication_outputs_exist():
    out = ROOT / "results" / "publication_final"
    en = out / "figures_professional" / "en"
    pl = out / "figures_professional" / "pl"

    required_tables = [
        "Table1_baseline_states.csv",
        "Table2_state_specific_4PL_fits.csv",
        "Table3_functional_target_concentrations.csv",
        "Table4_therapeutic_window_grid.csv",
        "Table5_best_effect_with_healthy_preservation_ge80.csv",
        "Table9_baseline_inferential_statistics.csv",
        "Table10_concentration_global_statistics.csv",
        "Table11_pathology_state_contrast_summary.csv",
        "TableS3_sensitivity_oat_final_ec50.csv",
        "TableS4_concentration_vs_zero_wilcoxon.csv",
        "TableS5_pathology_state_suppression_contrasts.csv",
    ]
    for name in required_tables:
        assert (out / name).exists(), name

    figure_stems = [
        "Figure1_study_design_professional",
        "Figure2_perampanel_professional",
        "Figure3_memantine_professional",
        "Figure4_diazepam_professional",
        "Figure5_therapeutic_window_professional",
        "FigureS1_sensitivity_professional",
    ]
    for directory in (en, pl):
        for stem in figure_stems:
            for extension in ("png", "svg", "pdf"):
                path = directory / f"{stem}.{extension}"
                assert path.exists(), path


def test_ec50_is_not_reported_when_half_suppression_not_observed():
    path = ROOT / "results" / "publication_final" / "Table2_state_specific_4PL_fits.csv"
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    for row in rows:
        max_suppression = float(row["max_observed_suppression"])
        if max_suppression < 0.50:
            ec50 = row["functional_EC50_nM"].strip().lower()
            assert ec50 in {"", "nan"}
