#!/usr/bin/env python3
"""Inferential statistics for the final paired n=20 manuscript experiment.

The analysis preserves the paired stochastic design by matching seed IDs across
states and concentrations. It writes publication-ready CSV tables to
``results/publication_final`` by default.

Important: the n=20 units are computational stochastic replicates, not biological
samples. P-values quantify consistency across the specified simulation ensemble and
must not be interpreted as biological population inference.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import friedmanchisquare, wilcoxon

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REFERENCE = ROOT / "results" / "reference"
DEFAULT_OUTPUT = ROOT / "results" / "publication_final"

STATES = ["healthy", "input_hyperexcitable", "nmda_excitotoxic"]
DRUGS = ["perampanel", "memantine", "diazepam"]


def holm_adjust(pvalues: list[float] | np.ndarray) -> np.ndarray:
    """Holm step-down family-wise error correction."""
    p = np.asarray(pvalues, dtype=float)
    m = len(p)
    if m == 0:
        return np.array([], dtype=float)
    order = np.argsort(p)
    adjusted = np.empty(m, dtype=float)
    running = 0.0
    for rank, idx in enumerate(order):
        value = (m - rank) * p[idx]
        value = max(value, running)
        running = value
        adjusted[idx] = min(value, 1.0)
    return adjusted


def matched_rank_biserial(x: np.ndarray, y: np.ndarray) -> float:
    """Matched-pairs rank-biserial effect size for paired samples."""
    d = np.asarray(x, dtype=float) - np.asarray(y, dtype=float)
    d = d[d != 0]
    if len(d) == 0:
        return 0.0
    ranks = pd.Series(np.abs(d)).rank(method="average").to_numpy()
    w_plus = ranks[d > 0].sum()
    w_minus = ranks[d < 0].sum()
    denom = len(d) * (len(d) + 1) / 2
    return float((w_plus - w_minus) / denom)


def paired_wilcoxon(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Return Wilcoxon statistic and p-value, handling all-zero differences."""
    d = np.asarray(x, dtype=float) - np.asarray(y, dtype=float)
    if np.allclose(d, 0.0):
        return 0.0, 1.0
    result = wilcoxon(
        x,
        y,
        alternative="two-sided",
        zero_method="wilcox",
        method="auto",
    )
    return float(result.statistic), float(result.pvalue)


def baseline_statistics(df: pd.DataFrame) -> pd.DataFrame:
    # Concentration-zero rows are duplicated across drugs; use one drug only.
    base = (
        df[(df.drug == "perampanel") & (df.concentration_nM == 0)]
        .pivot(index="seed", columns="state", values="firing_hz")
        .dropna()
    )
    arrays = [base[state].to_numpy(float) for state in STATES]
    statistic, p_value = friedmanchisquare(*arrays)
    kendall_w = statistic / (len(base) * (len(STATES) - 1))

    rows: list[dict] = [
        {
            "analysis": "overall",
            "comparison": "3 states",
            "n": len(base),
            "statistic": statistic,
            "df": len(STATES) - 1,
            "p_raw": p_value,
            "p_holm": np.nan,
            "effect_size": "Kendall W",
            "effect": kendall_w,
            "median_difference_hz": np.nan,
        }
    ]

    pairs = [
        ("healthy", "input_hyperexcitable"),
        ("healthy", "nmda_excitotoxic"),
        ("input_hyperexcitable", "nmda_excitotoxic"),
    ]
    temporary = []
    pvalues = []
    for a, b in pairs:
        x = base[a].to_numpy(float)
        y = base[b].to_numpy(float)
        stat, p = paired_wilcoxon(x, y)
        pvalues.append(p)
        temporary.append(
            (
                a,
                b,
                stat,
                p,
                matched_rank_biserial(x, y),
                float(np.median(x - y)),
            )
        )

    adjusted = holm_adjust(pvalues)
    for (a, b, stat, p, effect, median_delta), p_holm in zip(temporary, adjusted):
        rows.append(
            {
                "analysis": "posthoc",
                "comparison": f"{a} vs {b}",
                "n": len(base),
                "statistic": stat,
                "df": np.nan,
                "p_raw": p,
                "p_holm": p_holm,
                "effect_size": "matched rank-biserial r",
                "effect": effect,
                "median_difference_hz": median_delta,
            }
        )
    return pd.DataFrame(rows)


def concentration_statistics(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    global_rows: list[dict] = []
    posthoc_rows: list[dict] = []

    for drug in DRUGS:
        for state in STATES:
            group = df[(df.drug == drug) & (df.state == state)]
            pivot = (
                group.pivot(index="seed", columns="concentration_nM", values="firing_hz")
                .sort_index(axis=1)
                .dropna()
            )
            concentrations = [float(c) for c in pivot.columns]
            arrays = [pivot[c].to_numpy(float) for c in pivot.columns]
            statistic, p_value = friedmanchisquare(*arrays)
            kendall_w = statistic / (len(pivot) * (len(concentrations) - 1))

            nonzero = [c for c in concentrations if c > 0]
            temporary = []
            pvalues = []
            for concentration in nonzero:
                x = pivot[concentration].to_numpy(float)
                y = pivot[0.0].to_numpy(float)
                stat, p = paired_wilcoxon(x, y)
                pvalues.append(p)
                temporary.append(
                    (
                        concentration,
                        stat,
                        p,
                        matched_rank_biserial(x, y),
                        float(np.median(x - y)),
                        float(np.mean(x - y)),
                    )
                )

            adjusted = holm_adjust(pvalues)
            lowest_significant = None
            for (concentration, stat, p, effect, median_delta, mean_delta), p_holm in zip(
                temporary, adjusted
            ):
                if p_holm < 0.05 and lowest_significant is None:
                    lowest_significant = concentration
                posthoc_rows.append(
                    {
                        "drug": drug,
                        "state": state,
                        "concentration_nM": concentration,
                        "n": len(pivot),
                        "wilcoxon_W": stat,
                        "p_raw": p,
                        "p_holm": p_holm,
                        "rank_biserial_r": effect,
                        "median_delta_hz": median_delta,
                        "mean_delta_hz": mean_delta,
                    }
                )

            global_rows.append(
                {
                    "drug": drug,
                    "state": state,
                    "n": len(pivot),
                    "k_concentrations": len(concentrations),
                    "friedman_chi2": statistic,
                    "df": len(concentrations) - 1,
                    "p": p_value,
                    "kendall_W": kendall_w,
                    "lowest_significant_vs_zero_nM": lowest_significant,
                }
            )

    return pd.DataFrame(global_rows), pd.DataFrame(posthoc_rows)


def pathology_contrasts(pm: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    contrast_rows: list[dict] = []
    summary_rows: list[dict] = []

    for drug in DRUGS:
        input_state = pm[
            (pm.drug == drug) & (pm.pathology_state == "input_hyperexcitable")
        ]
        nmda_state = pm[(pm.drug == drug) & (pm.pathology_state == "nmda_excitotoxic")]
        merged = input_state.merge(
            nmda_state,
            on=["drug", "concentration_nM", "seed"],
            suffixes=("_input", "_nmda"),
        )
        concentrations = sorted(float(c) for c in merged.concentration_nM.unique() if c > 0)

        temporary = []
        pvalues = []
        for concentration in concentrations:
            group = merged[merged.concentration_nM == concentration]
            x = group.pathological_suppression_fraction_nmda.to_numpy(float)
            y = group.pathological_suppression_fraction_input.to_numpy(float)
            stat, p = paired_wilcoxon(x, y)
            pvalues.append(p)
            temporary.append(
                (
                    concentration,
                    len(group),
                    stat,
                    p,
                    matched_rank_biserial(x, y),
                    float(np.median(x - y)),
                    float(np.mean(x - y)),
                )
            )

        adjusted = holm_adjust(pvalues)
        first_significant = None
        drug_rows = []
        for (concentration, n, stat, p, effect, median_delta, mean_delta), p_holm in zip(
            temporary, adjusted
        ):
            if p_holm < 0.05 and first_significant is None:
                first_significant = concentration
            row = {
                "drug": drug,
                "concentration_nM": concentration,
                "n": n,
                "wilcoxon_W": stat,
                "p_raw": p,
                "p_holm": p_holm,
                "rank_biserial_r": effect,
                "median_delta_suppression_nmda_minus_input": median_delta,
                "mean_delta_suppression_nmda_minus_input": mean_delta,
            }
            contrast_rows.append(row)
            drug_rows.append(row)

        highest = max(drug_rows, key=lambda row: row["concentration_nM"])
        summary_rows.append(
            {
                "drug": drug,
                "first_significant_state_difference_nM": first_significant,
                "highest_concentration_nM": highest["concentration_nM"],
                "highest_mean_delta_percent": 100
                * highest["mean_delta_suppression_nmda_minus_input"],
                "highest_p_holm": highest["p_holm"],
                "highest_rank_biserial_r": highest["rank_biserial_r"],
            }
        )

    return pd.DataFrame(contrast_rows), pd.DataFrame(summary_rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-dir", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.reference_dir / "healthy_two_pathologies_drug_replicates.csv")
    pm = pd.read_csv(args.reference_dir / "healthy_two_pathologies_paired_metrics.csv")

    baseline = baseline_statistics(df)
    global_stats, concentration_posthoc = concentration_statistics(df)
    pathology_posthoc, pathology_summary = pathology_contrasts(pm)

    outputs = {
        "Table9_baseline_inferential_statistics.csv": baseline,
        "Table10_concentration_global_statistics.csv": global_stats,
        "Table11_pathology_state_contrast_summary.csv": pathology_summary,
        "TableS4_concentration_vs_zero_wilcoxon.csv": concentration_posthoc,
        "TableS5_pathology_state_suppression_contrasts.csv": pathology_posthoc,
    }
    for filename, frame in outputs.items():
        path = args.output_dir / filename
        frame.to_csv(path, index=False)
        print(f"wrote {path.relative_to(ROOT) if path.is_relative_to(ROOT) else path}")


if __name__ == "__main__":
    main()
