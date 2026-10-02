#!/usr/bin/env python3
"""Reproduce the final manuscript analysis pipeline from the repository root.

Default mode re-analyzes the version-controlled paired n=20 reference datasets and
regenerates the final 4PL analysis, inferential statistics, and professional figures.

Use ``--simulate`` to regenerate the paired healthy/input-driven/NMDA-driven
concentration-response datasets before analysis. Use ``--sensitivity`` to rerun the
final OAT sensitivity analysis at the final healthy-state functional EC50 values before
regenerating the professional Supplementary Figure S1.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def run(*args: str) -> None:
    cmd = [sys.executable, *args]
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=ROOT, check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--simulate",
        action="store_true",
        help="Regenerate the paired n=20 concentration-response datasets (slow).",
    )
    parser.add_argument(
        "--sensitivity",
        action="store_true",
        help="Rerun the final n=20 OAT sensitivity analysis (slow).",
    )
    parser.add_argument(
        "--tiff",
        action="store_true",
        help="Also export large 800-dpi TIFF versions of the professional figures.",
    )
    args = parser.parse_args()

    if args.simulate:
        run("scripts/compare_two_pathologies_drugs_fast.py")

    # Canonical main analysis: state-specific 4PL fits, 500 paired-seed bootstrap
    # resamples, healthy-preserving window, Tables 1-5 and legacy Figures 1-5.
    run("scripts/final_publication_analysis.py")

    # Inferential analysis of the paired computational replicates.
    run("scripts/run_stat_analysis.py")

    if args.sensitivity:
        run("scripts/run_sensitivity_oat_final_fast.py")

    # Journal-ready vector/PDF plus 800-dpi PNG figures. These use the committed
    # sensitivity table unless --sensitivity regenerated it in this run.
    figure_args = ["scripts/make_professional_figures.py"]
    if args.tiff:
        figure_args.append("--tiff")
    run(*figure_args)

    print("\nCanonical manuscript outputs: results/publication_final/")
    print("Inferential statistics: Tables 9-11 and Supplementary Tables S4-S5.")
    print("Professional figures: results/publication_final/figures_professional/")
    if not args.simulate:
        print("Input datasets were read from results/reference/.")
    if not args.sensitivity:
        print("Sensitivity Table S3 was not regenerated; add --sensitivity to rerun it.")
    if not args.tiff:
        print("TIFF export was skipped; add --tiff for 800-dpi TIFF files.")


if __name__ == "__main__":
    main()
