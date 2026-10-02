#!/usr/bin/env python3
"""Reproduce the final manuscript analysis pipeline from the repository root.

Default mode re-runs the final publication analysis from the version-controlled
paired n=20 reference datasets. The canonical analysis uses 500 paired-seed
bootstrap resamples with seed 20261001, as implemented in
``scripts/final_publication_analysis.py``.

Use ``--simulate`` to regenerate the paired healthy/input-driven/NMDA-driven
concentration-response datasets before analysis. Use ``--sensitivity`` to rerun
the final OAT sensitivity analysis at the final healthy-state functional EC50
values and regenerate Supplementary Figure S1 and Table S3.
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
    p = argparse.ArgumentParser()
    p.add_argument(
        "--simulate",
        action="store_true",
        help="Regenerate the paired n=20 concentration-response datasets (slow).",
    )
    p.add_argument(
        "--sensitivity",
        action="store_true",
        help="Rerun the final n=20 OAT sensitivity analysis (slow).",
    )
    args = p.parse_args()

    if args.simulate:
        run("scripts/compare_two_pathologies_drugs_fast.py")

    # Final manuscript analysis: state-specific 4PL fits, 500 paired-seed
    # bootstrap resamples, healthy-preserving window, Tables 1-5 and Figures 1-5.
    run("scripts/final_publication_analysis.py")

    if args.sensitivity:
        run("scripts/run_sensitivity_oat_final_fast.py")

    print("\nCanonical manuscript outputs: results/publication_final/")
    if not args.simulate:
        print("Input datasets were read from results/reference/.")
    if not args.sensitivity:
        print("Sensitivity outputs were not regenerated; add --sensitivity to rerun them.")
    print("Use --simulate --sensitivity for full regeneration from the model.")


if __name__ == "__main__":
    main()
