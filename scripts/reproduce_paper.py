#!/usr/bin/env python3
"""Reproduce the manuscript analysis pipeline from the repository root.

Default mode re-analyzes the version-controlled n=20 reference replicate dataset quickly.
Use --simulate to regenerate the 20x concentration-response simulations before fitting.
Use --sensitivity to additionally rerun the computationally expensive OAT analysis.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "results" / "reference"


def run(*args: str) -> None:
    cmd = [sys.executable, *args]
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=ROOT, check=True)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--simulate", action="store_true", help="Regenerate n=20 simulations (slow).")
    p.add_argument("--sensitivity", action="store_true", help="Rerun n=20 OAT sensitivity (slow).")
    p.add_argument("--bootstrap", type=int, default=200)
    p.add_argument("--outdir", default="results/reproduced")
    args = p.parse_args()

    out = ROOT / args.outdir
    out.mkdir(parents=True, exist_ok=True)

    replicates = out / "publication_n20_replicates.csv"
    if args.simulate:
        run("scripts/run_publication_n20.py")
        shutil.move(str(ROOT / "publication_n20_replicates.csv"), replicates)
    else:
        shutil.copy2(REF / "publication_n20_replicates.csv", replicates)

    run("scripts/fit_4pl.py", str(replicates), "--output", str(out / "publication_n20_fits.csv"))
    run(
        "scripts/bootstrap_4pl.py", str(replicates),
        "--bootstrap", str(args.bootstrap),
        "--fits-output", str(out / "publication_n20_bootstrap_fits.csv"),
        "--targets-output", str(out / "publication_n20_bootstrap_targets.csv"),
    )

    if args.sensitivity:
        run(
            "scripts/run_sensitivity_oat.py",
            "--fits", str(out / "publication_n20_fits.csv"),
            "--output", str(out / "publication_sensitivity_oat_n20.csv"),
        )
    else:
        shutil.copy2(REF / "publication_sensitivity_oat_n20.csv", out / "publication_sensitivity_oat_n20.csv")

    print(f"\nReproduction outputs: {out.relative_to(ROOT)}")
    print("Use --simulate and --sensitivity for full regeneration from the model.")


if __name__ == "__main__":
    main()
