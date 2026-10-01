#!/usr/bin/env python3
"""Paired-seed bootstrap of 4PL concentration-response curves.

Resamples replicate/seed identities with replacement so the paired design is preserved
across concentrations for each drug. Produces bootstrap EC50/Hill/Emax summaries and
concentrations predicted to yield 25%, 50%, 75%, and 90% firing suppression.
"""
from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit


def four_pl(x, bottom, top, ec50, hill):
    x = np.asarray(x, dtype=float)
    return bottom + (top - bottom) * (x ** hill) / (ec50 ** hill + x ** hill)


def target_concentration(pars, target: float) -> float:
    bottom, top, ec50, hill = pars
    if not (bottom < target < top):
        return float("nan")
    ratio = (target - bottom) / (top - target)
    return ec50 * ratio ** (1.0 / hill)


def fit_curve(x, y):
    nonzero = x[x > 0]
    p0 = [0.0, min(1.0, float(np.max(y))), float(np.median(nonzero)), 2.0]
    pars, _ = curve_fit(
        four_pl, x, y, p0=p0,
        bounds=([0, .2, .001, .1], [.2, 1.05, 5000, 20]),
        maxfev=50000,
    )
    return pars


def q(v, p):
    a = np.asarray([x for x in v if np.isfinite(x)], dtype=float)
    return float(np.percentile(a, p)) if len(a) else float("nan")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("csv")
    p.add_argument("--bootstrap", type=int, default=200)
    p.add_argument("--seed", type=int, default=20261001)
    p.add_argument("--fits-output", default="publication_n20_bootstrap_fits.csv")
    p.add_argument("--targets-output", default="publication_n20_bootstrap_targets.csv")
    args = p.parse_args()

    with open(args.csv, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    rng = np.random.default_rng(args.seed)
    drugs = sorted({r["drug"] for r in rows})
    fit_rows = []
    target_rows = []

    for drug in drugs:
        dr = [r for r in rows if r["drug"] == drug]
        reps = sorted({int(r["replicate"]) for r in dr})
        concentrations = sorted({float(r["concentration_nM"]) for r in dr})
        lookup = {(int(r["replicate"]), float(r["concentration_nM"])): float(r["suppression_fraction"]) for r in dr}
        boot_pars = []
        boot_targets = defaultdict(list)

        for _ in range(args.bootstrap):
            sampled = rng.choice(reps, size=len(reps), replace=True)
            x = np.asarray(concentrations, dtype=float)
            y = np.asarray([
                np.mean([lookup[(int(rep), c)] for rep in sampled]) for c in concentrations
            ], dtype=float)
            try:
                pars = fit_curve(x, y)
            except Exception:
                continue
            boot_pars.append(pars)
            for target in (0.25, 0.50, 0.75, 0.90):
                boot_targets[target].append(target_concentration(pars, target))

        arr = np.asarray(boot_pars)
        fit_rows.append({
            "drug": drug,
            "bootstrap_n": len(arr),
            "EC50_median_nM": float(np.median(arr[:, 2])),
            "EC50_CI2.5_nM": q(arr[:, 2], 2.5),
            "EC50_CI97.5_nM": q(arr[:, 2], 97.5),
            "Hill_median": float(np.median(arr[:, 3])),
            "Emax_median": float(np.median(arr[:, 1])),
        })
        for target in (0.25, 0.50, 0.75, 0.90):
            vals = boot_targets[target]
            finite = [v for v in vals if np.isfinite(v)]
            target_rows.append({
                "drug": drug,
                "suppression_target": target,
                "n_boot": len(finite),
                "concentration_median_nM": float(np.median(finite)) if finite else float("nan"),
                "CI2.5_nM": q(finite, 2.5),
                "CI97.5_nM": q(finite, 97.5),
            })

    for path, records in ((args.fits_output, fit_rows), (args.targets_output, target_rows)):
        with Path(path).open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=records[0].keys()); w.writeheader(); w.writerows(records)
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
