#!/usr/bin/env python3
"""Fit four-parameter logistic concentration-response curves to replicate CSV data."""

from __future__ import annotations

import argparse
import csv
import numpy as np
from scipy.optimize import curve_fit


def four_pl(x, bottom, top, ec50, hill):
    x = np.asarray(x, dtype=float)
    return bottom + (top - bottom) * (x ** hill) / (ec50 ** hill + x ** hill)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("csv")
    p.add_argument("--output", default="curve_fits.csv")
    args = p.parse_args()

    with open(args.csv, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    drugs = sorted({r["drug"] for r in rows})
    fits = []
    for drug in drugs:
        dr = [r for r in rows if r["drug"] == drug]
        concentrations = sorted({float(r["concentration_nM"]) for r in dr})
        x = np.array(concentrations)
        y = np.array([
            np.mean([float(r["suppression_fraction"]) for r in dr if float(r["concentration_nM"]) == c])
            for c in concentrations
        ])
        nonzero = x[x > 0]
        p0 = [0.0, min(1.0, float(np.max(y))), float(np.median(nonzero)), 2.0]
        pars, _ = curve_fit(four_pl, x, y, p0=p0, bounds=([0, .2, .001, .1], [.2, 1.05, 5000, 20]), maxfev=50000)
        pred = four_pl(x, *pars)
        r2 = 1 - np.sum((y - pred) ** 2) / np.sum((y - y.mean()) ** 2)
        fits.append(dict(drug=drug, bottom=pars[0], top_Emax=pars[1], EC50_nM=pars[2], hill=pars[3], R2=r2))

    with open(args.output, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fits[0].keys()); w.writeheader(); w.writerows(fits)
    for fit in fits:
        print(fit)


if __name__ == "__main__":
    main()
