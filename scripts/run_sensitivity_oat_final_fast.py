#!/usr/bin/env python3
"""Final OAT sensitivity analysis at the manuscript healthy-state functional EC50 values.

The EC50 values are read directly from the final state-specific 4PL table, avoiding
hard-coded carry-over from earlier calibration runs. The analysis uses the Numba core
from compare_two_pathologies_drugs_fast.py and the same n=20 paired seeds.
"""
from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from compare_two_pathologies_drugs_fast import (  # noqa: E402
    BASE_SEED,
    N_REPLICATES,
    gains_for,
    make_stream,
    pars,
    simulate_fast,
)

SCENARIOS = {
    "baseline": {},
    "EPSPd_-10%": {0: 0.90},
    "EPSPd_+10%": {0: 1.10},
    "IPSP_mag_-10%": {1: 0.90},
    "IPSP_mag_+10%": {1: 1.10},
    "NMDA_scale_-20%": {3: 0.80},
    "NMDA_scale_+20%": {3: 1.20},
    "CaMT_-2mV": {"add": (2, -2.0)},
    "CaMT_+2mV": {"add": (2, +2.0)},
    "Threshold_-2mV": {"add": (6, -2.0)},
    "Threshold_+2mV": {"add": (6, +2.0)},
}

LABELS = {
    "baseline": "Reference",
    "EPSPd_-10%": "AMPA EPSP -10%",
    "EPSPd_+10%": "AMPA EPSP +10%",
    "IPSP_mag_-10%": "GABA-A IPSP magnitude -10%",
    "IPSP_mag_+10%": "GABA-A IPSP magnitude +10%",
    "NMDA_scale_-20%": "NMDA PSP scale -20%",
    "NMDA_scale_+20%": "NMDA PSP scale +20%",
    "CaMT_-2mV": "NMDA activation threshold -2 mV",
    "CaMT_+2mV": "NMDA activation threshold +2 mV",
    "Threshold_-2mV": "Spike threshold -2 mV",
    "Threshold_+2mV": "Spike threshold +2 mV",
}


def load_final_healthy_ec50s(path: Path) -> dict[str, float]:
    df = pd.read_csv(path)
    h = df[(df["state"] == "healthy") & (df["EC50_estimable_within_tested_range"] == True)]  # noqa: E712
    out = {str(r.drug): float(r.functional_EC50_nM) for r in h.itertuples()}
    required = {"perampanel", "memantine", "diazepam"}
    if set(out) != required:
        raise RuntimeError(f"Expected final healthy EC50 values for {sorted(required)}, got {out}")
    return out


def modified_params(base: tuple, spec: dict) -> tuple:
    p = list(base)
    for key, value in spec.items():
        if isinstance(key, int):
            p[key] *= value
    if "add" in spec:
        idx, delta = spec["add"]
        p[idx] += delta
    return tuple(p)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--fits",
        default=str(ROOT / "results/publication_final/Table2_state_specific_4PL_fits.csv"),
    )
    ap.add_argument(
        "--output",
        default=str(ROOT / "results/publication_final/TableS3_sensitivity_oat_final_ec50.csv"),
    )
    ap.add_argument(
        "--figure",
        default=None,
        help="Optional path for the legacy diagnostic sensitivity plot. The canonical Supplementary Figure S1 is generated separately by make_professional_figures.py.",
    )
    args = ap.parse_args()

    ec50s = load_final_healthy_ec50s(Path(args.fits))
    streams = {rep: make_stream(BASE_SEED + rep, "healthy") for rep in range(N_REPLICATES)}
    base = pars()
    rows: list[dict] = []

    for drug in ("perampanel", "memantine", "diazepam"):
        concentration = ec50s[drug]
        ampa, nmda, gaba, _ = gains_for(drug, concentration)
        for scenario, spec in SCENARIOS.items():
            p = modified_params(base, spec)
            controls, treated, suppression = [], [], []
            for rep in range(N_REPLICATES):
                ex, inh = streams[rep]
                c = simulate_fast(ex, inh, *p, 1.0, 1.0, 1.0, 1.0)[0]
                d = simulate_fast(ex, inh, *p, ampa, nmda, gaba, 1.0)[0]
                controls.append(c)
                treated.append(d)
                if c > 0:
                    suppression.append(1.0 - d / c)

            rows.append({
                "drug": drug,
                "perturbation": LABELS[scenario],
                "scenario": scenario,
                "concentration_nM": concentration,
                "n_valid": len(suppression),
                "control_hz_mean": statistics.mean(controls),
                "control_hz_sd": statistics.stdev(controls) if len(controls) > 1 else 0.0,
                "drug_hz_mean": statistics.mean(treated),
                "drug_hz_sd": statistics.stdev(treated) if len(treated) > 1 else 0.0,
                "suppression_percent": 100.0 * statistics.mean(suppression) if suppression else math.nan,
                "suppression_sd_percent": 100.0 * statistics.stdev(suppression) if len(suppression) > 1 else (0.0 if suppression else math.nan),
            })

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    if args.figure:
        df = pd.DataFrame(rows)
        order = list(SCENARIOS)
        short = [
            "Ref", "AMPA\n-10%", "AMPA\n+10%", "GABA\n-10%", "GABA\n+10%",
            "NMDA scale\n-20%", "NMDA scale\n+20%", "NMDA thr.\n-2 mV", "NMDA thr.\n+2 mV",
            "Spike thr.\n-2 mV", "Spike thr.\n+2 mV",
        ]
        titles = {
            "perampanel": f"Perampanel ({ec50s['perampanel']:.2f} nM)",
            "memantine": f"Memantine ({ec50s['memantine']:.2f} nM)",
            "diazepam": f"Diazepam ({ec50s['diazepam']:.2f} nM)",
        }
        fig, axs = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
        for ax, drug in zip(axs, ("perampanel", "memantine", "diazepam")):
            g = df[df.drug == drug].set_index("scenario").loc[order]
            x = np.arange(len(order))
            y = g["suppression_percent"].to_numpy(float)
            e = g["suppression_sd_percent"].to_numpy(float)
            valid = np.isfinite(y)
            ax.errorbar(x[valid], y[valid], yerr=e[valid], marker="o", linestyle="none", capsize=3)
            ax.axhline(float(g.loc["baseline", "suppression_percent"]), linewidth=1, linestyle="--")
            ax.text(1, 5, "NE\n(control silent)", ha="center", va="bottom", fontsize=8)
            n_thr = int(g.loc["Threshold_+2mV", "n_valid"])
            ax.text(10, min(98, float(g.loc["Threshold_+2mV", "suppression_percent"])), f"n={n_thr}", ha="center", va="top", fontsize=8)
            ax.set_ylim(0, 105)
            ax.set_ylabel("Firing suppression (%)")
            ax.set_title(titles[drug])
            ax.grid(axis="y", alpha=0.2)
        axs[-1].set_xticks(np.arange(len(order)), short)
        fig.suptitle("One-at-a-time sensitivity analysis at final healthy-state functional EC50 values", fontsize=13)
        fig.tight_layout(rect=[0, 0, 1, 0.97])
        fig.savefig(args.figure, dpi=300, bbox_inches="tight")
        plt.close(fig)
    
    print(f"wrote {out}")
    if args.figure:
        print(f"wrote {args.figure}")


if __name__ == "__main__":
    main()
