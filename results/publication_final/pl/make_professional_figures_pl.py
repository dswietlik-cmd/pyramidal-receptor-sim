#!/usr/bin/env python3
"""Generate publication-ready figures for the final manuscript analysis.

Outputs vector SVG/PDF plus 800-dpi PNG files. TIFF export is optional because the
files are large. Figure uncertainty bands use a deterministic paired-replicate bootstrap
of the mean for visualization only; inferential tests are produced separately by
``scripts/run_stat_analysis.py``.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REFERENCE = ROOT / "results" / "reference"
DEFAULT_PUBLICATION = ROOT / "results" / "publication_final"
DEFAULT_OUTPUT = DEFAULT_PUBLICATION / "figures_professional"

DRUGS = ["perampanel", "memantine", "diazepam"]
STATES = ["healthy", "input_hyperexcitable", "nmda_excitotoxic"]
LABELS = {
    "healthy": "Fizjologiczny",
    "input_hyperexcitable": "Input-driven",
    "nmda_excitotoxic": "NMDA-driven",
    "perampanel": "Perampanel",
    "memantine": "Memantyna",
    "diazepam": "Diazepam",
}
COLORS = {
    "healthy": "#0072B2",
    "input_hyperexcitable": "#D55E00",
    "nmda_excitotoxic": "#009E73",
    "perampanel": "#0072B2",
    "memantine": "#CC79A7",
    "diazepam": "#009E73",
}
MARKERS = {"healthy": "o", "input_hyperexcitable": "s", "nmda_excitotoxic": "^"}
LINESTYLES = {"healthy": "-", "input_hyperexcitable": "--", "nmda_excitotoxic": "-."}

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.titlesize": 10,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "axes.linewidth": 0.8,
        "lines.linewidth": 1.5,
        "savefig.transparent": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    }
)


def four_pl(x: np.ndarray, bottom: float, top: float, ec50: float, hill: float) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    return bottom + (top - bottom) * (x**hill) / (ec50**hill + x**hill)


def clean_axis(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.55, alpha=0.7)
    ax.set_axisbelow(True)


def bootstrap_mean_ci(values: np.ndarray, rng: np.random.Generator, n_boot: int) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    idx = rng.integers(0, len(values), size=(n_boot, len(values)))
    means = values[idx].mean(axis=1)
    return np.percentile(means, [2.5, 97.5])


def save_figure(fig, output_dir: Path, stem: str, dpi: int, write_tiff: bool) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(output_dir / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(output_dir / f"{stem}.png", dpi=dpi, bbox_inches="tight")
    if write_tiff:
        fig.savefig(output_dir / f"{stem}.tif", dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def figure1(output_dir: Path, dpi: int, write_tiff: bool) -> None:
    fig, ax = plt.subplots(figsize=(7.3, 4.6))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    state_specs = [
        (0.03, 0.72, 0.27, 0.18, "Stan fizjologiczny", "napęd 1,00×; NMDA 1,00×", "#EAF3F8"),
        (0.365, 0.72, 0.27, 0.18, "Input-driven", "napęd pobudzający 1,65×", "#FCEFE8"),
        (0.70, 0.72, 0.27, 0.18, "NMDA-driven", "mnożnik patologii NMDA 4,10×", "#E9F6F1"),
    ]
    for x, y, w, h, title, subtitle, facecolor in state_specs:
        ax.add_patch(
            FancyBboxPatch(
                (x, y),
                w,
                h,
                boxstyle="round,pad=0.012,rounding_size=0.018",
                facecolor=facecolor,
                edgecolor="#555555",
                linewidth=0.9,
            )
        )
        ax.text(x + w / 2, y + h * 0.65, title, ha="center", va="center", fontsize=10, fontweight="bold")
        ax.text(x + w / 2, y + h * 0.33, subtitle, ha="center", va="center", fontsize=8)

    x0, y0, w0, h0 = 0.08, 0.39, 0.84, 0.17
    ax.add_patch(
        FancyBboxPatch(
            (x0, y0),
            w0,
            h0,
            boxstyle="round,pad=0.015,rounding_size=0.018",
            facecolor="#F7F7F7",
            edgecolor="#444444",
            linewidth=0.9,
        )
    )
    for j, (drug, target) in enumerate(
        [("Perampanel", "AMPA ↓"), ("Memantyna", "NMDA ↓"), ("Diazepam", "GABA-A ↑")]
    ):
        xx = x0 + (j + 0.5) * w0 / 3
        ax.text(xx, y0 + h0 * 0.66, drug, ha="center", va="center", fontweight="bold", fontsize=9)
        ax.text(xx, y0 + h0 * 0.32, target, ha="center", va="center", fontsize=9)
        if j < 2:
            ax.plot([x0 + (j + 1) * w0 / 3] * 2, [y0 + 0.03, y0 + h0 - 0.03], color="#C7C7C7", lw=0.7)
    ax.text(0.5, 0.616, "Warstwa farmakodynamiczna", ha="center", va="bottom", fontsize=8.4, fontweight="bold", bbox=dict(facecolor="white", edgecolor="none", alpha=0.96, pad=1.1), zorder=6)
    ax.text(0.5, 0.584, "stężenie → efektywna aktywność receptora", ha="center", va="bottom", fontsize=7.9, bbox=dict(facecolor="white", edgecolor="none", alpha=0.96, pad=0.8), zorder=6)

    x1, y1, w1, h1 = 0.12, 0.09, 0.76, 0.17
    ax.add_patch(
        FancyBboxPatch(
            (x1, y1),
            w1,
            h1,
            boxstyle="round,pad=0.015,rounding_size=0.018",
            facecolor="#F4F1FA",
            edgecolor="#444444",
            linewidth=0.9,
        )
    )
    ax.text(0.5, y1 + h1 * 0.70, "Integracja synaptyczna i generowanie wyładowań", ha="center", va="center", fontweight="bold", fontsize=9.2)
    ax.text(0.5, y1 + h1 * 0.40, "AMPA PSP + NMDA PSP + GABA-A PSP → Vm", ha="center", va="center", fontsize=7.9)
    ax.text(0.5, y1 + h1 * 0.20, "→ częstość wyładowań / ISI / PAN / utrzymanie aktywności fizjologicznej", ha="center", va="center", fontsize=7.6)
    for sx in [0.165, 0.5, 0.835]:
        ax.add_patch(FancyArrowPatch((sx, 0.715), (0.5, 0.555), arrowstyle="-|>", mutation_scale=11, lw=0.85, color="#666666"))
    ax.add_patch(FancyArrowPatch((0.5, 0.385), (0.5, 0.265), arrowstyle="-|>", mutation_scale=11, lw=0.9, color="#555555"))
    ax.text(0.985, 0.02, "n = 20 sparowanych realizacji; ziarna 20261200–20261219", ha="right", va="bottom", fontsize=7.5, color="#555555")
    save_figure(fig, output_dir, "Figure1_study_design_professional", dpi, write_tiff)


def drug_figures(df: pd.DataFrame, fits: pd.DataFrame, output_dir: Path, dpi: int, write_tiff: bool, n_boot: int, seed: int) -> None:
    rng = np.random.default_rng(seed)
    for index, drug in enumerate(DRUGS, start=2):
        fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.35), gridspec_kw={"wspace": 0.30})
        drug_data = df[df.drug == drug]
        for state in STATES:
            group = drug_data[drug_data.state == state]
            concentrations = sorted(group.concentration_nM.unique())
            means, lo, hi, suppression, suppression_lo, suppression_hi = [], [], [], [], [], []
            for concentration in concentrations:
                values = group[group.concentration_nM == concentration].firing_hz.to_numpy(float)
                means.append(values.mean())
                ci = bootstrap_mean_ci(values, rng, n_boot)
                lo.append(ci[0])
                hi.append(ci[1])
                s_values = group[group.concentration_nM == concentration].within_state_suppression_fraction.to_numpy(float) * 100
                suppression.append(s_values.mean())
                ci_s = bootstrap_mean_ci(s_values, rng, n_boot)
                suppression_lo.append(ci_s[0])
                suppression_hi.append(ci_s[1])

            x = np.asarray(concentrations, dtype=float)
            means = np.asarray(means)
            axes[0].plot(x, means, marker=MARKERS[state], ms=4.2, color=COLORS[state], linestyle=LINESTYLES[state], label=LABELS[state])
            axes[0].fill_between(x, np.asarray(lo), np.asarray(hi), color=COLORS[state], alpha=0.13, linewidth=0)

            suppression = np.asarray(suppression)
            axes[1].errorbar(
                x,
                suppression,
                yerr=[suppression - np.asarray(suppression_lo), np.asarray(suppression_hi) - suppression],
                marker=MARKERS[state],
                ms=4.0,
                color=COLORS[state],
                linestyle="none",
                capsize=2.3,
                elinewidth=0.8,
                label=LABELS[state],
            )
            fit = fits[(fits.drug == drug) & (fits.state == state)].iloc[0]
            nonzero = [v for v in concentrations if v > 0]
            xfine = np.linspace(max(0.001, min(nonzero) * 0.2), max(concentrations), 350)
            params = [fit["bottom"], fit["Emax_4PL"], fit["extrapolated_4PL_EC50_nM"], fit["Hill"]]
            axes[1].plot(xfine, four_pl(xfine, *params) * 100, color=COLORS[state], linestyle=LINESTYLES[state], lw=1.15, alpha=0.85)
            if bool(fit["EC50_estimable_within_tested_range"]):
                axes[1].axvline(float(fit["functional_EC50_nM"]), color=COLORS[state], lw=0.75, alpha=0.55, linestyle=":")

        axes[0].set_xlabel("Stężenie (nM)")
        axes[0].set_ylabel("Częstość wyładowań (Hz)")
        axes[1].set_xlabel("Stężenie (nM)")
        axes[1].set_ylabel("Supresja względem stanu bazowego (%)")
        axes[1].axhline(50, color="#777777", lw=0.75, linestyle=":", alpha=0.8)
        axes[1].set_ylim(-3, 105)
        clean_axis(axes[0])
        clean_axis(axes[1])
        axes[0].text(-0.16, 1.05, "A", transform=axes[0].transAxes, fontweight="bold", fontsize=11)
        axes[1].text(-0.16, 1.05, "B", transform=axes[1].transAxes, fontweight="bold", fontsize=11)
        axes[0].set_title("")
        axes[1].set_title("")
        axes[0].legend(frameon=False, loc="upper right")
        fig.subplots_adjust(left=0.09, right=0.99, bottom=0.17, top=0.89, wspace=0.31)
        save_figure(fig, output_dir, f"Figure{index}_{drug}_professional", dpi, write_tiff)


def therapeutic_window_figure(pm: pd.DataFrame, best: pd.DataFrame, output_dir: Path, dpi: int, write_tiff: bool) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.5), gridspec_kw={"wspace": 0.33})
    path_marker = {"input_hyperexcitable": "s", "nmda_excitotoxic": "^"}
    path_label = {"input_hyperexcitable": "Input-driven", "nmda_excitotoxic": "NMDA-driven"}

    for drug in DRUGS:
        for pathology in ["input_hyperexcitable", "nmda_excitotoxic"]:
            group = pm[(pm.drug == drug) & (pm.pathology_state == pathology)]
            rows = []
            for concentration in sorted(group.concentration_nM.unique()):
                if concentration == 0:
                    continue
                z = group[group.concentration_nM == concentration]
                healthy_preservation = (z.healthy_drug_hz / z.healthy_baseline_hz).mean() * 100
                pathological_suppression = z.pathological_suppression_fraction.mean() * 100
                rows.append((concentration, healthy_preservation, pathological_suppression))
            arr = np.asarray(rows, dtype=float)
            axes[0].plot(arr[:, 1], arr[:, 2], color=COLORS[drug], lw=1.0, alpha=0.75)
            axes[0].scatter(
                arr[:, 1],
                arr[:, 2],
                color=COLORS[drug],
                marker=path_marker[pathology],
                s=22,
                edgecolor="white",
                linewidth=0.4,
                label=f"{LABELS[drug]} - {path_label[pathology]}",
            )
    axes[0].axvline(80, color="#666666", lw=0.8, linestyle="--")
    axes[0].axvspan(80, 102, color="#EAF4EA", alpha=0.55, zorder=-5)
    axes[0].set_xlim(-2, 102)
    axes[0].set_ylim(-3, 105)
    axes[0].set_xlabel("Zachowanie aktywności fizjologicznej (%)")
    axes[0].set_ylabel("Supresja patologiczna (%)")
    axes[0].set_title("")
    clean_axis(axes[0])
    axes[0].legend(frameon=False, fontsize=6.6, loc="upper left")

    x = np.arange(len(DRUGS))
    for j, pathology in enumerate(["input_hyperexcitable", "nmda_excitotoxic"]):
        suppression_values, pan_values = [], []
        for drug in DRUGS:
            row = best[(best.drug == drug) & (best.pathology_state == pathology)].iloc[0]
            suppression_values.append(row.pathological_suppression_fraction * 100)
            pan_values.append(row.PAN_mean * 100)
        xpos = x + (j - 0.5) * 0.28
        axes[1].scatter(
            xpos,
            suppression_values,
            s=45,
            marker=path_marker[pathology],
            color="#222222",
            facecolors="none" if j == 0 else "#222222",
            label=f"Supresja – {path_label[pathology]}",
        )
        axes[1].scatter(
            xpos + 0.07,
            pan_values,
            s=95,
            marker="_",
            color=COLORS["memantine"] if j == 0 else COLORS["diazepam"],
            linewidths=1.7,
            label=f"PAN – {path_label[pathology]}",
        )
    axes[1].set_xticks(x, [LABELS[d] for d in DRUGS])
    axes[1].set_ylim(-2, 30)
    axes[1].set_ylabel("Efekt przy zachowaniu ≥80%\naktywności fizjologicznej (%)")
    axes[1].set_title("")
    clean_axis(axes[1])
    axes[1].legend(frameon=False, fontsize=6.5, loc="upper right")
    axes[0].text(-0.16, 1.05, "A", transform=axes[0].transAxes, fontweight="bold", fontsize=11)
    axes[1].text(-0.16, 1.05, "B", transform=axes[1].transAxes, fontweight="bold", fontsize=11)
    fig.subplots_adjust(left=0.09, right=0.99, bottom=0.18, top=0.88, wspace=0.33)
    save_figure(fig, output_dir, "Figure5_therapeutic_window_professional", dpi, write_tiff)


def sensitivity_figure(sensitivity: pd.DataFrame, output_dir: Path, dpi: int, write_tiff: bool) -> None:
    order = [
        "baseline",
        "EPSPd_-10%",
        "EPSPd_+10%",
        "IPSP_mag_-10%",
        "IPSP_mag_+10%",
        "NMDA_scale_-20%",
        "NMDA_scale_+20%",
        "CaMT_-2mV",
        "CaMT_+2mV",
        "Threshold_-2mV",
        "Threshold_+2mV",
    ]
    ylabels = [
        "Referencja",
        "AMPA EPSP -10%",
        "AMPA EPSP +10%",
        "GABA-A IPSP -10%",
        "GABA-A IPSP +10%",
        "Skala PSP NMDA -20%",
        "Skala PSP NMDA +20%",
        "Próg aktywacji NMDA -2 mV",
        "Próg aktywacji NMDA +2 mV",
        "Próg wyładowania -2 mV",
        "Próg wyładowania +2 mV",
    ]
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 4.45), sharey=True, gridspec_kw={"wspace": 0.08})
    for ax, drug in zip(axes, DRUGS):
        group = sensitivity[sensitivity.drug == drug].set_index("scenario").reindex(order)
        y = np.arange(len(order))[::-1]
        values = group.suppression_percent.to_numpy(float)
        errors = group.suppression_sd_percent.to_numpy(float)
        valid = np.isfinite(values)
        ax.errorbar(values[valid], y[valid], xerr=errors[valid], fmt="o", ms=4.2, color=COLORS[drug], ecolor="#777777", elinewidth=0.7, capsize=2)
        ax.axvline(float(group.loc["baseline", "suppression_percent"]), color="#555555", linestyle="--", lw=0.8)
        for yi, (scenario, value) in enumerate(zip(order, values)):
            yy = y[yi]
            if not np.isfinite(value):
                ax.text(3, yy, "NE", va="center", ha="left", fontsize=7.2, color="#555555")
            if scenario == "Threshold_+2mV":
                n_valid = int(group.loc[scenario, "n_valid"])
                ax.text(min(111, (value if np.isfinite(value) else 5) + 3), yy, f"n={n_valid}", va="center", fontsize=7.0)
        ax.set_xlim(0, 115)
        ax.set_title(LABELS[drug], fontweight="bold")
        ax.grid(axis="x", color="#D9D9D9", lw=0.5, alpha=0.7)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].set_yticks(np.arange(len(order))[::-1], ylabels)
    for ax in axes[1:]:
        ax.tick_params(labelleft=False)
    fig.supxlabel("Supresja częstości wyładowań (%)", y=0.035, fontsize=8.5)
    fig.subplots_adjust(left=0.27, right=0.99, bottom=0.16, top=0.91, wspace=0.08)
    save_figure(fig, output_dir, "FigureS1_sensitivity_professional", dpi, write_tiff)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-dir", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--publication-dir", type=Path, default=DEFAULT_PUBLICATION)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--dpi", type=int, default=800)
    parser.add_argument("--bootstrap-n", type=int, default=2000)
    parser.add_argument("--bootstrap-seed", type=int, default=20261002)
    parser.add_argument("--tiff", action="store_true", help="Also write large 800-dpi TIFF files.")
    args = parser.parse_args()

    df = pd.read_csv(args.reference_dir / "healthy_two_pathologies_drug_replicates.csv")
    paired = pd.read_csv(args.reference_dir / "healthy_two_pathologies_paired_metrics.csv")
    fits = pd.read_csv(args.publication_dir / "Table2_state_specific_4PL_fits.csv")
    best = pd.read_csv(args.publication_dir / "Table5_best_effect_with_healthy_preservation_ge80.csv")
    sensitivity = pd.read_csv(args.publication_dir / "TableS3_sensitivity_oat_final_ec50.csv")

    figure1(args.output_dir, args.dpi, args.tiff)
    drug_figures(df, fits, args.output_dir, args.dpi, args.tiff, args.bootstrap_n, args.bootstrap_seed)
    therapeutic_window_figure(paired, best, args.output_dir, args.dpi, args.tiff)
    sensitivity_figure(sensitivity, args.output_dir, args.dpi, args.tiff)
    print(f"wrote professional figures to {args.output_dir}")


if __name__ == "__main__":
    main()
