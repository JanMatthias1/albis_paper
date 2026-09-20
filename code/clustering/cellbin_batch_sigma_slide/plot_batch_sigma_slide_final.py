#!/usr/bin/env python
"""
Final (manuscript) version of plot_batch_sigma_slide.py: just the three
modalities at their winning Phase-1 BANKSY config, one line each -- Cell,
Bin 16um, and Spot at its corrected dispersion / strongmix run (the stale-
dispersion and real-domain-mix spot lines from the diagnostic version were
validation checks, not part of the final figure).
"""

from __future__ import annotations

import csv
import statistics
import json
import re
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from manuscript_style import MODALITY_LOOKUP, apply_style
apply_style()

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[2]
# Clustering results live here; raw simulations may live under data/noisy/.
SLIDE_ROOT = SIM_PAPER_DIR / "data" / "figure_3" / "cellbin_batch_sigma_slide"

# Reused from the published artifact summarizing this same figure, so the
# plot and its write-up share one color identity rather than two unrelated
# palettes for the same three series.
MODALITIES = {
    "cell": {"label": "Cell (λ=0.5, k_geom=60)", "color": MODALITY_LOOKUP["cell"], "canonical": 1.5},
    # --modality passed to the pipeline is "bin" (only cell/bin/spot are
    # valid), even though the folder is "bin16um" -- true_mod fixes the
    # ari_summary_<true_mod>.json lookup below.
    "bin16um": {"label": "Bin (16 µm) (λ=0.5, k_geom=100)", "color": MODALITY_LOOKUP["bin16um"], "canonical": 0.7, "true_mod": "bin"},
    "spot": {"label": "Spot (λ=0.1, k_geom=8)", "color": MODALITY_LOOKUP["spot"], "canonical": 0.3},
}

INK = "#000000"
INK_SOFT = "#000000"
INK_MUTED = "#000000"
GRIDLINE = "#DDDDDD"
BASELINE = "#c8cdd8"

TITLE_SIZE = 17
LABEL_SIZE = 16
TICK_SIZE = 12
LEGEND_SIZE = 12

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": LABEL_SIZE,
    "text.color": INK,
    "axes.titlesize": TITLE_SIZE,
    "axes.titleweight": "bold",
    "axes.labelsize": LABEL_SIZE,
    "axes.labelcolor": INK_SOFT,
    "axes.edgecolor": BASELINE,
    "xtick.labelsize": TICK_SIZE,
    "ytick.labelsize": TICK_SIZE,
    "xtick.color": INK_SOFT,
    "ytick.color": INK_SOFT,
    "legend.fontsize": LEGEND_SIZE,
})


def discover_points(mod: str) -> list[dict]:
    """Pool baseline (seed 2025), seed101 and seed202 by numeric batch sigma."""
    points = {}
    true_mod = MODALITIES[mod].get("true_mod", mod)
    mod_dir = SLIDE_ROOT / mod
    if not mod_dir.is_dir():
        return []
    for run_dir in sorted(mod_dir.iterdir()):
        match = re.fullmatch(r"bs([0-9.]+)(?:_seed(101|202|2025))?", run_dir.name)
        if not match:
            continue
        bs = float(match[1])
        seed = int(match[2] or 2025)
        ari_path = run_dir / "ari" / f"ari_summary_{true_mod}.json"
        if not ari_path.is_file():
            continue
        rows = json.loads(ari_path.read_text())
        domain = next(row for row in rows if row["ground_truth"] == "domain_true")
        values = points.setdefault(bs, {})
        if seed in values:
            raise ValueError(f"Duplicate seed {seed} for {mod} batch_sigma={bs}")
        values[seed] = float(domain["ari"])
    result = []
    for bs, values in sorted(points.items()):
        result.append(dict(batch_sigma=bs, mean=statistics.mean(values.values()),
                           sd=statistics.stdev(values.values()) if len(values) > 1 else None,
                           n=len(values), seeds=",".join(map(str, sorted(values))),
                           missing_seeds=",".join(map(str, sorted({2025, 101, 202} - values.keys())))))
    return result


def main() -> None:
    all_points = {mod: discover_points(mod) for mod in MODALITIES}

    fig, ax = plt.subplots(figsize=(9, 6), facecolor="white")
    ax.set_facecolor("white")
    ax.yaxis.grid(True, color=GRIDLINE, linewidth=0.9, zorder=0)
    ax.set_axisbelow(True)

    for mod, cfg in MODALITIES.items():
        points = all_points[mod]
        if not points:
            continue
        xs = [p["batch_sigma"] for p in points]
        ys = [p["mean"] for p in points]
        ax.plot(xs, ys, marker="o", label=cfg["label"], color=cfg["color"],
                 linewidth=2.4, markersize=6, markeredgecolor="white",
                 markeredgewidth=0.8, zorder=3)

        replicated = [p for p in points if p["sd"] is not None]
        if replicated:
            means = np.array([p["mean"] for p in replicated])
            sds = np.array([p["sd"] for p in replicated])
            # ARI has a practical floor near 0 (a symmetric SD doesn't know
            # that) -- clip only the drawn lower whisker so the chart never
            # implies negative ARI is a real state; the true SD is unchanged
            # in the CSV this script also writes.
            lower = np.minimum(sds, means)
            ax.errorbar([p["batch_sigma"] for p in replicated], means,
                        yerr=[lower, sds], fmt="none",
                        color=cfg["color"], capsize=3.5, linewidth=1.3,
                        alpha=0.75, zorder=2)

        for p in points:
            if p["n"] < 3:
                ax.annotate(f"n={p['n']}", (p["batch_sigma"], p["mean"]),
                            xytext=(0, 7), textcoords="offset points",
                            fontsize=8, color=INK_MUTED, ha="center")
        ax.axvline(cfg["canonical"], color=cfg["color"], linestyle=":",
                    linewidth=1.3, alpha=0.45, zorder=1)

    ax.set_xlabel("Batch effect magnitude (batch_sigma)", labelpad=10)
    ax.set_ylabel("Domain ARI", labelpad=10)
    ax.set_title("Domain recovery vs. batch effect magnitude", color=INK, pad=14)
    ax.axhline(0, color=BASELINE, linewidth=1.0, zorder=1)
    ax.set_ylim(-0.03, 0.85)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(BASELINE)
    ax.tick_params(length=0)
    ax.legend(loc="upper right", frameon=False)
    fig.tight_layout()

    out_path = SLIDE_ROOT / "batch_sigma_slide_domain_ari_final.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[save] {out_path}")

    csv_path = SLIDE_ROOT / "batch_sigma_slide_domain_ari_final.csv"
    with csv_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["modality", "batch_sigma", "mean", "sd", "n", "seeds", "missing_seeds"])
        writer.writeheader()
        for mod, points in all_points.items():
            for p in points:
                writer.writerow(dict(modality=mod, **p))
                print(f"{mod}: batch_sigma={p['batch_sigma']:g} mean={p['mean']:.4f} "
                      f"sd={p['sd']} n={p['n']} missing={p['missing_seeds'] or 'none'}")
    print(f"[save] {csv_path}")


if __name__ == "__main__":
    main()
