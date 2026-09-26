#!/usr/bin/env python
"""
Manuscript panel: BANKSY vs. plain gene-expression embedding, each scored on
BOTH spatial-domain and cell-type recovery, across cell / bin16um / spot, on
the weak (Figure 2) and strong domain-mix data.

    rows     strong mix, weak mix
    columns  BANKSY -> PCA -> Harmony   (batch_sigma = 0)
             genes  -> PCA -> Harmony   (batch_sigma = 0)
             genes  -> PCA -> Harmony   (each modality's tuned batch_sigma:
                                         cell 1.5 / bin16um 0.7 / spot 0.3)

The first two columns are the controlled comparison (same data, only the
embedding differs). BANKSY is not shown under batch: its neighbour
aggregation keeps the per-slice batch shift and removes all cross-slice
overlap, so Harmony cannot correct it and Leiden clusters by slice (found
2026-09-24; archived notes: data/misc/banksy_cell_20260925/README.md). The third column shows the plain
pipeline doing its actual job, correcting a real batch effect.

Sources (all resolution-matched Leiden ARI from 02_leiden_resolution_sweep.py):
  BANKSY     figure_3/<mix>_domain_mix/banksy/{domain,cell_type}/<mod>/lam<λ>/ari/
             final runs at the λ chosen on the strong mix and reused for the weak
             mix; λ is read from code/clustering/<mix>_domain_mix/banksy/
             {domain,cell_type}/final_tasks.tsv (the tables that ran the jobs), so
             the figure cannot drift from the code. Choice + reasons:
             code/clustering/strong_domain_mix/banksy/README.md
  genes σ=0  figure_3/<mix>_domain_mix/pca_harmony/<mod>/bs0/ari_recovery_qc/
  genes σ    figure_3/<mix>_domain_mix/pca_harmony/<mod>/bs<σ>/ari_recovery_qc/
             (weak mix at its usual σ = the Figure 2 data)
Every bar is the mean over simulation seeds 2025/101/202 (seed 2025 folders
have no suffix, the others end in _seed<n>), with a sample-SD error bar and
one dot per seed; a bar with fewer than 3 finished seeds is labelled n=1/n=2.
Missing inputs are drawn as "pending" instead of failing.

Usage:
    conda activate /dcs04/hicks/data/Jan/sim_project/sim_paper/env/albis-tutorial
    python sim_paper/code/clustering/ari_recovery_summary/plot_banksy_vs_pca_recovery.py
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from manuscript_style import MODALITY_LOOKUP, apply_style
apply_style()

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[2]
PROJECT_DIR = SIM_PAPER_DIR.parent
FIG3_DIR = SIM_PAPER_DIR / "data" / "figure_3"
CODE_DIR = SCRIPT_DIR.parent

MODALITIES = ["cell", "bin16um", "spot"]
MODALITY_DISPLAY = {"cell": "Cell", "bin16um": "Bin (16 µm)", "spot": "Spot"}
# modality name baked into ari_summary_<...>.json (bin16um ran as --modality bin)
TRUE_MOD = {"cell": "cell", "bin16um": "bin", "spot": "spot"}
BATCH_SIGMA = {"cell": "1.5", "bin16um": "0.7", "spot": "0.3"}
TARGETS = ["domain_true", "cell_type_true"]
SEEDS = [2025, 101, 202]
# ground truth -> BANKSY final-run folder
TARGET_DIR = {"domain_true": "domain", "cell_type_true": "cell_type"}
ROWS = {"strong": ("Strong domain mix", "strong_domain_mix"), "weak": ("Weak domain mix", "weak_domain_mix")}
COLUMNS = ["banksy", "genes_bs0", "genes_bs"]
COLUMN_TITLE = {
    "banksy": "BANKSY → PCA → Harmony\nσ = 0",
    "genes_bs0": "Genes → PCA → Harmony\nσ = 0",
    "genes_bs": "Genes → PCA → Harmony\ntuned σ",
}

INK = "#000000"
INK_MUTED = "#3d3d3d"
GRIDLINE = "#DDDDDD"
BASELINE = "#c3c2b7"

# Same type scale as plot_domain_vs_celltype.py / Figure 2.
TITLE_SIZE = 16
LABEL_SIZE = 15
TICK_SIZE = 12
ANNOT_SIZE = 10


def read_ari(path: Path, ground_truth: str) -> float:
    """ARI for one ground truth, or NaN if the run hasn't finished yet."""
    if not path.exists():
        return np.nan
    data = json.loads(path.read_text())
    return next(r["ari"] for r in data if r["ground_truth"] == ground_truth)


def final_runs(mix: str, target: str) -> list[tuple[str, str, int, Path]]:
    """[(modality, λ, seed, ari_summary path)] from that mix's final_tasks.tsv."""
    tsv = CODE_DIR / mix / "banksy" / TARGET_DIR[target] / "final_tasks.tsv"
    return [(r["modality"], r["lambda"], int(r["seed"]),
             PROJECT_DIR / r["out_dir"] / "ari" / f"ari_summary_{TRUE_MOD[r['modality']]}.json")
            for r in csv.DictReader(tsv.open(), delimiter="\t")]


def genes_path(row: str, m: str, batched: bool, seed: int) -> Path:
    point = (f"bs{BATCH_SIGMA[m]}" if batched else "bs0") + ("" if seed == 2025 else f"_seed{seed}")
    return (FIG3_DIR / ROWS[row][1] / "pca_harmony" / m / point / "ari_recovery_qc"
            / f"ari_summary_{TRUE_MOD[m]}.json")


def load() -> tuple[dict, dict]:
    """({(row, column, mod, target): {seed: ari}}, {(mod, target): λ}).
    Unfinished runs are left out of the seed dict."""
    data, lambdas = {}, {}
    def put(key, seed, value):
        if not np.isnan(value):
            data.setdefault(key, {})[seed] = value
    for row, (_, mix) in ROWS.items():
        for target in TARGETS:
            for m, lam, seed, path in final_runs(mix, target):
                # each BANKSY bar comes from the run tuned for that target
                put((row, "banksy", m, target), seed, read_ari(path, target))
                if lambdas.setdefault((m, target), lam) != lam:
                    raise ValueError(f"{mix} {target} {m}: λ {lam} differs from the "
                                     f"other mix's {lambdas[(m, target)]} (both mixes must share λ)")
        for m in MODALITIES:
            for col, batched in (("genes_bs0", False), ("genes_bs", True)):
                for target in TARGETS:
                    for seed in SEEDS:
                        put((row, col, m, target), seed, read_ari(genes_path(row, m, batched, seed), target))
    return data, lambdas


def summarize(values: dict) -> tuple[float, float, int]:
    """(mean, sample SD, n) over finished seeds; NaN mean if none."""
    v = np.array(list(values.values()))
    if v.size == 0:
        return np.nan, np.nan, 0
    return float(v.mean()), float(v.std(ddof=1)) if v.size > 1 else np.nan, int(v.size)


def style_axis(ax) -> None:
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(BASELINE)
    ax.yaxis.grid(True, color=GRIDLINE, linewidth=0.9, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(axis="both", colors=INK_MUTED, labelsize=TICK_SIZE)
    ax.axhline(0, color=BASELINE, linewidth=1.0, zorder=1)


def plot_panel(ax, data: dict, lambdas: dict, row: str, col: str) -> None:
    x = np.arange(len(MODALITIES))
    width = 0.38
    colors = [MODALITY_LOOKUP[m] for m in MODALITIES]
    for i, target in enumerate(TARGETS):
        seeds = [data.get((row, col, m, target), {}) for m in MODALITIES]
        stats = [summarize(s) for s in seeds]
        means = np.array([s[0] for s in stats])
        sds = np.array([s[1] for s in stats])
        pos = x + (i - 0.5) * width
        # domain = solid, cell type = same modality colour, hatched on white
        kwargs = (dict(color=colors, edgecolor="white") if target == "domain_true"
                  else dict(color="white", edgecolor=colors, hatch="///"))
        bars = ax.bar(pos, np.nan_to_num(means), width, linewidth=1.2, zorder=2, **kwargs)
        ax.errorbar(pos, np.nan_to_num(means), yerr=np.nan_to_num(sds), fmt="none",
                    ecolor=INK, elinewidth=1.1, capsize=3, zorder=4)
        for xp, s in zip(pos, seeds):
            ax.scatter(np.full(len(s), xp), list(s.values()), s=12, color=INK_MUTED,
                       edgecolor="white", linewidth=0.5, zorder=5)
        # + 0.0 turns a rounded -0.00 into 0.00; flag bars missing seeds
        labels = ["pending" if n == 0 else f"{round(mu, 2) + 0.0:.2f}" + ("" if n == len(SEEDS) else f"\nn={n}")
                  for mu, _, n in stats]
        tops = np.nan_to_num(means) + np.nan_to_num(sds)
        for xp, top, label in zip(pos, tops, labels):
            ax.text(xp, max(top, 0) + 0.02, label, ha="center", va="bottom",
                    fontsize=ANNOT_SIZE, color=INK, fontweight="bold")
    ticks = [MODALITY_DISPLAY[m] for m in MODALITIES]
    if col == "genes_bs":
        ticks = [f"{t}\nσ={BATCH_SIGMA[m]}" for t, m in zip(ticks, MODALITIES)]
    elif col == "banksy":
        ticks = [f"{t}\nλ {lambdas[(m, 'domain_true')]} / {lambdas[(m, 'cell_type_true')]}"
                 for t, m in zip(ticks, MODALITIES)]
    ax.set_xticks(x)
    ax.set_xticklabels(ticks, fontsize=TICK_SIZE, color=INK)
    ax.set_ylim(-0.05, 1.08)
    style_axis(ax)


def main() -> None:
    data, lambdas = load()

    fig, axes = plt.subplots(len(ROWS), len(COLUMNS), figsize=(16, 9.6), sharey=True)
    for r, (row, (row_title, _)) in enumerate(ROWS.items()):
        for c, col in enumerate(COLUMNS):
            ax = axes[r, c]
            plot_panel(ax, data, lambdas, row, col)
            if r == 0:
                ax.set_title(COLUMN_TITLE[col], fontsize=TITLE_SIZE, fontweight="bold",
                             color=INK, pad=14)
            if c == 0:
                ax.set_ylabel(f"{row_title}\n\nAdjusted Rand Index", fontsize=LABEL_SIZE, color=INK)
    legend = [Patch(facecolor="#808080", edgecolor="white", label="Spatial domain"),
              Patch(facecolor="white", edgecolor="#808080", hatch="///", label="Cell type")]
    fig.legend(handles=legend, loc="lower center", ncol=2, frameon=False,
               fontsize=TICK_SIZE, bbox_to_anchor=(0.5, -0.03))
    fig.text(0.01, -0.03, "Bars: mean over simulation seeds 2025/101/202 (± SD, dots = seeds). "
             "BANKSY λ shown as domain / cell type.", fontsize=ANNOT_SIZE,
             color=INK_MUTED, ha="left", va="center")
    fig.tight_layout(h_pad=2.5)

    out_path = FIG3_DIR / "ari_recovery_summary" / "banksy_vs_pca_recovery.png"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {out_path}")

    for row, (row_title, _) in ROWS.items():
        print(f"\n{row_title}:")
        print(f"  {'':12s} {'λ dom/ct':>9s}   " + "  ".join(
            f"{c + ' ' + t.split('_')[0]:>20s}" for c in COLUMNS for t in TARGETS) + "   (mean ± SD, n)")
        for m in MODALITIES:
            lam = f"{lambdas[(m, 'domain_true')]}/{lambdas[(m, 'cell_type_true')]}"
            cells = []
            for col in COLUMNS:
                for t in TARGETS:
                    mu, sd, n = summarize(data.get((row, col, m, t), {}))
                    cells.append(f"{mu:6.3f} ± {sd:5.3f} ({n})")
            print(f"  {MODALITY_DISPLAY[m]:12s} {lam:>9s}   " + "  ".join(f"{c:>20s}" for c in cells))


if __name__ == "__main__":
    main()
