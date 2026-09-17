#!/usr/bin/env python
"""
Final (manuscript) version of plot_batch_sigma_slide.py: just the three
modalities at their winning Phase-1 BANKSY config, one line each -- Cell,
Bin 16um, and Spot at its corrected dispersion / strongmix run (the stale-
dispersion and real-domain-mix spot lines from the diagnostic version were
validation checks, not part of the final figure).
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[2]
# 2026-09-17: consolidated so EVERY point (including cell's bs0/bs1.5 and
# bin16um's bs0/bs0.7, which used to live only under the separate
# banksy_batch_compare/ and Figure 4 sim_data/ trees respectively) now has
# its raw h5ad + banksy_matrix/ari directly under
# data/figure_3/cellbin_batch_sigma_slide/<mod>/bs<value>/ -- one tree, no
# more special-cased endpoint lookup. banksy_batch_compare/ is retired
# (its cell/spot pieces were absorbed here or found stale, see
# project_figure3_banksy_domain_sweep memory 2026-09-17).
SLIDE_ROOT = SIM_PAPER_DIR / "data" / "figure_3" / "cellbin_batch_sigma_slide"

MODALITIES = {
    "cell": {"label": "Cell (λ=0.5, k_geom=200)", "color": "C0", "canonical": 1.5},
    # --modality passed to the pipeline is "bin" (only cell/bin/spot are
    # valid), even though the folder is "bin16um" -- true_mod fixes the
    # ari_summary_<true_mod>.json lookup below.
    "bin16um": {"label": "Bin 16µm (λ=0.5, k_geom=100)", "color": "C1", "canonical": 0.7, "true_mod": "bin"},
    "spot": {"label": "Spot (λ=0.1, k_geom=8)", "color": "C3", "canonical": 0.3},
}


def find_leak(score_path: Path) -> float:
    sd = json.loads(score_path.read_text())
    return sd["levels"]["domain"]["slice_id_leakage_ari"]


def discover_points(mod: str) -> list[tuple[float, float, float, bool]]:
    points: dict[float, tuple[float, float, bool]] = {}
    true_mod = MODALITIES[mod].get("true_mod", mod)

    mod_dir = SLIDE_ROOT / mod
    if mod_dir.is_dir():
        for run_dir in sorted(mod_dir.iterdir()):
            m = re.match(r"bs([0-9.]+)$", run_dir.name)
            if not m:
                continue
            bs = float(m.group(1))
            ari_path = run_dir / "ari" / f"ari_summary_{true_mod}.json"
            score_path = SLIDE_ROOT / "scores" / f"composition_recovery_{mod}_{run_dir.name}.json"
            if not ari_path.is_file() or not score_path.is_file():
                continue
            ari_data = json.loads(ari_path.read_text())
            domain = next(r for r in ari_data if r["ground_truth"] == "domain_true")
            leak = find_leak(score_path)
            points[bs] = (domain["ari"], leak, abs(bs - MODALITIES[mod]["canonical"]) < 1e-9)

    return [(bs, ari, leak, canon) for bs, (ari, leak, canon) in sorted(points.items())]


def main() -> None:
    all_points = {mod: discover_points(mod) for mod in MODALITIES}

    fig, ax = plt.subplots(figsize=(8, 5.5))
    for mod, cfg in MODALITIES.items():
        points = all_points[mod]
        if not points:
            continue
        xs = [p[0] for p in points]
        ys = [p[1] for p in points]
        ax.plot(xs, ys, marker="o", label=cfg["label"], color=cfg["color"], linewidth=2, markersize=5)
        for bs, ari, leak, canon in points:
            if canon:
                ax.axvline(bs, color=cfg["color"], linestyle=":", alpha=0.35)

    ax.set_xlabel("batch_sigma")
    ax.set_ylabel("Domain ARI (BANKSY + Harmony)")
    ax.set_title("Domain recovery vs. batch effect magnitude")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylim(-0.05, 1.0)
    ax.legend(loc="upper right")
    fig.text(0.5, -0.02, "dotted vline = each modality's canonical (tuned) batch_sigma",
              ha="center", fontsize=7.5, color="dimgray")
    fig.tight_layout()

    out_path = SLIDE_ROOT / "batch_sigma_slide_domain_ari_final.png"
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {out_path}")

    print("\nSummary:")
    for mod, points in all_points.items():
        print(f"  {mod}:")
        for bs, ari, leak, canon in points:
            tag = " (canonical)" if canon else ""
            print(f"    batch_sigma={bs:<6} ARI={ari:.4f}  leak={leak:.4f}{tag}")


if __name__ == "__main__":
    main()
