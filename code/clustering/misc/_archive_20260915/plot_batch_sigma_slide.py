#!/usr/bin/env python
"""
Domain-recovery ARI vs. batch_sigma, at each modality's own Phase-1 winning
BANKSY config. Auto-discovers every bs<value> directory under
data/figure_3/cellbin_batch_sigma_slide/<mod>/ plus the two original
banksy_batch_compare endpoints for cell (bs=0, bs=1.5, at the SAME
lambda=0.5/k_geom=200 config, never rebuilt separately). Shows whether
recovery degrades smoothly or falls off a cliff as the per-slice batch
effect strengthens, and where each modality's canonical (tuned) batch_sigma
sits relative to that.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[2]
SLIDE_ROOT = SIM_PAPER_DIR / "data" / "figure_3" / "cellbin_batch_sigma_slide"
COMPARE_ROOT = SIM_PAPER_DIR / "data" / "figure_3" / "banksy_batch_compare"

MODALITIES = {
    "cell": {"label": "Cell (λ=0.5, k_geom=200)", "color": "C0", "canonical": 1.5},
    "bin": {"label": "Bin 16µm (λ=0.5, k_geom=100)", "color": "C1", "canonical": 0.7},
    "spot": {"label": "Spot (λ=0.1, k_geom=8) — stale dispersion", "color": "C2", "canonical": 0.3},
    # Re-run 2026-09-14 with spot's ACTUAL calibrated dispersion (log_mu=-2.0,
    # theta=0.25, jitter=0.10) instead of generate_strongmix.sh's stale
    # defaults (log_mu=-2.5, theta=2.0, jitter=1.0) -- see
    # spot_batch_sigma_slide_corrected.sh. Directory name differs from the
    # ari_summary_<modality>.json filename (still "spot"), hence true_mod.
    "spot_corrected": {"label": "Spot (λ=0.1, k_geom=8) — corrected dispersion, strongmix", "color": "C3", "canonical": 0.3, "true_mod": "spot"},
    # Same corrected dispersion as spot_corrected but WITHOUT --strong-domain-mix
    # -- the real domain_type_mix spot_celltype_panel.sh's actual data uses.
    "spot_weakmix": {"label": "Spot (λ=0.1, k_geom=8) — corrected dispersion, real domain mix", "color": "C4", "canonical": 0.3, "true_mod": "spot"},
}

# cell's bs=0/1.5 endpoints were only ever built at the original
# banksy_batch_compare paths (k_geom=200, unchanged) -- not duplicated under
# cellbin_batch_sigma_slide/cell/.
COMPARE_ENDPOINTS = {
    "cell": {0.0: "prebatch", 1.5: "tuned"},
    "spot": {0.0: "prebatch", 0.3: "tuned"},
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

    for bs, batch_label in COMPARE_ENDPOINTS.get(mod, {}).items():
        if bs in points:
            continue
        run_dir = COMPARE_ROOT / mod / batch_label
        ari_path = run_dir / "ari" / f"ari_summary_{mod}.json"
        score_path = COMPARE_ROOT / "scores" / f"composition_recovery_{mod}_{batch_label}.json"
        if ari_path.is_file() and score_path.is_file():
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
            if leak > 0.3:
                ax.annotate("†", (bs, ari), textcoords="offset points", xytext=(0, 6), ha="center", color="firebrick", fontsize=11)
            if canon:
                ax.axvline(bs, color=cfg["color"], linestyle=":", alpha=0.35)

    ax.set_xlabel("batch_sigma")
    ax.set_ylabel("Domain ARI (BANKSY + Harmony)")
    ax.set_title("Domain recovery vs. batch effect magnitude")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylim(-0.05, 1.0)
    ax.legend(loc="upper right")
    fig.text(0.5, -0.02, "† = slice_id leakage ARI > 0.3 | dotted vline = each modality's canonical (tuned) batch_sigma",
              ha="center", fontsize=7.5, color="dimgray")
    fig.tight_layout()

    out_path = SLIDE_ROOT / "batch_sigma_slide_domain_ari.png"
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
