#!/usr/bin/env python
"""
Build the task list for batch_sigma_slide_seeds.sh -- the seed-replication
follow-up to cellbin_batch_sigma_slide.sh / batch_sigma_slide_fine.sh.

Those two scripts each swept batch_sigma at a SINGLE seed (generate_simulation_
noisy.py had no --seed flag; every point on batch_sigma_slide_domain_ari_final.png
is one random draw). This adds:

  1. Two extra seeds (101, 202) for every batch_sigma point that already
     exists on disk under data/figure_3/strong_domain_mix/batch_sigma_slide/<mod>/, for
     all three modalities -- so each existing point gets a 3-way replicate
     (original + 2 new) to show whether the cliff position/plateau is stable
     under a different random tissue draw, not an artifact of one seed.
  2. Two new spot batch_sigma points beyond its canonical (0.3): 0.4 and 0.6,
     each built at 3 seeds from scratch (no prior baseline exists for these).

Deterministic -- rerunning this script reproduces a byte-identical task file.

Usage:
    python gen_batch_sigma_slide_seed_tasks.py
    -> writes batch_sigma_slide_seed_tasks.tsv next to this script
"""

from __future__ import annotations

from pathlib import Path

# Original grid; preserve row order for already submitted array jobs.
EXISTING_BS = {
    "cell": [0, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 1.5],
    "bin": [0, 0.25, 0.30, 0.35, 0.40, 0.45, 0.7],
    "spot": [0, 0.05, 0.1, 0.12, 0.13, 0.14, 0.15, 0.18, 0.2, 0.22, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3],
}

# New batch_sigma points, beyond canonical, spot only (per 2026-09-16 decision --
# cell/bin16 were left as-is). No baseline exists yet, so these need a "" (no
# --seed flag, i.e. the library default seed=2025) run in addition to the 2
# extra seeds, to seed a matching baseline point on the plot.
NEW_BS = {
    "spot": [0.4, 0.6],
}

EXTRA_SEEDS = [101, 202]

OUT_PATH = Path(__file__).resolve().parent / "batch_sigma_slide_seed_tasks.tsv"


def build_tasks() -> list[tuple[str, float, str]]:
    tasks: list[tuple[str, float, str]] = []
    for mod, bs_list in EXISTING_BS.items():
        for bs in bs_list:
            for seed in EXTRA_SEEDS:
                tasks.append((mod, bs, str(seed)))
    for mod, bs_list in NEW_BS.items():
        for bs in bs_list:
            tasks.append((mod, bs, ""))  # baseline, no --seed flag
            for seed in EXTRA_SEEDS:
                tasks.append((mod, bs, str(seed)))
    # Append controls so existing array indices 0-67 remain unchanged.
    for mod in ("cell", "bin"):
        for seed in EXTRA_SEEDS:
            tasks.append((mod, 0.05, str(seed)))
    return tasks


def main() -> None:
    tasks = build_tasks()
    lines = ["modality\tbatch_sigma\tseed"]
    for mod, bs, seed in tasks:
        lines.append(f"{mod}\t{bs:g}\t{seed}")
    OUT_PATH.write_text("\n".join(lines) + "\n")
    print(f"[write] {OUT_PATH}  ({len(tasks)} tasks)")

    by_mod: dict[str, int] = {}
    for mod, _bs, _seed in tasks:
        by_mod[mod] = by_mod.get(mod, 0) + 1
    for mod, n in sorted(by_mod.items()):
        print(f"  {mod}: {n} tasks")


if __name__ == "__main__":
    main()
