#!/usr/bin/env python
"""
Generate the final Figure 2 count-distribution comparison plots for the
current best sim-vs-real config per modality, each run four ways:

  - full panel   -> mean_variance/mean_dropout (the dispersion-fit
    comparison we've tuned carefully; HVG-matching would artificially
    inflate real's apparent overdispersion here, since HVG selection
    explicitly picks the most-variable genes -- see 2026-08-24 discussion).
  - HVG-matched  -> total_counts/genes_per_cell/sparsity_summary (where the
    real confound is panel size -- 556 sim genes vs. 18,085-36,601 real --
    and matching real down to its top-N highly-variable genes, N = sim's
    gene count, gives a fair like-for-like instead of a random subsample).
  - QC-filtered  -> mean_variance/mean_dropout/total_counts using sim data
    with 00_qc_filter.py's minimal filter applied (drops fully-empty /
    near-empty observations), matching how real data is already prepared
    (10x in_tissue + SpotSweeper QC) before ANY of its plots are made.
    Deliberately NOT used for sparsity_summary/empty-row-rate -- that stat
    exists specifically to report the zero-inflation problem honestly;
    pre-filtering the exact observations it's meant to count would hide the
    thing it's measuring, not fix it. See 2026-08-24 discussion in
    figure.md.
  - QC + HVG together -> total_counts/genes_per_cell specifically (best
    version of just those two -- both corrections at once). NOT a safe
    source for mean_variance/mean_dropout (still has the HVG dispersion
    distortion) or sparsity_summary (still hides the empty-row rate) --
    those caveats are each about one side of the combination and don't
    cancel by combining them.

All modes still write every plot count_distribution.py produces; which
panel actually goes in the paper for each mode is a manual pick once these
land, per the paragraph in FIGURE2_DRAFT.md.

Edit CONFIGS below as the pending marker_foldchange sweep (cell) lands --
this is a snapshot of the current best candidates as of 2026-08-24, not a
fixed final answer.

Expected environment:
    conda activate /dcs04/hicks/data/Jan/sim_project/albis/env/albis-tutorial

Example:
    python sim_paper/code/count_distribution/generate_figure2_final.py
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_PAPER_DIR = SCRIPT_DIR.parents[1]
COUNT_DISTRIBUTION_PY = SCRIPT_DIR / "count_distribution.py"
QC_FILTER_PY = SCRIPT_DIR.parent / "clustering" / "00_qc_filter.py"
NOISY_DIR = SIM_PAPER_DIR / "data" / "noisy"
REAL_QC_DIR = SIM_PAPER_DIR / "data" / "real_data_qc"
OUTPUT_ROOT = SIM_PAPER_DIR / "data" / "count_distribution" / "figure2_final"

# Each entry: one (modality, sim config, real reference) comparison to run.
# sim_tag is the data/noisy/<sim_tag>/simulation_<modality>_z.h5ad directory name.
CONFIGS = [
    dict(modality="cell", sim_tag="log_mu_-2.5_theta_0.25", real_label="non_diseased_lung"),
    dict(modality="cell", sim_tag="log_mu_-2.5_theta_0.25", real_label="lung_cancer"),
    dict(modality="bin", sim_tag="packing_pf0p04_log_mu_-1.5", real_label="breast_cancer_visium_hd"),
    dict(modality="bin", sim_tag="packing_pf0p04_log_mu_-1.5", real_label="human_pancreas_visium_hd"),
    dict(modality="spot", sim_tag="packing_pf0p04", real_label="breast_cancer_visium"),
]


def run_one(modality: str, sim_input: Path, real_input: Path, real_label: str, output_dir: Path, match_panel_size: bool) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable, str(COUNT_DISTRIBUTION_PY),
        "--modality", modality,
        "--input", str(sim_input),
        "--compare-input", str(real_input),
        "--compare-label", real_label,
        "--output-dir", str(output_dir),
    ]
    if match_panel_size:
        cmd.append("--match-panel-size")

    print(f"\n{'='*80}\n[run] {' '.join(cmd)}\n{'='*80}")
    subprocess.run(cmd, check=True)


def ensure_qc_filtered(modality: str, sim_tag: str) -> Path:
    qc_path = NOISY_DIR / sim_tag / f"simulation_{modality}_z_qc.h5ad"
    if qc_path.is_file():
        print(f"[qc] already exists: {qc_path}")
        return qc_path

    cmd = [
        sys.executable, str(QC_FILTER_PY),
        "--modality", modality,
        "--packing-tag", sim_tag,
    ]
    print(f"\n{'='*80}\n[run] {' '.join(cmd)}\n{'='*80}")
    subprocess.run(cmd, check=True)
    return qc_path


def main() -> None:
    for cfg in CONFIGS:
        modality, sim_tag, real_label = cfg["modality"], cfg["sim_tag"], cfg["real_label"]
        sim_input = NOISY_DIR / sim_tag / f"simulation_{modality}_z.h5ad"
        real_input = REAL_QC_DIR / real_label / f"{real_label}_qc.h5ad"

        if not sim_input.is_file():
            raise SystemExit(f"Missing sim input: {sim_input}")
        if not real_input.is_file():
            raise SystemExit(f"Missing real input: {real_input}")

        combo_dir = OUTPUT_ROOT / f"{modality}_vs_{real_label}"
        run_one(modality, sim_input, real_input, real_label, combo_dir / "full_panel", match_panel_size=False)
        run_one(modality, sim_input, real_input, real_label, combo_dir / "hvg_matched", match_panel_size=True)

        qc_input = ensure_qc_filtered(modality, sim_tag)
        run_one(modality, qc_input, real_input, real_label, combo_dir / "qc_filtered", match_panel_size=False)

        # Both corrections together -- best version specifically for total_counts/
        # genes_per_cell (gets QC'd sim + panel-size-matched real at once). NOT a
        # safe source for mean_variance/mean_dropout (HVG-matching still distorts
        # real's apparent dispersion, unaffected by QC'ing sim) or sparsity_summary
        # (QC-filtering sim still hides the empty-row rate that stat exists to
        # report, unaffected by HVG-matching real) -- those two caveats are each
        # about one side of this combination and don't cancel by combining them.
        run_one(modality, qc_input, real_input, real_label, combo_dir / "qc_and_hvg_matched", match_panel_size=True)

    print(f"\n[done] All comparisons written under {OUTPUT_ROOT}")


if __name__ == "__main__":
    main()
