#!/bin/bash
#SBATCH --job-name=tune_bin_cell_radius
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs/tune_bin_cell_radius_%A_%a.out
#SBATCH --time=02:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#SBATCH --array=0-2
#
# One-off follow-up to the interactive cell_r_mean=7.5-vs-14 smoke test
# (2026-08-25, see figure.md/DATA_VERSIONS.md): that test found enlarging
# cell radius (molecule-spillover radius per cell, --cell-r-mean) shrinks
# bin's empty-bin fraction (20.6%->7.0%) and low-count spike, but only
# partially -- the genes-per-bin distribution was still a monotonic decay,
# not real data's unimodal peak. This sweeps further out (20/35/50um) to see
# whether the shape actually bends toward unimodal as radius keeps rising,
# or just keeps compressing without ever forming a peak.
#
# All 3 tasks fix n_cells=50000, sphere_r_um=900, capture_window_um=975 (same
# as the baseline smoke test, so bin count/inter-cell spacing don't move --
# isolates radius as the only changed variable). --allow-cell-overlap is
# required at these radii (non-overlapping placement is geometrically
# infeasible once cell-body volume exceeds the sphere's own volume at fixed
# n_cells/sphere_r_um) -- validated elsewhere in this project as not
# changing bin-level count statistics at the original r_mean=7.5, but NOT
# yet re-validated at these larger radii, so treat that as an open caveat
# when interpreting results.
#
# Each task also computes and saves genes-per-bin summary stats + a
# comparison histogram against the r_mean=7.5 baseline and r_mean=14 result,
# so results are ready to review without needing another interactive session.

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/logs
mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/data/misc/logs

ENV_PREFIX="/dcs04/hicks/data/Jan/sim_project/sim_app/env/sim-app-tutorial"
PYTHON_BIN="${ENV_PREFIX}/bin/python"
export PYTHONUNBUFFERED=1

if ! command -v conda >/dev/null 2>&1; then
    module load conda 2>/dev/null || true
fi
if command -v conda >/dev/null 2>&1; then
    CONDA_BASE="$(conda info --base)"
    # shellcheck source=/dev/null
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate "${ENV_PREFIX}"
else
    echo "WARNING: conda was not found; using ${PYTHON_BIN} directly." >&2
fi

cd /dcs04/hicks/data/Jan/sim_project

case "${SLURM_ARRAY_TASK_ID}" in
  0) R_MEAN=20 ;;
  1) R_MEAN=35 ;;
  2) R_MEAN=50 ;;
esac
TAG="radius_smoketest_r${R_MEAN}_overlap"

echo "[generate] cell_r_mean=${R_MEAN}, sphere_r_um=900 (fixed), capture_window_um=975 (fixed), allow_cell_overlap"
"${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
    --modality bin --n-cells 50000 --sphere-r-um 900 --capture-window-um 975 \
    --cell-r-mean "${R_MEAN}" --allow-cell-overlap \
    --out-tag "${TAG}"

echo "[compare] genes-per-bin summary vs. the two earlier smoke-test configs"
"${PYTHON_BIN}" - "${R_MEAN}" "${TAG}" <<'PYEOF'
import json
import sys
import numpy as np
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

r_mean = sys.argv[1]
tag = sys.argv[2]

configs = [
    ("radius_smoketest_r7p5", "cell_r_mean=7.5 (baseline)"),
    ("radius_smoketest_r14_overlap", "cell_r_mean=14, overlap"),
    (tag, f"cell_r_mean={r_mean}, overlap"),
]

summary = {}
fig, ax = plt.subplots(figsize=(9, 6))
bins = np.logspace(0, np.log10(600), 50)
for run_tag, label in configs:
    path = f"sim_paper/data/noisy/{run_tag}/simulation_bin_z.h5ad"
    try:
        adata = sc.read_h5ad(path)
    except FileNotFoundError:
        print(f"  [skip] {path} not found yet")
        continue
    genes_per_bin = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    n_genes_panel = adata.n_vars
    zero_frac = float((genes_per_bin == 0).mean())
    nonzero = genes_per_bin[genes_per_bin > 0]
    low = ((genes_per_bin >= 1) & (genes_per_bin <= 6)).mean()
    high = (genes_per_bin >= 0.5 * n_genes_panel).mean()
    mid = 1.0 - zero_frac - low - high
    summary[run_tag] = {
        "n_bins": int(adata.n_obs),
        "empty_frac": zero_frac,
        "low_1_6_frac": float(low),
        "mid_frac": float(mid),
        "high_ge50pct_panel_frac": float(high),
        "median_nonempty": float(np.median(nonzero)) if nonzero.size else None,
    }
    print(f"  {run_tag}: empty={zero_frac:.2%} low={low:.2%} mid={mid:.2%} high={high:.2%} "
          f"median_nonempty={summary[run_tag]['median_nonempty']}")
    ax.hist(nonzero, bins=bins, alpha=0.45, label=f"{label}  (empty={zero_frac:.1%})", density=True)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Genes detected per bin (non-empty bins only)")
ax.set_ylabel("Density (log scale)")
ax.set_title("cell_r_mean sweep: genes-per-bin shape vs. baseline")
ax.legend(fontsize=8)
fig.tight_layout()

out_dir = f"sim_paper/data/noisy/{tag}"
fig.savefig(f"{out_dir}/genes_per_bin_radius_compare.png", dpi=150)
with open(f"{out_dir}/genes_per_bin_radius_compare_summary.json", "w") as f:
    json.dump(summary, f, indent=2)
print(f"  [save] {out_dir}/genes_per_bin_radius_compare.png")
print(f"  [save] {out_dir}/genes_per_bin_radius_compare_summary.json")
PYEOF

echo "[done] ${TAG}"
