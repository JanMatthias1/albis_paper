#!/bin/bash
#SBATCH --job-name=replot_pca_harmony_qc
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/replot_pca_harmony_qc_%A_%a.out
#SBATCH --array=0-2
#SBATCH --time=12:00:00
#SBATCH --mem=320G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# One-off (2026-08-26): regenerate the pca_harmony.py before/after-Harmony
# UMAP plots on the FULL dataset (--no-umap-sample), not the default 50k
# subsample -- so the "Before/After Harmony" panel is computed on the same
# point set as the other Figure 3 panel it sits next to
# (leiden_pca_qc_celltype_matched/umap_true_vs_predicted_*.png, which
# 03_clustering_plots.py already computes on all observations, no
# subsampling). At 50k/716k (bin) and 50k/600k (cell) the two panels were
# each a real, independently-fit UMAP embedding, similar overall shape but
# not the same layout -- not a rendering artifact, just two different inputs.
# Spot (n_obs=5131) was already under the 50k threshold so this changes
# nothing there, included anyway for completeness/coherence.
# Time bumped 4h->12h since full-scale neighbors+UMAP (twice, before+after)
# is much slower than on a 50k subsample, especially for bin's 716k obs.
# Uses --plots-only against the already-computed figure_3/pca_harmony_single_cell
# h5ad for each modality -- does not touch X_pca_harmony or any downstream
# ari_recovery_qc/leiden_pca_qc_celltype_matched output.

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

MODALITIES=(spot bin cell)
MODALITY="${MODALITIES[${SLURM_ARRAY_TASK_ID:-0}]}"
H5AD="sim_paper/data/figure_3/pca_harmony_single_cell/${MODALITY}/pca_harmony_qc/simulation_${MODALITY}_z_pca_harmony_qc.h5ad"

echo "Modality: ${MODALITY}"
"${PYTHON_BIN}" sim_paper/code/clustering/pca_harmony.py \
    --modality "${MODALITY}" \
    --input "${H5AD}" \
    --output "${H5AD}" \
    --plots-only \
    --no-umap-sample
