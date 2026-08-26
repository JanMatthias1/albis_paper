#!/bin/bash
#SBATCH --job-name=replot_celltype_panels
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/replot_celltype_panels_%A_%a.out
#SBATCH --array=0-4
#SBATCH --time=04:00:00
#SBATCH --mem=320G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# One-off: regenerate cell-type-panel plots after the PANEL_FIGSIZE/PANEL_RECT
# aspect-ratio fix (pca_harmony.py's before/after-Harmony plot and
# clustering_leiden_louvain.py's true-vs-predicted plot now share one fixed
# canvas size) and the earlier title/label formatting fix. --plots-only in
# both scripts re-renders from the already-computed h5ad without rerunning
# PCA/Harmony/clustering.
#
# Array indices 0-1: pca_harmony.py --plots-only for bin/cell (spot already
# done interactively). Array indices 2-4: clustering_leiden_louvain.py
# --plots-only for the celltype_matched Leiden output, spot/bin/cell.
#
# NOTE on --modality below for bin/cell (idx 3, 4): bin_celltype_panel.sh /
# cell_celltype_panel.sh never pass --modality to clustering_leiden_louvain.py,
# so it silently defaulted to "spot" -- data is correct (loaded from the
# right --input), but the OUTPUT FILENAME is wrong
# (simulation_spot_z_leiden_pca_res*.h5ad inside the bin/cell folders). Using
# --modality spot here matches that existing on-disk filename; the panel
# script bug itself is a separate follow-up, not fixed by this one-off.

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

ROOT="sim_paper/data/figure_3/pca_harmony_single_cell"

case "${SLURM_ARRAY_TASK_ID:-0}" in
  0)
    MOD=bin
    H5AD="${ROOT}/${MOD}/pca_harmony_qc/simulation_${MOD}_z_pca_harmony_qc.h5ad"
    "${PYTHON_BIN}" sim_paper/code/clustering/pca_harmony.py \
        --modality "${MOD}" --input "${H5AD}" --output "${H5AD}" --plots-only
    ;;
  1)
    MOD=cell
    H5AD="${ROOT}/${MOD}/pca_harmony_qc/simulation_${MOD}_z_pca_harmony_qc.h5ad"
    "${PYTHON_BIN}" sim_paper/code/clustering/pca_harmony.py \
        --modality "${MOD}" --input "${H5AD}" --output "${H5AD}" --plots-only
    ;;
  2)
    "${PYTHON_BIN}" sim_paper/code/clustering/clustering_leiden_louvain.py \
        --modality spot --algorithm leiden --pipeline pca_harmony --resolution 0.52390625 \
        --output-dir "${ROOT}/spot/leiden_pca_qc_celltype_matched" --plots-only
    ;;
  3)
    "${PYTHON_BIN}" sim_paper/code/clustering/clustering_leiden_louvain.py \
        --modality spot --algorithm leiden --pipeline pca_harmony --resolution 0.15015625 \
        --output-dir "${ROOT}/bin/leiden_pca_qc_celltype_matched" --plots-only
    ;;
  4)
    "${PYTHON_BIN}" sim_paper/code/clustering/clustering_leiden_louvain.py \
        --modality spot --algorithm leiden --pipeline pca_harmony --resolution 0.38375 \
        --output-dir "${ROOT}/cell/leiden_pca_qc_celltype_matched" --plots-only
    ;;
esac
