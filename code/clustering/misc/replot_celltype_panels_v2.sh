#!/bin/bash
#SBATCH --job-name=replot_celltype_v2
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/replot_celltype_v2_%A_%a.out
#SBATCH --array=0-5
#SBATCH --time=04:00:00
#SBATCH --mem=320G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Second replot pass: after switching plot_umap_before_after (01.2_pca_harmony.py) and
# plot_umap_true_vs_predicted (step03_cluster_and_plot.py) from tight_layout/
# bbox_inches="tight" to fixed PANEL_MARGINS via subplots_adjust, so the two panel
# types share pixel-identical axes-box geometry, not just the same canvas size --
# plus pretty_label formatting (domain_true -> Domain True, cluster_label ->
# Cluster Label) now also applied in step03_cluster_and_plot.py. --plots-only
# re-renders from the already-computed h5ad without rerunning PCA/Harmony/clustering.
# Filenames for bin/cell's celltype_matched dirs now use the corrected --modality
# (see the earlier --modality panel-script fix and file rename).

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

ROOT="sim_paper/data/figure_3/pca_harmony_single_cell"

case "${SLURM_ARRAY_TASK_ID:-0}" in
  0)
    MOD=bin
    H5AD="${ROOT}/${MOD}/pca_harmony_qc/simulation_${MOD}_z_pca_harmony_qc.h5ad"
    "${PYTHON_BIN}" sim_paper/code/clustering/01.2_pca_harmony.py \
        --modality "${MOD}" --input "${H5AD}" --output "${H5AD}" --plots-only
    ;;
  1)
    MOD=cell
    H5AD="${ROOT}/${MOD}/pca_harmony_qc/simulation_${MOD}_z_pca_harmony_qc.h5ad"
    "${PYTHON_BIN}" sim_paper/code/clustering/01.2_pca_harmony.py \
        --modality "${MOD}" --input "${H5AD}" --output "${H5AD}" --plots-only
    ;;
  2)
    "${PYTHON_BIN}" sim_paper/code/clustering/step03_cluster_and_plot.py \
        --modality spot --algorithm leiden --pipeline pca_harmony --resolution 0.52390625 \
        --output-dir "${ROOT}/spot/leiden_pca_qc_celltype_matched" --plots-only
    ;;
  3)
    "${PYTHON_BIN}" sim_paper/code/clustering/step03_cluster_and_plot.py \
        --modality bin --algorithm leiden --pipeline pca_harmony --resolution 0.15015625 \
        --output-dir "${ROOT}/bin/leiden_pca_qc_celltype_matched" --plots-only
    ;;
  4)
    "${PYTHON_BIN}" sim_paper/code/clustering/step03_cluster_and_plot.py \
        --modality cell --algorithm leiden --pipeline pca_harmony --resolution 0.38375 \
        --output-dir "${ROOT}/cell/leiden_pca_qc_celltype_matched" --plots-only
    ;;
  5)
    "${PYTHON_BIN}" sim_paper/code/clustering/step03_cluster_and_plot.py \
        --modality cell --algorithm leiden --pipeline pca_harmony --resolution 0.5 \
        --output-dir "${ROOT}/cell/leiden_pca_qc_res0p5" --plots-only
    ;;
esac
