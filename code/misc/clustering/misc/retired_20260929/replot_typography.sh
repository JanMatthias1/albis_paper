#!/usr/bin/env bash
# Re-render Figure 3's weak_domain_mix/pca_harmony plots (PCA/UMAP before-after-
# Harmony, ground-truth-vs-predicted UMAP, contingency heatmaps) with Figure
# 2's typography (bold 17pt titles / 16pt labels / 12pt ticks+legend, see
# 01.2_pca_harmony.py's and step03_cluster_and_plot.py's new rcParams block)
# -- no data/embedding recompute, --plots-only on both scripts reloads the
# existing pca_harmony_qc.h5ad / leiden res*.h5ad files directly. Requested
# 2026-09-18 alongside the same typography pass on
# ari_recovery_summary/plot_domain_vs_celltype.py.
#
# Array over the 4 modalities already in figure_3/weak_domain_mix/pca_harmony/:
# cell, bin (8um), bin16um, spot. bin/bin16um both pass --modality bin (no
# separate bin16um choice in either script -- see bin16um_celltype_panel.sh);
# only --output/--output-dir (CLUSTER_ROOT) differ.
#SBATCH --job-name=replot_typography
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/weak_domain_mix/pca_harmony/logs/replot_typography_%A_%a.out
#SBATCH --array=0-3
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/weak_domain_mix/pca_harmony/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

# folder, --modality, leiden --resolution (exact string from the existing
# leiden_pca_qc_celltype_matched/*.h5ad filename, so step03's plots-only mode
# reconstructs the identical output path)
FOLDERS=(cell/bs1.5 bin8um/bs0.7 bin16um/bs0.7 spot/bs0.3)
MODALITIES=(cell bin bin spot)
RESOLUTIONS=(0.7575 0.1866552734375 0.547265625 0.9210156249999999)

FOLDER="${FOLDERS[${SLURM_ARRAY_TASK_ID:-0}]}"
MOD="${MODALITIES[${SLURM_ARRAY_TASK_ID:-0}]}"
RES="${RESOLUTIONS[${SLURM_ARRAY_TASK_ID:-0}]}"
CLUSTER_ROOT="sim_paper/data/figure_3/weak_domain_mix/pca_harmony/${FOLDER}"

echo "[task ${SLURM_ARRAY_TASK_ID:-0}] folder=${FOLDER} modality=${MOD} resolution=${RES}"

echo "[replot] 01.2_pca_harmony.py --plots-only"
"${PYTHON_BIN}" sim_paper/code/clustering/01.2_pca_harmony.py \
    --modality "${MOD}" \
    --input "${CLUSTER_ROOT}/pca_harmony_qc/simulation_${MOD}_z_pca_harmony_qc.h5ad" \
    --output "${CLUSTER_ROOT}/pca_harmony_qc/simulation_${MOD}_z_pca_harmony_qc.h5ad" \
    --plots-only

echo "[replot] step03_cluster_and_plot.py --plots-only"
"${PYTHON_BIN}" sim_paper/code/clustering/step03_cluster_and_plot.py \
    --modality "${MOD}" \
    --input "${CLUSTER_ROOT}/pca_harmony_qc/simulation_${MOD}_z_pca_harmony_qc.h5ad" \
    --output-dir "${CLUSTER_ROOT}/leiden_pca_qc_celltype_matched" \
    --pipeline pca_harmony --resolution "${RES}" \
    --plots-only

echo "[done] task ${SLURM_ARRAY_TASK_ID:-0} -> ${CLUSTER_ROOT}"
