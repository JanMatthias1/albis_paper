#!/bin/bash
#SBATCH --job-name=sweep_diag
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/sweep_diag_%A_%a.out
#SBATCH --array=0-15%16
#SBATCH --time=04:00:00
#SBATCH --mem=150G
#SBATCH --cpus-per-task=8
#SBATCH --partition=shared
#
# Diagnostic UMAP/classification plots at representative points along the
# strongmix batch_sigma sweep curve (cellbin_batch_sigma_slide/), to
# illustrate "where it breaks" alongside the existing ARI-vs-batch_sigma
# summary (plot_batch_sigma_slide.py). All banksy_matrix + ari_recovery
# h5ad inputs already exist from the original sweep jobs -- this only adds
# the qualitative UMAP/contingency layer, reusing stored PCA embeddings.
#
# - spot_corrected (λ=0.1, k_geom=8, corrected dispersion, strongmix): the
#   notch's endpoints/edges -- bs=0.12 (plateau, last clean point before the
#   collapse), 0.18 (deep in the flat zero trough), 0.22 (trough, still
#   flat), 0.28 (still collapsed, one step before the jump back up), 0.29
#   (recovery). bs=0.0/0.3 already have full before/after-UMAP +
#   classification via the official banksy_batch_compare/spot/{prebatch,tuned}.
#   Both umap (before/after Harmony) and classification (true vs predicted)
#   are new here -- spot's sweep tree has neither yet, only PC1/PC2 scatter.
# - bin (λ=0.5, k_geom=100): bs=0/0.25/0.45/0.7 already have before/after
#   UMAP (job 35673475) but no classification plot anywhere for bin (bin
#   isn't in banksy_batch_compare/ any more -- archived stale kg200).
# - cell (λ=0.5, k_geom=200): bs=0.5/1.0 already have before/after UMAP
#   (job 35673475); bs=0/1.5 already have full diagnostics via the official
#   banksy_batch_compare/cell/{prebatch,tuned}. Only classification is new
#   for the two intermediate points.
set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

SLIDE_ROOT="sim_paper/data/figure_3/cellbin_batch_sigma_slide"

# task -> (step, modality-for---modality, sweep-subdir, bs-label)
STEPS=(umap umap umap umap umap classification classification classification classification classification classification classification classification classification classification classification)
MODS=(spot  spot  spot  spot  spot  spot          spot          spot          spot          spot          bin           bin           bin           bin           cell          cell)
SUBDIRS=(spot_corrected spot_corrected spot_corrected spot_corrected spot_corrected spot_corrected spot_corrected spot_corrected spot_corrected spot_corrected bin bin bin bin cell cell)
BS=(bs0.12 bs0.18 bs0.22 bs0.28 bs0.29 bs0.12 bs0.18 bs0.22 bs0.28 bs0.29 bs0 bs0.25 bs0.45 bs0.7 bs0.5 bs1.0)

i="${SLURM_ARRAY_TASK_ID:-0}"
STEP="${STEPS[$i]}"
MOD="${MODS[$i]}"
SUBDIR="${SUBDIRS[$i]}"
LABEL="${BS[$i]}"

RUN_DIR="${SLIDE_ROOT}/${SUBDIR}/${LABEL}"

if [[ "${STEP}" == "umap" ]]; then
    H5AD="${RUN_DIR}/banksy_matrix/simulation_${MOD}_z_banksy_pca_harmony_qc.h5ad"
    echo "[task ${i}] umap  modality=${MOD} subdir=${SUBDIR} label=${LABEL}"
    "${PYTHON_BIN}" sim_paper/code/clustering/misc/plot_banksy_batch_compare_umap.py \
        --modality "${MOD}" --input "${H5AD}" --label "${SUBDIR}_${LABEL}"
else
    H5AD="${RUN_DIR}/ari/simulation_${MOD}_z_ari_recovery.h5ad"
    echo "[task ${i}] classification  modality=${MOD} subdir=${SUBDIR} label=${LABEL}"
    "${PYTHON_BIN}" sim_paper/code/clustering/misc/plot_banksy_batch_compare_true_vs_pred.py \
        --modality "${MOD}" --input "${H5AD}" --label "${SUBDIR}_${LABEL}"
fi
echo "[done] task ${i}"
