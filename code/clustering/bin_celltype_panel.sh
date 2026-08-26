#!/bin/bash
#SBATCH --job-name=bin_celltype_panel
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/bin_celltype_panel_%j.out
#SBATCH --time=12:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Figure 3, cell-type recovery panel, BIN modality (8um Visium-HD-like grid).
# Fully self-contained: generates the sim data if not already present,
# QC-filters, runs plain PCA->Harmony->Leiden, reports ARI against
# cell_type_true, and produces the ground-truth-vs-predicted UMAP +
# contingency heatmap plots.
#
# Uses the WEAK (manuscript-baseline) domain_type_mix -- separate dataset
# from the domain panel (bin_domain_panel.sh) by design, see
# cell_celltype_panel.sh for why the two panels don't share one dataset.
#
# batch_sigma=0.5 (NOT the manuscript default 0.22): at 0.22, bin's small
# 8um aggregation footprint retains so much real signal that batch effect
# never visibly separates slices even before Harmony runs -- no meaningful
# before/after correction story to demonstrate. 0.5 gives a striking
# "flower petal" separation pre-Harmony and a near-complete merge after
# (0.33 was too weak, 0.8 pushed too far and left a visible residual
# cluster uncorrected). See figure.md, 2026-08-25, batch_sigma bracket.
# Pipeline: plain PCA+Harmony, NOT BANKSY -- BANKSY collapses bin's
# cell_type_true ARI (0.28 plain -> 0.004 BANKSY on the domain-panel
# dataset); every BANKSY parameter tried so far (k_geom, nbr_weight_decay,
# max_m) has failed to stabilize it for bin.

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

SIM_TAG="packing_pf0p04_bsigma05"
MODALITY="bin"
SIM_QC="sim_paper/data/noisy/${SIM_TAG}/simulation_${MODALITY}_z_qc.h5ad"
# Final-config output lives under figure_3/ (2026-08-25 reorg, same convention
# as the figure_2/ move) -- not the generic clustering_<tag>/ sweep location.
CLUSTER_ROOT="sim_paper/data/figure_3/pca_harmony_single_cell/${MODALITY}"

if [[ ! -f "sim_paper/data/noisy/${SIM_TAG}/simulation_${MODALITY}_z.h5ad" ]]; then
    echo "[generate] ${SIM_TAG} not found, generating"
    "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality "${MODALITY}" --sphere-r-um 2050 --batch-sigma 0.5 --out-tag "${SIM_TAG}"
fi

if [[ ! -f "${SIM_QC}" ]]; then
    echo "[qc] ${SIM_QC} not found, running 00_qc_filter.py"
    "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py --modality "${MODALITY}" --packing-tag "${SIM_TAG}"
fi

if [[ ! -f "${CLUSTER_ROOT}/pca_harmony_qc/simulation_${MODALITY}_z_pca_harmony_qc.h5ad" ]]; then
    echo "[pca_harmony] running"
    "${PYTHON_BIN}" sim_paper/code/clustering/pca_harmony.py \
        --modality "${MODALITY}" --input "${SIM_QC}" \
        --output "${CLUSTER_ROOT}/pca_harmony_qc/simulation_${MODALITY}_z_pca_harmony_qc.h5ad"
fi

echo "[ari] resolution-matched ARI recovery"
"${PYTHON_BIN}" sim_paper/code/clustering/ari_vs_ground_truth.py \
    --modality "${MODALITY}" \
    --packing-tag "${SIM_TAG}" \
    --input "${CLUSTER_ROOT}/pca_harmony_qc/simulation_${MODALITY}_z_pca_harmony_qc.h5ad" \
    --output-dir "${CLUSTER_ROOT}/ari_recovery_qc"

# Use the SAME resolution ari_vs_ground_truth.py's binary search already found
# for cell_type_true (achieves the true category count exactly), rather than a
# fixed guess -- otherwise this qualitative plot's predicted-cluster count can
# drift from the true count and look like a mismatch that isn't really there
# (found 2026-08-25 on spot: fixed res=0.5 landed on 7 clusters vs. 8 true
# types, while the matched res=0.524 hits 8/8).
RESOLUTION=$("${PYTHON_BIN}" -c "
import json
with open('${CLUSTER_ROOT}/ari_recovery_qc/ari_summary_${MODALITY}.json') as f:
    summary = json.load(f)
print(next(r['resolution'] for r in summary if r['ground_truth'] == 'cell_type_true'))
")
echo "[leiden] cell_type_true-matched-resolution qualitative plots (resolution=${RESOLUTION}, UMAP true-vs-predicted, contingency heatmap)"
"${PYTHON_BIN}" sim_paper/code/clustering/clustering_leiden_louvain.py \
    --input "${CLUSTER_ROOT}/pca_harmony_qc/simulation_${MODALITY}_z_pca_harmony_qc.h5ad" \
    --output-dir "${CLUSTER_ROOT}/leiden_pca_qc_celltype_matched" \
    --pipeline pca_harmony --resolution "${RESOLUTION}"

echo "[done] cell_type_true ARI -> ${CLUSTER_ROOT}/ari_recovery_qc/ari_summary_${MODALITY}.json"
