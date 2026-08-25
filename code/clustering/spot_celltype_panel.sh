#!/bin/bash
#SBATCH --job-name=spot_celltype_panel
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/spot_celltype_panel_%j.out
#SBATCH --time=08:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Figure 3, cell-type recovery panel, SPOT modality (Visium-like, 100um
# spacing / 27.5um capture radius).
# Fully self-contained: generates the sim data if not already present,
# QC-filters, runs plain PCA->Harmony->Leiden, reports ARI against
# cell_type_true, and produces the ground-truth-vs-predicted UMAP +
# contingency heatmap plots.
#
# Uses the WEAK (manuscript-baseline) domain_type_mix -- separate dataset
# from the domain panel (spot_domain_panel.sh), see cell_celltype_panel.sh
# for why the two panels don't share one dataset.
#
# batch_sigma=0.15 (NOT the manuscript default 0.22): spot's large capture
# radius pools far more cells per observation than bin/cell, diluting real
# biological signal relative to the same fixed-size batch shift -- at 0.22,
# Harmony completely fails to correct it (stays as 10 isolated per-slice
# islands). 0.15 gives a real pre-Harmony separation that Harmony fully
# resolves. See figure.md, 2026-08-24/25, spot batch_sigma diagnosis.
# Pipeline: plain PCA+Harmony (this panel; domain panel uses BANKSY instead,
# see spot_domain_panel.sh -- spot is the one modality where the two panels'
# best pipelines actually differ).

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

SIM_TAG="packing_pf0p04_bsigma015"
MODALITY="spot"
SIM_QC="sim_paper/data/noisy/${SIM_TAG}/simulation_${MODALITY}_z_qc.h5ad"
# Final-config output lives under figure_3/ (2026-08-25 reorg, same convention
# as the figure_2/ move) -- not the generic clustering_<tag>/ sweep location.
CLUSTER_ROOT="sim_paper/data/figure_3/pca_harmony_single_cell/${MODALITY}"

if [[ ! -f "sim_paper/data/noisy/${SIM_TAG}/simulation_${MODALITY}_z.h5ad" ]]; then
    echo "[generate] ${SIM_TAG} not found, generating"
    "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality "${MODALITY}" --sphere-r-um 2050 --batch-sigma 0.15 --out-tag "${SIM_TAG}"
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

echo "[leiden] fixed-resolution qualitative plots (UMAP true-vs-predicted, contingency heatmap)"
"${PYTHON_BIN}" sim_paper/code/clustering/clustering_leiden_louvain.py \
    --input "${CLUSTER_ROOT}/pca_harmony_qc/simulation_${MODALITY}_z_pca_harmony_qc.h5ad" \
    --output-dir "${CLUSTER_ROOT}/leiden_pca_qc_res0p5" \
    --pipeline pca_harmony --resolution 0.5

echo "[done] cell_type_true ARI -> ${CLUSTER_ROOT}/ari_recovery_qc/ari_summary_${MODALITY}.json"
