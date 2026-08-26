#!/bin/bash
#SBATCH --job-name=spot_domain_panel
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/spot_domain_panel_%j.out
#SBATCH --time=08:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Figure 3, domain recovery panel, SPOT modality (Visium-like, 100um
# spacing / 27.5um capture radius).
# Fully self-contained: generates the sim data if not already present,
# QC-filters, builds the BANKSY matrix (lambda=0.8, domain-segmentation
# mode), runs Harmony+Leiden, reports ARI against domain_true, and produces
# the ground-truth-vs-predicted UMAP + contingency heatmap plots.
#
# Uses the STRONG domain_type_mix + batch_sigma=0.15 (spot's final value,
# see spot_celltype_panel.sh).
#
# Pipeline: BANKSY, not plain -- this is the one modality/panel where
# BANKSY clearly wins: domain_true ARI = 0.611 (BANKSY) vs. 0.369 (plain)
# on this exact dataset, nearly double, and cell_type_true stays reasonably
# intact under BANKSY here too (0.275 vs. plain's 0.311 -- not destroyed,
# unlike bin/cell where BANKSY collapses cell-type recovery). This is the
# best domain-recovery result across every modality/pipeline/dataset
# combination tested. See figure.md, 2026-08-25, combined final test.

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs

# Two envs, matching 01_build_banksy_matrix.py's cross-env handoff design:
# the BANKSY build step needs the isolated sim-app-banksy stack; everything
# else (generate/QC/ari/leiden) runs in the normal sim-app-tutorial env.
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env_banksy.sh
BANKSY_PYTHON="${PYTHON_BIN}"
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
TUTORIAL_PYTHON="${PYTHON_BIN}"

cd /dcs04/hicks/data/Jan/sim_project

SIM_TAG="packing_pf0p04_bsigma015_strongdomainmix"
MODALITY="spot"
SIM_QC="sim_paper/data/noisy/${SIM_TAG}/simulation_${MODALITY}_z_qc.h5ad"
CLUSTER_ROOT="sim_paper/data/clustering_${SIM_TAG}/${MODALITY}"

if [[ ! -f "sim_paper/data/noisy/${SIM_TAG}/simulation_${MODALITY}_z.h5ad" ]]; then
    echo "[generate] ${SIM_TAG} not found, generating"
    "${TUTORIAL_PYTHON}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality "${MODALITY}" --sphere-r-um 2050 --batch-sigma 0.15 \
        --strong-domain-mix \
        --out-tag "${SIM_TAG}"
fi

if [[ ! -f "${SIM_QC}" ]]; then
    echo "[qc] ${SIM_QC} not found, running 00_qc_filter.py"
    "${TUTORIAL_PYTHON}" sim_paper/code/clustering/00_qc_filter.py --modality "${MODALITY}" --packing-tag "${SIM_TAG}"
fi

if [[ ! -f "${CLUSTER_ROOT}/banksy_pca_harmony_qc/simulation_${MODALITY}_z_banksy_pca_harmony_qc.h5ad" ]]; then
    echo "[banksy] building BANKSY matrix (lambda=0.8, domain-segmentation mode)"
    "${BANKSY_PYTHON}" sim_paper/code/clustering/01_build_banksy_matrix.py \
        --modality "${MODALITY}" --packing-tag "${SIM_TAG}" --input "${SIM_QC}"
fi

echo "[ari] resolution-matched ARI recovery"
"${TUTORIAL_PYTHON}" sim_paper/code/clustering/ari_vs_ground_truth.py \
    --modality "${MODALITY}" \
    --packing-tag "${SIM_TAG}" \
    --input "${CLUSTER_ROOT}/banksy_pca_harmony_qc/simulation_${MODALITY}_z_banksy_pca_harmony_qc.h5ad" \
    --output-dir "${CLUSTER_ROOT}/banksy_ari_recovery"

# Use the SAME resolution ari_vs_ground_truth.py's binary search already found
# for domain_true (achieves the true category count exactly), rather than a
# fixed guess -- same fix as the celltype panels (see cell_celltype_panel.sh),
# 2026-08-25: fixed res=0.5 gave 9 predicted clusters vs. 6 true domains here,
# while the matched res=0.057 hits 6/6 (ARI=0.611).
RESOLUTION=$("${TUTORIAL_PYTHON}" -c "
import json
with open('${CLUSTER_ROOT}/banksy_ari_recovery/ari_summary_${MODALITY}.json') as f:
    summary = json.load(f)
print(next(r['resolution'] for r in summary if r['ground_truth'] == 'domain_true'))
")
echo "[leiden] domain_true-matched-resolution qualitative plots (resolution=${RESOLUTION}, UMAP true-vs-predicted, contingency heatmap)"
"${TUTORIAL_PYTHON}" sim_paper/code/clustering/clustering_leiden_louvain.py \
    --input "${CLUSTER_ROOT}/banksy_pca_harmony_qc/simulation_${MODALITY}_z_banksy_pca_harmony_qc.h5ad" \
    --output-dir "${CLUSTER_ROOT}/leiden_pca_banksy_domain_matched" \
    --pipeline pca_harmony --resolution "${RESOLUTION}"

echo "[done] domain_true ARI -> ${CLUSTER_ROOT}/banksy_ari_recovery/ari_summary_${MODALITY}.json"
