#!/bin/bash
#SBATCH --job-name=bin_domain_panel
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/bin_domain_panel_%j.out
#SBATCH --time=12:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Figure 3, domain recovery panel, BIN modality (8um Visium-HD-like grid).
# Fully self-contained: generates the sim data if not already present,
# QC-filters, builds the BANKSY matrix (lambda=0.8, domain-segmentation
# mode), runs Harmony+Leiden, reports ARI against domain_true, and produces
# the ground-truth-vs-predicted UMAP + contingency heatmap plots.
#
# Uses the STRONG domain_type_mix + batch_sigma=0.5 (bin's final value, see
# bin_celltype_panel.sh).
#
# Pipeline: BANKSY -- decided 2026-08-25 to commit to BANKSY for all three
# modalities' domain panels (BANKSY is a standard method for Visium-HD-like
# bin data; the default-settings collapse seen earlier was treated as a
# tuning gap, not evidence BANKSY doesn't apply here). k_geom/lambda/
# nbr_weight_decay below are the current best candidate as tuning continues
# (reciprocal/max_m=0/k_geom=30 sweep, figure.md 2026-08-25) -- update the
# --k-geom/--nbr-weight-decay/--max-m flags on the build step below once a
# variant is confirmed to beat the others, without changing the pipeline
# choice itself.

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs

# Two envs, matching 01_build_banksy_matrix.py's cross-env handoff design:
# the BANKSY build step needs the isolated sim-app-banksy stack; everything
# else (generate/QC/ari/leiden) runs in the normal albis-tutorial env.
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env_banksy.sh
BANKSY_PYTHON="${PYTHON_BIN}"
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
TUTORIAL_PYTHON="${PYTHON_BIN}"

cd /dcs04/hicks/data/Jan/sim_project

SIM_TAG="packing_pf0p04_bsigma05_strongdomainmix"
MODALITY="bin"
SIM_QC="sim_paper/data/noisy/${SIM_TAG}/simulation_${MODALITY}_z_qc.h5ad"
CLUSTER_ROOT="sim_paper/data/clustering_${SIM_TAG}/${MODALITY}"

if [[ ! -f "sim_paper/data/noisy/${SIM_TAG}/simulation_${MODALITY}_z.h5ad" ]]; then
    echo "[generate] ${SIM_TAG} not found, generating"
    "${TUTORIAL_PYTHON}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality "${MODALITY}" --sphere-r-um 2050 --batch-sigma 0.5 \
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
"${TUTORIAL_PYTHON}" sim_paper/code/clustering/02_leiden_resolution_sweep.py \
    --modality "${MODALITY}" \
    --packing-tag "${SIM_TAG}" \
    --input "${CLUSTER_ROOT}/banksy_pca_harmony_qc/simulation_${MODALITY}_z_banksy_pca_harmony_qc.h5ad" \
    --output-dir "${CLUSTER_ROOT}/banksy_ari_recovery"

# Use the SAME resolution 02_leiden_resolution_sweep.py's binary search already found
# for domain_true (achieves the true category count exactly), rather than a
# fixed guess -- same fix as the celltype panels (see cell_celltype_panel.sh
# and spot_domain_panel.sh's 2026-08-25 note).
RESOLUTION=$("${TUTORIAL_PYTHON}" -c "
import json
with open('${CLUSTER_ROOT}/banksy_ari_recovery/ari_summary_${MODALITY}.json') as f:
    summary = json.load(f)
print(next(r['resolution'] for r in summary if r['ground_truth'] == 'domain_true'))
")
echo "[leiden] domain_true-matched-resolution qualitative plots (resolution=${RESOLUTION}, UMAP true-vs-predicted, contingency heatmap)"
"${TUTORIAL_PYTHON}" sim_paper/code/clustering/step03_cluster_and_plot.py \
    --modality "${MODALITY}" \
    --input "${CLUSTER_ROOT}/banksy_pca_harmony_qc/simulation_${MODALITY}_z_banksy_pca_harmony_qc.h5ad" \
    --output-dir "${CLUSTER_ROOT}/leiden_pca_banksy_domain_matched" \
    --pipeline pca_harmony --resolution "${RESOLUTION}"

echo "[done] domain_true ARI -> ${CLUSTER_ROOT}/banksy_ari_recovery/ari_summary_${MODALITY}.json"
