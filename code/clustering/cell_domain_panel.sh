#!/bin/bash
#SBATCH --job-name=cell_domain_panel
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/cell_domain_panel_%j.out
#SBATCH --time=08:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Figure 3, domain recovery panel, CELL modality.
# Fully self-contained: generates the sim data if not already present,
# QC-filters, builds the BANKSY matrix (lambda=0.8, domain-segmentation
# mode), runs Harmony+Leiden, reports ARI against domain_true, and produces
# the ground-truth-vs-predicted UMAP + contingency heatmap plots.
#
# Uses the STRONG domain_type_mix (every domain gets two real ~30%-enriched
# signature cell types instead of the manuscript baseline's near-uniform
# rows) -- built specifically to test whether domain structure is
# recoverable AT ALL when the compositional signal itself is genuinely
# strong. See figure.md, 2026-08-25, "domain_type_mix compositional
# distinguishability".
#
# Caveat: this dataset also uses log_mu=-2.5 (vs. cell_celltype_panel.sh's
# default 0.7) -- inherited from earlier Figure-2 expression tuning, not
# deliberately matched to the cell-type panel's config.
#
# Pipeline: BANKSY -- decided 2026-08-25 to commit to BANKSY for all three
# modalities' domain panels. Default settings collapsed to exactly 0.0/0.0
# here, but that's being treated as a tuning gap: BANKSY's neighbor-
# averaging is the only mechanism that gives a single cell any local
# compositional context at all (a lone cell's own type carries zero
# information about its neighborhood mixture on its own), so it's
# specifically plausible BANKSY matters most here, not least. Tuning sweep
# (reciprocal/max_m=0/k_geom=30, figure.md 2026-08-25) in progress -- update
# the --k-geom/--nbr-weight-decay/--max-m flags on the build step below once
# a variant is confirmed to beat the others.

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

SIM_TAG="log_mu_-2.5_theta_0.25_strongdomainmix"
MODALITY="cell"
SIM_QC="sim_paper/data/noisy/${SIM_TAG}/simulation_${MODALITY}_z_qc.h5ad"
CLUSTER_ROOT="sim_paper/data/clustering_${SIM_TAG}/${MODALITY}"

if [[ ! -f "sim_paper/data/noisy/${SIM_TAG}/simulation_${MODALITY}_z.h5ad" ]]; then
    echo "[generate] ${SIM_TAG} not found, generating"
    "${TUTORIAL_PYTHON}" sim_paper/code/data/generate_simulation_noisy.py \
        --modality "${MODALITY}" \
        --base-gene-lognormal -2.5 0.7 --theta 0.25 \
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

echo "[leiden] fixed-resolution qualitative plots (UMAP true-vs-predicted, contingency heatmap)"
"${TUTORIAL_PYTHON}" sim_paper/code/clustering/clustering_leiden_louvain.py \
    --input "${CLUSTER_ROOT}/banksy_pca_harmony_qc/simulation_${MODALITY}_z_banksy_pca_harmony_qc.h5ad" \
    --output-dir "${CLUSTER_ROOT}/leiden_pca_banksy_res0p5" \
    --pipeline pca_harmony --resolution 0.5

echo "[done] domain_true ARI -> ${CLUSTER_ROOT}/banksy_ari_recovery/ari_summary_${MODALITY}.json"
