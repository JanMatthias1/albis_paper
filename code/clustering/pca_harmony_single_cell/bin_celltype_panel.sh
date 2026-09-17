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
# 2026-08-26: switched to the EXACT SAME dataset as Figure 2's
# figure_2/bin_vs_breast_cancer_visium_hd comparison
# (packing_pf0p04_log_mu_0.0_bsigma08) instead of this panel's own
# previously-separate packing_pf0p04_bsigma05 tag -- the only difference
# between the two was log_mu (0.7 default here vs. Figure 2's real-data-
# tuned 0.0); batch_sigma=0.5 already matched. Per user decision, Figure 3's
# cell-typing panels now read Figure 2's data directly (SIM_RAW/SIM_QC point
# at data/figure_2/<tag>/, not a separate data/noisy/<tag>/ copy) so both
# figures are guaranteed to describe the same simulated bin dataset. If that
# file isn't there yet, the fallback below generates+QCs it with the
# identical config Figure 2 uses, then moves it into data/figure_2/ itself
# (same move-based pattern as
# code/count_distribution/bin_vs_breast_cancer_visium_hd_with_batch_tuned.sh,
# avoiding the data/noisy-vs-data/figure_2 path mismatch documented in
# figure.md 2026-08-26).
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

SIM_TAG="packing_pf0p04_log_mu_0.0_bsigma08"
MODALITY="bin"
SIM_RAW="sim_paper/data/figure_2/smaller_sphere/data/${SIM_TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="sim_paper/data/figure_2/smaller_sphere/data/${SIM_TAG}/simulation_${MODALITY}_z_qc.h5ad"
NOISY_DIR="sim_paper/data/noisy/${SIM_TAG}"
FIG2_DIR="sim_paper/data/figure_2/smaller_sphere/data/${SIM_TAG}"
# Final-config output lives under figure_3/ (2026-08-25 reorg, same convention
# as the figure_2/ move) -- not the generic clustering_<tag>/ sweep location.
CLUSTER_ROOT="sim_paper/data/figure_3/pca_harmony_single_cell/${MODALITY}"

if [[ ! -f "${SIM_QC}" ]]; then
    if [[ ! -f "${SIM_RAW}" ]]; then
        echo "[generate] ${SIM_TAG} not found under data/figure_2/, generating (same config as Figure 2)"
        "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
            --modality "${MODALITY}" --sphere-r-um 2050 \
            --base-gene-lognormal 0.0 0.7 --batch-sigma 0.8 \
            --out-tag "${SIM_TAG}"
        echo "[qc] ${SIM_TAG}"
        "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py --modality "${MODALITY}" --packing-tag "${SIM_TAG}"
        echo "[move] ${NOISY_DIR} -> ${FIG2_DIR}"
        mkdir -p sim_paper/data/figure_2/smaller_sphere/data
        mv "${NOISY_DIR}" "${FIG2_DIR}"
    else
        echo "[qc] ${SIM_QC} not found but raw exists, running 00_qc_filter.py directly against figure_2/"
        "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py --modality "${MODALITY}" \
            --input "${SIM_RAW}" --output "${SIM_QC}"
    fi
fi

if [[ ! -f "${CLUSTER_ROOT}/pca_harmony_qc/simulation_${MODALITY}_z_pca_harmony_qc.h5ad" ]]; then
    echo "[pca_harmony] running"
    "${PYTHON_BIN}" sim_paper/code/clustering/01.2_pca_harmony.py \
        --modality "${MODALITY}" --input "${SIM_QC}" \
        --no-umap-sample \
        --output "${CLUSTER_ROOT}/pca_harmony_qc/simulation_${MODALITY}_z_pca_harmony_qc.h5ad"
fi

echo "[ari] resolution-matched ARI recovery"
"${PYTHON_BIN}" sim_paper/code/clustering/02_leiden_resolution_sweep.py \
    --modality "${MODALITY}" \
    --packing-tag "${SIM_TAG}" \
    --input "${CLUSTER_ROOT}/pca_harmony_qc/simulation_${MODALITY}_z_pca_harmony_qc.h5ad" \
    --output-dir "${CLUSTER_ROOT}/ari_recovery_qc"

# Use the SAME resolution 02_leiden_resolution_sweep.py's binary search already found
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
"${PYTHON_BIN}" sim_paper/code/clustering/step03_cluster_and_plot.py \
    --modality "${MODALITY}" \
    --input "${CLUSTER_ROOT}/pca_harmony_qc/simulation_${MODALITY}_z_pca_harmony_qc.h5ad" \
    --output-dir "${CLUSTER_ROOT}/leiden_pca_qc_celltype_matched" \
    --pipeline pca_harmony --resolution "${RESOLUTION}"

echo "[done] cell_type_true ARI -> ${CLUSTER_ROOT}/ari_recovery_qc/ari_summary_${MODALITY}.json"
