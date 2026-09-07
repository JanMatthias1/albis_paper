#!/bin/bash
#SBATCH --job-name=cell_celltype_panel
#SBATCH --output=/dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs/cell_celltype_panel_%j.out
#SBATCH --time=08:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=shared
#
# Figure 3, cell-type recovery panel, CELL modality.
# Fully self-contained: generates the sim data if not already present,
# QC-filters, runs plain PCA->Harmony->Leiden, reports ARI against
# cell_type_true (the ground truth this panel targets), and produces the
# ground-truth-vs-predicted UMAP + contingency heatmap plots.
#
# Uses the WEAK (manuscript-baseline) domain_type_mix -- this is a separate
# dataset from the domain panel (cell_domain_panel.sh) on purpose: combining
# a strengthened domain_type_mix with Harmony's slice-batch-correction
# collapsed cell_type_true ARI from ~1.0 to ~0.26 (likely because slices now
# carry real, not just batch, compositional differences that Harmony can't
# distinguish from technical noise). Rather than chase that interaction,
# each panel gets the dataset built for its own question -- see figure.md,
# 2026-08-25, "two-dataset strategy".
#
# 2026-08-26: switched to the EXACT SAME dataset as Figure 2's
# figure_2/cell_vs_non_diseased_lung comparison instead of this panel's own
# previously-separate packing_pf0p04 tag. The old tag used sphere_r_um=2050
# (~4% packing, borrowed from bin/spot) with default dispersion
# (theta=2.0/jitter=1.0/log_mu=0.7); the shared tag uses cell's own baseline
# sphere_r_um=6000 default (~0.16% packing -- cell has no aggregation to need
# the bin/spot packing fix). This panel reads Figure 2's data directly
# (SIM_RAW/SIM_QC point at data/figure_2/<tag>/); the fallback below
# regenerates the identical config if it's missing.
#
# 2026-08-27: cell config retuned to log_mu=-2.3, theta=0.40, theta_jitter=0.15,
# batch_sigma=1.5 (tag log_mu_-2.3_theta_0.40_jitter0.15_bsigma15). batch_sigma
# had to go to 1.5 for a visible pre/post-Harmony demo, but at the earlier
# log_mu=-2.5/theta=0.25 that dropped cell_type_true ARI to 0.47 (from 0.54
# pre-batch) and cell's Figure 2 theta_hat to 0.056 (real ~0.16). Raising theta
# to 0.40 pre-compensates the batch-effect dispersion drop: ARI recovers to
# 0.64, Figure 2 theta_hat to ~0.083, genes/cell 30->42. Tradeoff: median
# total_counts overshoots (~327 vs real ~90). See figure.md.
# Pipeline: plain PCA+Harmony, NOT BANKSY -- BANKSY is neutral-to-destructive
# for cell_type_true recovery at cell resolution in every test run so far.

set -euo pipefail

mkdir -p /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/logs
source /dcs04/hicks/data/Jan/sim_project/sim_paper/code/clustering/_env.sh
cd /dcs04/hicks/data/Jan/sim_project

SIM_TAG="log_mu_-2.3_theta_0.40_jitter0.15_bsigma15"
MODALITY="cell"
SIM_RAW="sim_paper/data/figure_2/${SIM_TAG}/simulation_${MODALITY}_z.h5ad"
SIM_QC="sim_paper/data/figure_2/${SIM_TAG}/simulation_${MODALITY}_z_qc.h5ad"
NOISY_DIR="sim_paper/data/noisy/${SIM_TAG}"
FIG2_DIR="sim_paper/data/figure_2/${SIM_TAG}"
# Final-config output lives under figure_3/ (2026-08-25 reorg, same convention
# as the figure_2/ move) -- not the generic clustering_<tag>/ sweep location.
CLUSTER_ROOT="sim_paper/data/figure_3/pca_harmony_single_cell/${MODALITY}"

if [[ ! -f "${SIM_QC}" ]]; then
    if [[ ! -f "${SIM_RAW}" ]]; then
        echo "[generate] ${SIM_TAG} not found under data/figure_2/, generating (same config as Figure 2)"
        "${PYTHON_BIN}" sim_paper/code/data/generate_simulation_noisy.py \
            --modality "${MODALITY}" \
            --base-gene-lognormal -2.3 0.7 --theta 0.40 --theta-jitter 0.15 --batch-sigma 1.5 \
            --out-tag "${SIM_TAG}"
        echo "[qc] ${SIM_TAG}"
        "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py --modality "${MODALITY}" --packing-tag "${SIM_TAG}"
        echo "[move] ${NOISY_DIR} -> ${FIG2_DIR}"
        mkdir -p sim_paper/data/figure_2
        mv "${NOISY_DIR}" "${FIG2_DIR}"
    else
        echo "[qc] ${SIM_QC} not found but raw exists, running 00_qc_filter.py directly against figure_2/"
        "${PYTHON_BIN}" sim_paper/code/clustering/00_qc_filter.py --modality "${MODALITY}" \
            --input "${SIM_RAW}" --output "${SIM_QC}"
    fi
fi

if [[ ! -f "${CLUSTER_ROOT}/pca_harmony_qc/simulation_${MODALITY}_z_pca_harmony_qc.h5ad" ]]; then
    echo "[pca_harmony] running"
    "${PYTHON_BIN}" sim_paper/code/clustering/pca_harmony.py \
        --modality "${MODALITY}" --input "${SIM_QC}" \
        --no-umap-sample \
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
    --modality "${MODALITY}" \
    --input "${CLUSTER_ROOT}/pca_harmony_qc/simulation_${MODALITY}_z_pca_harmony_qc.h5ad" \
    --output-dir "${CLUSTER_ROOT}/leiden_pca_qc_celltype_matched" \
    --pipeline pca_harmony --resolution "${RESOLUTION}"

echo "[done] cell_type_true ARI -> ${CLUSTER_ROOT}/ari_recovery_qc/ari_summary_${MODALITY}.json"
